#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 13:10:48 2020

@author: ichaudr

"""


from Bio import pairwise2, Seq, Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Blast import NCBIWWW, NCBIXML
import datetime
import time
import json
import traceback
import csv
import random
import sys


#****INPUT JSON FILE PATH****#
INPUT_JSON = "input_VC1080_allGenes.json" 
#****************************#

#Entrez request parameters
REQUEST_LIMIT = 5
SLEEP_TIME = .5
EMAIL = None
E_API = None

#BLAST parameters
taxonomic_limit = None
e_val  = None
coverage_min = None
max_hits = 50

#Nucleotide sequence parameters:
    # if max_percent_similarity = -1, no filtering is applied
isolate_promoters = False
max_percent_similarity = .8
upstream_adj = 250
downstream_adj = 50
capture_with_n = False
get_centroid = True
min_len = 1

#Output:
    #Type:
        #nuc -> returns the DNA sequence for the BLAST hits or the promoter sequences depending on isolate_promoters
        #prot -> returns the protein accession numbers for the BLAST hits
    #ftype: fasta or csv. CSV is the default for the prot output type
output_type = "nuc"
output_ftype = "fasta"
file_name = "output"

#Input:
    #input_records holds the set of protein accession keys as strings to go through the BLAST search
input_records = None


def get_percent_matches(seq1_in, seq2_in, seq_type="dna"):
    '''
    Performs a pairwise alignment between 2 sequences and returns the percent
    similarity 
    
    % Similarity is calculated as: #matches / (#matches + #mismataches); gaps are not included
    
    Parameters
    ----------
    seq1_str : string
        Sequence 1.
    seq2_str : string
        Sequence 2.
    seq_type : string, optional
        "dna" or "prot" for defining the sequence type. The default is "dna".

    Returns
    -------
    average_percent_similar : double
        The average % gapless matches between the set of optimal alignments
        generated for the two sequences.

    '''
    
    seq1 = Seq.Seq(seq1_in)
    seq2 = Seq.Seq(seq2_in)
    
    #Alignments will hold all global optimal alignments
    if seq_type == "dna":
        
        # Perform DNA alignment with the following scores:
            # Match = +2
            # Mismatch = -1
            # Gap-opening = -1
            # Gap-extending = -.5
        
        alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -1, -.5)
    elif seq_type == "prot":
        alignments = pairwise2.align.globaldx(seq1, seq2, blosum62)
    
    #If there is more than one global optimal alignment, the % identity is 
    # calculated for each alignemnt and averaged. The average is compared to
    # the threshold. 
    #Holds the percent matches for each alignment 
    percent_matches = []
       
    for align in alignments:
        
        #The format of an alignment:
            # [0] : first sequence aligned
            # [1] : second sequenced aligned
            # [2] : alignment score
            # [3] : starting position of the alignment
            # [4] : ending position of the alignment 
        
        seq1_aligned = align[0]
        seq2_aligned = align[1]
        matches = 0
        
        #For a global alignment, the start position is always 0
        align_size = align[4]
        size_adj = 0
           
        for x in range(align_size):
            
            #size_adj is subtracted from the length of the alignment to remove
            # gapped positions. 
            if seq1_aligned == "-" or seq2_aligned == "-":
                size_adj += 1
                continue
            
            if seq1_aligned[x] == seq2_aligned[x]:
                matches += 1
                
        #The size of the alignment is adjusted for the gapped positions
        current_percent_matches = (matches/(align_size-size_adj))
        
        percent_matches.append(current_percent_matches)
        #print("The percent match: ", current_percent_matches)
    
    #The average percent match is calculated for the ensemble of alignments 
    # that was returned. 
    total_match_percentage = 0
    
    for match_percentage in percent_matches:
        total_match_percentage += match_percentage
    if len(percent_matches) != 0:
        average_percent_similar = total_match_percentage / len(percent_matches)
    else:
        average_percent_similar = 0
    
    #print("The average percent match is: ", average_percent_similar)
    
    return average_percent_similar


def are_similar(seq1_str, seq2_str, percent_ident, seq_type="dna"):    
    '''
    Performs a pairwise alignment and returns true if the two sequnces are
    <X% identical based off the ensemble of alignments.
    % Identity is calculated as: #matches / (#matches + #mismataches); gaps are not included
    
    
    Parameters
    ----------
    seq1_str : string
        Sequence 1.
    seq2_str : string
        Sequence 2.
    percent_ident : double
        The maximum similarity allowed between the Sequence 1 and Sequence 2.
    seq_type : string, optional
        "dna" or "prot" for defining the sequence type. The default is "dna".

    Returns
    -------
    bool
        Returns TRUE if the two sequences have a more similar than the
        threshold passed in.

    '''
    
    return (get_percent_matches(seq1_str, seq2_str, seq_type) >= percent_ident)


def get_cent(sequences, seq_type="dna"):
    '''
    Returns the "centroid" of a set of sequences. I.e the sequence that has 
    the highest average simalrity to the other sequences in the set.

    Parameters
    ----------
    sequences : list[strings]
        A list containing the set of strings to evaluate as strings.
    seq_type : string, optional
        "dna" or "prot" for defining the sequence type. The default is "dna".


    Returns
    -------
    string
        Sequence of the centroid sequence in the set.
    
    OR
    
    Bio.Seq.Seq
        Sequence object of the centroid of the set

    '''
    
    
    #If the set of sequences only contains 1 or 2 elements, the first element
    # is returned. 
    if len(sequences) == 0:
        raise Exception("No sequences entered, no centroid")
    elif len(sequences) <= 2:
        return sequences[0]
    
    
    #Check if the inputted records are Seq objects or strings 
    object_type = ""
    if type(sequences[0]) is SeqRecord:
        object_type = "seq"
    elif type(sequences[0]) is str:
        object_type = "str"
    else:
        print("\t|~> Incorrect object type for input sequences. Needs to be string or Bio.Seq.Seq object.")
        return None
    
    #This list holds tuples where the first element of each tuple is a sequence
    # from the set that is passed in, and the second element is the 
    # average similarity of that sequence to all the other sequences 
    # in the set passed in. 
    similarity_list = []
    
    for i in range(len(sequences)):
        
        #Cummulates the total % similarity of the i sequence to every
        # sequence in the set. This value will be used to determine the
        # average similarity of sequence i to every other sequence.
        cumul_sum = 0
        
        for j in range(len(sequences) - 1):
            
            #The % similarity is not considered for a given sequence against 
            # itself. 
            if i == j:
                continue
            
            #Calculate the percent similarity if sequences are passed in 
            # as strings
            if object_type == "str":
                cumul_sum += get_percent_matches(sequences[i], 
                                                 sequences[j + 1], 
                                                 seq_type)
            
            #Caclulate the percent similarity if the sequences are passed in 
            # as sequence objects. 
            elif object_type == "seq":
                cumul_sum += get_percent_matches(str(sequences[i].seq), 
                                                 str(sequences[j + 1].seq), 
                                                 seq_type)
        
        #A tuple containing the i sequence and the average similarity is 
        # appended to the similarity_list[].
        similarity_list.append((sequences[i], cumul_sum / (len(sequences) - 1) ))
    
    #The similarity_list is iterated through to find the sequence with the 
    # highest average similarity to the other seqeunces, and it is returned 
    # as the centroid of the set. 
    curr_max_sim = 0
    curr_cent = ""
    for t in similarity_list:
        if t[1] > curr_max_sim:
            curr_cent = t[0]
    
    return curr_cent
    

def sim_filter(input_seqs, percent_ident):
    '''
    Takes in a list of sequences and threshold identity. The list of sequences 
    is filtered so that no two sequences are more similar than the threshold. 

    Parameters
    ----------
    input_seqs : list[string] or list[Bio.Seq.Seq]
        List of input sequences to be filtered as strings or sequence objects.
    percent_ident : double
        The maximum similarity allowed between any two sequences in the
        filtered list.

    Returns
    -------
    filtered_list : list[string]
        List of sequences with no two sequences that are more similar than the
        threshold.

    '''
    
    print("|~> Filtering the sequences")
    
    #Check if the inputted records are Seq objects or strings 
    object_type = ""
    if type(input_seqs[0]) is SeqRecord:
        object_type = "seq"
    elif type(input_seqs[0]) is str:
        object_type = "str"
    else:
        print("\t|~> Incorrect object type for input sequences. Needs to be string or Bio.Seq.Seq object.")
        return None
    
    #Holds 'bins'. Each bin holds all the sequences in the input set that 
    # have a higher similarity than the threshold.
    #Format of the dictionary:
        # {bin0: [seq1, seq1, seq3 ...]}
    seq_bins = {}
    filtered_list = []
    
    i = 0
    while len(input_seqs) > 0:      
        
        #Setting up the next bin and populating the first element of the list 
        # at that key with the first element remaining in the input sequence 
        # list. 
        key = "bin" + str(i)
        seq_bins.update({key : [input_seqs.pop(0)]})
        print("\t|~> Making bin " + key)
    
        #Iterate through the rest of the sequences in the input sequence list
        # to find all sequences that are within the threshold similarity and
        # add them to the bin.
        to_remove = []
        for seq in input_seqs:
            
            #Handle the comparison if the input type is a string
            if object_type == "str":
                
                if are_similar(seq_bins[key][0], seq, percent_ident):
                    seq_bins[key].append(seq)
                    to_remove.append(seq)
            
            #Handle the comparision if a Sequence object was passed in
            elif object_type == "seq":
                
                if are_similar(str(seq_bins[key][0].seq), str(seq.seq), percent_ident):
                    seq_bins[key].append(seq)
                    to_remove.append(seq)
        
        #Handles removing the sequences if the input type was strings
        if object_type == "str":
            for seq in to_remove:
                input_seqs.remove(seq)
        
        #Handles removing the sequences if the input type was a sequence object
        elif object_type == "seq":
            sequences_to_remove = []
            
            for seq in to_remove:
                sequences_to_remove.append(str(seq.seq))
            
            temp = []
            for seq in input_seqs:
                if str(seq.seq) in sequences_to_remove:
                    continue
                temp.append(seq)
            
            input_seqs = temp

        print("\t\t" + str(len(seq_bins[key])) + " sequences successfully added to " + str(key))

        i = i + 1
        
    #After sorting all of the sequences into bins based off of their
    # similarity, the centroid of each bin is added to the final filtered list. 
    
    for sbin in seq_bins.items():
        
        #Determines the centroid for each bin. This is more computationally
        # taxing, but results in the more percise center of the set.
        if get_centroid:
            print("\t|~> Getting centroid for " + str(sbin[0]))
            filtered_list.append(get_cent(sbin[1]))
        
        #Selects a random element from each bin. Not as computationally 
        # taxing but does not guarentee the representative sequence. 
        else:
            print("\t|~> Getting random sequence from " + str(sbin[0]))
            filtered_list.append(random.choice(sbin[1]))
    
    return filtered_list
    

def search_blast(input_records, max_hits=50, e_cutoff=10E-10, tax_limit=None, min_cover=None):
    '''
    Performs blast search for a set of records. 

    Parameters
    ----------
    input_records : list[string]
        The list of accession numbers to conduct the BLAST search.
    max_hits : int, optional
        Tne max number of hits to return. The default is 50.
    e_cutoff : float, optional
        The threshold for the E-value of the hits. The default is 10E-10.
    tax_limit : [string], optional
        The taxonomic limits to search in. The default is None.
    min_cover: float, optional
        The minimum coverage of the hits 


    Returns
    -------
    hits : list[(input_record, string)]
        A list containing the associated input record and accession number for the hit.

    '''
    
    print("|~> BLAST search:")
    
    #check the legnth of the input_records. 
    if len(input_records) < 1:
        raise Exception("Need at least one protein record to conduct BLAST search.")
    
    #The final list of BLAST hits
    hits = []
    
    
    #Gets the accession numbers for all the hits in the BLAST search and
    # appends hits[] with every unique record. 
    for input_record in input_records:
        
        #Fetches the protein record based off the accession number 
        print("\t|~> Getting protein record")
        
        
        for i in range(REQUEST_LIMIT):
            
            try:
            
                handle = Entrez.efetch("protein", id=input_record, rettype="fasta",
                                       retmode="text")
                
                time.sleep(SLEEP_TIME)
                break
                
            except:
                
                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                
                if i == (REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
        
        
        #Fetches the protein sequence to be used in the BLAST search
        print("\t|~> Getting protein sequence")
        input_seq = (SeqIO.read(handle, "fasta")).seq
        
        print("\t|~> Performing BLAST search: " + str(input_record))
        #Performs the appropriate BLAST search based
        if tax_limit != None:

            #Holds the parameter for the taxonomic limitation set
            taxon = ""

            if len(tax_limit) == 1:

                taxon = "txid" + str(tax_limit[0]) + "[orgn]"

            else:

                #Goes through each of the taxa limits appends to the overall entrez_querey parameter
                for i in range(len(tax_limit) - 1):
                    taxon = taxon + "txid" + str(tax_limit[0]) + "[orgn]" + " AND "

                taxon = taxon + "txid" + str(tax_limit[-1]) + "[orgn]"

            
            
            
            for i in range(REQUEST_LIMIT):
            
                try:
                
                    result_handle = NCBIWWW.qblast("blastp", "nr" ,input_seq, 
                                               entrez_query=taxon, expect=e_cutoff,
                                               hitlist_size=max_hits)

                    #Parses the resulting hits as a list
                    print("\t|~> Getting records")
                    blast_records = list(NCBIXML.parse(result_handle))
                    
                    time.sleep(SLEEP_TIME)
                    break
                    
                except:
                    
                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                    
                    if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
        
        else:
            
            
            for i in range(REQUEST_LIMIT):
            
                try:
                
                    result_handle = NCBIWWW.qblast("blastp", "nr", input_seq, 
                                           expect=e_cutoff, hitlist_size=max_hits)
                    
                    time.sleep(SLEEP_TIME)

                    #Parses the resulting hits as a list
                    print("\t|~> Getting records")
                    blast_records = list(NCBIXML.parse(result_handle))

                    break
                    
                except:
                    
                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                    
                    if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
            
            
        '''
        #Parses the resulting hits as a list
        print("\t|~> Getting records")
        blast_records = list(NCBIXML.parse(result_handle))'''
        
        #Adds each unique accession number to hits[]
        for record in blast_records[0].alignments:
            
            curr_hit_rec = record.hit_id.split('|')[-2]
            print("\t\t|~> Analyzing hit " + str(curr_hit_rec))
            
            for hit in record.hsps:
                     
                #Checks if hit meets the minimum coverage if provided 
                if min_cover:
                    
                    
                    cov = (hit.query_end - hit.query_start + 1) / (len(input_seq))
                    
                    
                    ####WRONG METHOD FOR COVERAGE####
                    #cov = 0
                    #Length of the query
                    #q_len = float(record.length)
                    
                    #Length of the alignment
                    #align_len = float(hit.align_length)
                    
                    #Calculate coverage
                    #if q_len != 0:
                    #   cov = float(align_len / q_len)
                        
                    #################################
                    
                    if(cov >= min_cover):
                        
                        #Check if the hit is already in the return list
                        if len(hits) == 0:
                            print("\t\t|~> Adding first hit (Coverage = " + str(cov) + "): " + str(curr_hit_rec))
                            hits.append((input_record, curr_hit_rec))
                        elif (not (curr_hit_rec in list(zip(*hits))[1])):
                            print("\t\t|~> Adding hit (Coverage = " + str(cov) + "): " + str(curr_hit_rec))
                            hits.append((input_record, curr_hit_rec))
                            
                    #Prints error if the minimum coverage is not met    
                    else:
                        print("\t\t|~> Hit did not meet coverage (Coverage = " + str(cov) + ") requirement: " + str(curr_hit_rec))
                else:
                    
                    #Check if the hit is already in the return list
                    if len(hits) == 0:
                        print("\t\t|~> Adding first hit: " + str(curr_hit_rec))
                        hits.append((input_record, curr_hit_rec))
                    elif (not (curr_hit_rec in list(zip(*hits))[1])):
                        print("\t\t|~> Adding hit: " + str(curr_hit_rec))
                        hits.append((input_record, curr_hit_rec))
                
            
    print("\t|~> Returning " + str(len(hits)) + " unique hits")
    return hits


def get_nuc_rec_from_prot(prot_id):
    '''
    Gets the the most complete nucelotide record for a given protein accession.

    Parameters
    ----------
    prot_id : string
        Accession number for the protein record of interest.

    Returns
    -------
    max_p_record : string
        The accession for the nucleotide record with the most complete genome.

    '''
    
    print("|~> Getting nucelotide records for " + str(prot_id))
    
    #Fetch the protein record
    print("\t|~> Fetching the protein records")
    
    
    for i in range(REQUEST_LIMIT):
            
                try:
                
                   records = Entrez.read(Entrez.efetch(db="protein",id=prot_id, rettype="ipg", 
                                       retmode="xml"))
                    
                   time.sleep(SLEEP_TIME)
                   break
                    
                except:
                    
                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                    
                    if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
    
    
   #The priority scores for the types of gene records available
    p_scores = {"NC_": 7, "AC_": 7, 
                "AE": 6, "CP": 6, "CY": 6,
                "NZ_": 5, "NT_": 5, "NW_": 5,
                "AAAA-AZZZ": 4, 
                "U": 3, "AF": 3, "AY": 3, "DQ": 3}
    
    
    #Gene records are appended to this list
    genome_list = []
    
    print("\t|~> Recording the gene records with priorities")
    if 'ProteinList' in records['IPGReport'].keys():
        for idprotein in records['IPGReport']['ProteinList']:
            if 'CDSList' in idprotein.keys():
                for cds in idprotein['CDSList']:
                    cds_acc = cds.attributes['accver']
                    cds_start = cds.attributes['start']
                    cds_stop = cds.attributes['stop']
                    cds_strand = cds.attributes['strand']
                    cds_scr = 0
                    #assign priority
                    for key in p_scores:
                        if cds_acc.startswith(key):
                            cds_scr = p_scores[key]
                    #create and append record
                    cds_rec = {'acc':cds_acc, 'start':cds_start,
                               'stop':cds_stop, 'strand':cds_strand,
                               'p_score':cds_scr}
                    genome_list.append(cds_rec)
            else:
                continue
    else:
        print("\t|~> No gene records found for " + prot_id)
        return None
    
    #Finds the genome with the max p-score
    if len(genome_list) > 0:    
        max_p_record = genome_list[0]
        for genome in genome_list:
            if genome["p_score"] > max_p_record["p_score"]:
                max_p_record = genome
        
        print("\t|~> Returning gene record for " + prot_id + ". p-score: " + 
              str(max_p_record["p_score"]))
        
        return max_p_record
    else:
        print("\t|~> No gene records found for " + prot_id)
        return None


def get_nuc_seq(nuc_rec_in, start_adj=250, stop_adj=3, isolate_promoters=False):
    '''
    Returns the nucleotide sequence of nucelotide record. 

    Parameters
    ----------
    nuc_rec_in : list[(string, string, string)]
        First element in the tuple is the input protein record the hit came from.
        The second element is the protein accession from the BLAST search of the input protein.
        The third element is the nucleotide record for the hit.
    start_adj : int, optional
        The number of nucleotides to go upstream of the start site. The default is 250.
    stop_adj : int, optional
        The number of nucelotides to go downstream of the start site. The default is 3.
    isolate_promoters : bool, optional
        Set whether to return promoter sequences or the raw sequence. The default is False.

    Returns
    -------
    SeqRecord
        Holds the promoter or nucleotide sequence.

    '''
    
    nuc_rec = nuc_rec_in[2]
    
    
    if nuc_rec is None:
        print("\t|~>Record is not valid")
        return None
    
    print("|~> Get nucleotide sequence for " + str(nuc_rec['acc']))
    
    print("\t|~> Adjusting start and stop positions")
    
    if  nuc_rec['strand']=='+':
        s_start=int(nuc_rec['start'])-start_adj
        s_stop=int(nuc_rec['start'])+stop_adj
        s_strand=1
    else:
        s_stop=int(nuc_rec['stop'])+start_adj
        s_start=int(nuc_rec['stop'])-stop_adj
        s_strand=2
    
    
    if isolate_promoters:
        
        print("\t|~> Getting genbank record")
        
        #Fetch and read the annotated GenBank record
        
        for i in range(REQUEST_LIMIT):
            
                try:
                
                   handle = Entrez.efetch(db="nuccore",id=nuc_rec['acc'], strand=s_strand,
                               seq_start=s_start, seq_stop=s_stop,
                               rettype='gbwithparts', retmode="XML")
        
    
                   genome_record = Entrez.read(handle, "xml")
                    
                   time.sleep(SLEEP_TIME)    
                   break
                    
                except:
                    
                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                    
                    if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
        
        
        
        print("\t|~> Parsing intervals for coding regions")
        # Find all coding regions in the returned GenBank sequence. 
        coding_intervals = []
        
        sequence = genome_record[0]['GBSeq_sequence']
        
        for feature in genome_record[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'gene':
                if "GBInterval_from" in feature['GBFeature_intervals'][0]:
                    coding_start = feature['GBFeature_intervals'][0]['GBInterval_from']
                    coding_end = feature['GBFeature_intervals'][0]['GBInterval_to']
                    coding_intervals.append((coding_start, coding_end))
        
        #The FASTA ID for the promoter sequence is in the following format:
            # p_NucleotideRecord
        print("\t|~> Returning promoter sequence")
        return_id = "p_" + str(nuc_rec['acc'])

        #Setting up the description for the FASTA record
        return_description = "Original_Query_Protein " + str(nuc_rec_in[0]) + " BLAST_Hit_Accession " + str(nuc_rec_in[1])
        
        #If there is only one coding region in the selected sequence, then 
        # the sequence is returned unmodified. 
        if len(coding_intervals) == 1:
            #Appends information to record description
            print("\t\t|~>No-additional-coding-regions-found-Returning-full-sequence")
            return SeqRecord(Seq.Seq(sequence), id= return_id, description=return_description)
        
        #If no coding intervals are indentified, None is returned.
        elif len(coding_intervals) == 0:
            print("\t\t|~> No coding intervals found for record: " + str(nuc_rec['acc']) + ".")
            return None
        
        #The start of the promoter is set to the start/end of the upstream gene
        # based on the directionality. ( --> --> or <-- -->)
        promoter_start = max(int(coding_intervals[-2][0]), int(coding_intervals[-2][1]))
        
        
        #Everything upstream of the promoter start is clipped off the 
        # sequence and the substring is returned.
        return_seq = str(sequence[promoter_start : ])
        #Appends information to record description
        print("\t\t|~>Successfully-clipped-off-upstream-coding-regions")
        return SeqRecord(Seq.Seq(return_seq), id= return_id, description=return_description)
        
    #If promoters aren't being isolated
    else:
        
        print("\t|~> Getting FASTA record")
        #Fetch the requested nucleotide sequence and return without any 
        # modification. 
        
        
        for i in range(REQUEST_LIMIT):
            
                try:
                
                  handle = Entrez.efetch(db="nuccore",id=nuc_rec['acc'], strand=s_strand,
                           seq_start=s_start, seq_stop=s_stop,
                           rettype='fasta', retmode="txt")
                    
                  time.sleep(SLEEP_TIME)   
                  break
                    
                except:
                    
                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")
                    
                    if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
        
        
        print("\t|~> Returnig sequence")
        return SeqIO.read(handle, "fasta")


def load_json(j_file):
    '''
    This function loads all the data from the input JSON file into global
    variables. 
    
    Parameters
    ----------
    j_file : string, optional
        The name/path of the input JSON file. The default is "input.json".

    '''
    
    #Entrez request parameters
    global REQUEST_LIMIT
    global SLEEP_TIME
    global EMAIL 
    global E_API 
    
    #BLAST parameters
    global taxonomic_limit
    global e_val
    global coverage_min
    global max_hits
    
    #Nucleotide sequence parameters:
    global isolate_promoters 
    global max_percent_similarity
    global upstream_adj
    global downstream_adj 
    global capture_with_n
    global get_centroid
    global min_len
    
    #Output:
    global output_type
    global output_ftype 
    global file_name 
    
    #Input:
    global input_records 
    
    
    try:
    
        #Open the JSON config file
        file_reader = json.load(open(j_file))
        
        #Read in the input records
        input_records = file_reader["input_records"]
        
        #Read in the output format and data type (prot or dna)
        output_ftype = file_reader["output_parameters"][0]["file_type"]
        output_type = file_reader["output_parameters"][0]["output_type"]
        file_name = "output/" + file_reader["output_parameters"][0]["file_name"] + "_" + str(datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S"))
        
        #Read in Entrez parameters
        REQUEST_LIMIT = file_reader["entrez_parameters"][0]["request_limit"]
        SLEEP_TIME = file_reader["entrez_parameters"][0]["sleep_time"]
        EMAIL = file_reader["entrez_parameters"][0]["email"]
        E_API = file_reader["entrez_parameters"][0]["api_key"]
        
        #Read in BLAST parameters
        taxonomic_limit = file_reader["blast_parameters"][0]["tax_limit"]
        e_val  = file_reader["blast_parameters"][0]["e_val"]
        coverage_min = file_reader["blast_parameters"][0]["coverage_min"]
        max_hits = file_reader["blast_parameters"][0]["max_hits"]
            
        if output_type == "prot":
            output_ftype = "csv"
        
        if output_type == "nuc":
            #Read in parameters for nucleotide sequece ID
            isolate_promoters = file_reader["nucleotide_parameters"][0]["isolate_promoters"]
            max_percent_similarity = file_reader["nucleotide_parameters"][0]["max_sim"]
            upstream_adj = file_reader["nucleotide_parameters"][0]["upstream_adj"]
            downstream_adj = file_reader["nucleotide_parameters"][0]["downstream_adj"]
            capture_with_n = file_reader["nucleotide_parameters"][0]["capture_with_n"]
            get_centroid = file_reader["nucleotide_parameters"][0]["get_centroid"]
            min_len = file_reader["nucleotide_parameters"][0]["min_length"]
            
    except:
        print(traceback.format_exc())
        

def main_run():
    '''
    This is the main function that will parse the parameters from the JSON
    file and run through the appropriate pipeline. 
    
    '''
    #Takes note of the start time
    start_time = datetime.datetime.now()
       
    #Loads parameters from the JSON
    print("|~> Loading parameters from input JSON")
    if(len(sys.argv) == 2):
        load_json("input/" + sys.argv[1])
    else:
        load_json("input/" + INPUT_JSON)
    print("\t|~> Input: return " + output_type + " to a ." + output_ftype + " file: " + file_name)
    
    #Set Entrez profile    
    Entrez.email = EMAIL
    Entrez.apikey = E_API
    
    
    #Handles the pipeline if only the protein accesions from the BLAST search
    # are requested. The results are stored in an output CSV file.
    if output_type == "prot":
        
        #Conduct the BLAST search
        blast_records = search_blast(input_records, max_hits, e_val, 
                                taxonomic_limit, coverage_min)
        
        
        #Form the file name and open the csv writer 
        o_file = file_name + ".csv"
        o_file_stream = csv.writer(open(o_file, mode='w'))
        
        #Write the column header 
        o_file_stream.writerow(["Input record", "Hit Accession #"])
        
        #Write each individual hit to the output file. 
        print("|~> Writing to output file")
        for hit in blast_records:
            o_file_stream.writerow([hit[0], hit[1]])
        
    
    #Handles the pipeline is the nucelotide sequence is requested of the 
    # BLAST hits.
    elif output_type == "nuc":
        
        #Perform the BLAST search
        blast_hits = search_blast(input_records, max_hits, e_val, 
                                  taxonomic_limit, coverage_min)
        
        #Hold the nucleotide records for each of the BLAST results
        nucleotide_records = []
        
        for hit in blast_hits:
            nucleotide_records.append((hit[0], hit[1], get_nuc_rec_from_prot(hit[1])))
            
        
        #Holds the sequence records
        nuc_seq_rec = []
        
        #For each record, fetch  the sequence record that is requested
        for record in nucleotide_records:
            
            #Get the sequence record requested
            seq = get_nuc_seq(record, upstream_adj, downstream_adj, 
                              isolate_promoters)
            
            if seq is None:
                continue
            
            if len(seq.seq) > min_len:
                if not capture_with_n and ("N"in seq.seq or "n" in seq.seq):
                    continue
                nuc_seq_rec.append(seq)
            else:
                print("|~> Record did not meet min_len requirement")
        
        
        
        #If no filtering is set, sequences are written to the output file
        if max_percent_similarity == -1:
                    
            #Open the requested output file if writing to a csv
            if output_ftype == "csv":
                o_file = file_name + ".csv"
                o_file_stream = csv.writer(open(o_file, mode='w'))
                
                #Writer the column headers
                o_file_stream.writerow(["Accession # _ Input Record", "Sequence"])
                
                print("|~> Writing to output file")
                for record in nuc_seq_rec:
                    o_file_stream.writerow([ record.id , record.seq ])
               
                
            #Handle output if a fasta file is requested
            elif output_ftype == "fasta":
                
                #The output file name
                o_file = file_name + ".fasta"
                
                #Write the records to the output file
                print("|~> Writing to output file")
                SeqIO.write(nuc_seq_rec, o_file, "fasta")
        
        
        #Handles the filtering pipeline if that is set in the parameters
        elif max_percent_similarity > 0:
            
            '''
            seq_recs = []
            
            #Extracts all sequences from the sequence records to filter
            for rec in nuc_seq_rec:
                if rec is not None:
                    seq_recs.append(rec)
            '''
            
            #Filteres the list of sequences
            filtered_seqs = sim_filter(nuc_seq_rec, max_percent_similarity)
            
            #Write those sequences to the output file
            
            #Open the requested output file if writing to a csv
            if output_ftype == "csv":
                o_file = file_name + ".csv"
                o_file_stream = csv.writer(open(o_file, mode='w'))
                
                #Writer the column headers
                o_file_stream.writerow(["Accession #", "Requested Sequence"])
                
                #Write each indiival sequence to the output file
                print("|~> Writing to output file")
               
                for seq in filtered_seqs:
                    o_file_stream.writerow([ seq.id , seq.seq])
                   
                
                
            #Handle output if a fasta file is requested
            elif output_ftype == "fasta":
                
                #The output file name
                o_file = file_name + ".fasta"
                
                '''
                #Each of the raw sequences needed to be converted to a
                # sequence record before writing to a fasta file
                filtered_seq_records = []
                
                i = 0
                for seq in filtered_seqs:
                    filtered_seq_records.append(SeqRecord(Seq.Seq(seq), id=("filtered_seq " + str(i))))
                    i = i + 1
                '''
                #Write the records to the output file
                print("|~> Writing to output file")
                SeqIO.write(filtered_seqs, o_file, "fasta")
    
    #Takes note of the end time
    end_time = datetime.datetime.now()
    
    #Calculate total time taken
    delta_time = end_time - start_time
    
    print("*"*10 + " Finished in " + str(delta_time.total_seconds() / 60) + "min" + " " + "*"*10)

main_run()


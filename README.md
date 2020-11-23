<h1> sequence_filter </h1>

**Overview**

![Image of Yaktocat](https://github.com/ErillLab/I-BINFTools/blob/master/sequence_filter/archive-misc/seq-filter-pipeline.png)

A set of protein accesions are passed in to yield BLAST results or nucleotide sequences of the BLAST results. 

 **How it works** <br><br>
A BLAST search is conducted for each of the protein accessions passed in, and all unique hits are returned by the `search_blast()` function. The unique hits can be written to a CSV file or continue through the pipline depending on the configuration in the JSON file.
<br>
<br>
The set of unique protein records from the BLAST searches are used to pull the "most complete" nucleotide record for that protein (`get_nuc_rec_from_prot()`). The set of nucelotide records are used to pull the entire nucleotide sequences or the isolated promoter sequence based on the parameters in the JSON (```get_nuc_seq()```).
<br>
<br>
The returned nucelotide sequences can be filtered to a specific threshold similarity. For example, if the threshold is set to .8, no two sequence will be more than 80% similar in the filtered list. It works by sorting all sequences within the thrshold similarity into bins, and the centroid of each bin is selected to be added to the filtered list.
<br>
<br>
The output sequences are written to the output file (CSV or FASTA).

 **Usage** <br>
 
 The parameters and input sequences are set in a JSON file following the structure below. The JSON file must be placed in the `input` folder and the filename must be specified in the `filter_seq.py` file. 
 
 <p> Example JSON file structure and definitions: </p>
 <br>

```JSON
    {
        
        "entrez_parameters": [{
        
            "request_limit" : 5,
            "sleep_time": 0.5,
            "email" : "ichaudr1@umbc.edu",
            "api_key" : null
            
        }],
    
        "blast_parameters" : [{
    
        	"tax_limit" : null,
        	"e_val" : 10E-10,
        	"coverage_min": 0.8,
        	"max_hits": 500
    
        }],
    
        "nucleotide_parameters" : [{
    
        	"isolate_promoters": true,
        	"capture_with_n": false,
        	"get_centroid": true,
        	"min_length" : 10,
        	"max_sim" : 0.75,
        	"upstream_adj" : 250,
        	"downstream_adj" : 3
    
        }],
    
        "output_parameters" : [{
    
        	"output_type" : "nuc",
        	"file_type" : "fasta",
        	"file_name" : "rvvAE-A1552_BLAST-p500"
    
        }],
    
        "input_records" : ["WP_000373826.1", "WP_000173586.1", "WP_000778624.1"]
        
            
    }
```

Parameter | Definition 
---|---
entrez_parameters | 
`request_limit` | The number of attempts that can be made for each request to the NCBI databases. 
`sleep_time` | The delay between each request in seconds
`email` | Email to be associated with the Entrez requests. Should there be issues with your requests, NCBI will contact you before blocking access. 
`api_key` | NCBI API key
---|------ 
blast_parameters |
`tax_limit` | Limits the BLAST search to a specified taxon
`e-val` | Cutoff for the e-value for the BLAST hits that are returned
`coverage_min` | Minimum coverage for the BLAST hits that are returned 
`max_hits`| Max number of hits returned.
---|------
nucleotide_parameters |
`isolate_promoters` | To return the promoter region. If `false`, will return the entire nucleotide sequence
`capture_with_n` | To return sequences with ambiguous nucloetides (N). If `true`, all sequences will be returned. 
`get_centroid` | Applys to the filtering process. If `true`, the centroid for each bin is selected. This is more computationally expensive as it requires many pairwise alignments. If `false`, a random member from each bin is selected. This is faster, but not as exact. 
`min_length` | Minimum length of the nucleotide sequence returned. 
`max_sim` | The similarity threhold to be used for the filtering process. If set to `-1`, the results will not be filtered. 
`upstream_adj` | How far upstream from the ATG start to retreive for the nucleotide sequence. 
`downstream_adj` | How far downstream from the ATG start to retreive for the nucleotide sequence.
---|------
output_parameters | 
`output_type` | `nuc` - returns nucleotide sequence <br> `prot` - returns protein record
`file_type` | CSV or FASTA file type. For the protein records, the CSV is default. 
`file_name` | Name of output file. The final output will be written in the `output` folder with the name `file_name + DATE-TIME + file_type` (i.e. rvvAE-A1552_BLAST-p500_04-07-2020_19-23-26.fasta)

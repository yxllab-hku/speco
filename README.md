# speco
Short Peptide and Enzyme Co-Occurrence analysis pipeline was used in this paper: [Expanded sequence space of radical S-adenosylmethionine-dependent enzyme involved in post-translational macrocyclization](https://www.researchsquare.com/article/rs-1789925/v1)
Speco requires prodigal-short,hmmsearch, bigscape and fasttree. If <20000 bp sequence exists, multiple non-coding nucleotide "N" will be concatenated in the 3-end accordingly in order to be processed by prodigal-short. Usage example:python3 speco.py -p PF00067 -l 20 -m 100 -t orf -u 2 -d 2 -i /data/lab/ -o /data/lab/output

Change line 34,35 in speco_fold_in.py:

```
    pfam_acc_file = '/data/bbhe/Hmmscan/pfam34.0/pfam_34.0_acc.txt'  # 
    hmm_file_path = '/data/bbhe/Hmmscan/pfam34.0/Pfam-A.hmm'  # 
```
to your own path:
```
    pfam_acc_file = 'your path/pfam_34.0_acc.txt'  # 
    hmm_file_path = 'your path/Pfam-A.hmm'  # 
```
Pfam database (pfam34.0/Pfam-A.hmm) can be downloaded from https://www.ebi.ac.uk/interpro/download/Pfam/
pfam_acc_file can be downloaded from https://github.com/Bio-bbhe/speco/tree/main/pfam_file

Below is workflow of speco. A web version is available via http://biospeco.org/, which can only handel very small sequence data since it's hosted on my private cloud server.

<img src="https://user-images.githubusercontent.com/82441159/175210002-078ddb69-27ca-45e6-875f-4bf6a3117652.png" width="500" height="500">

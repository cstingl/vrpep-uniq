# vrpep-uniq
Tools for calculating the uniqueness and other characteristics of antibody variable region peptides. 
The following steps describe the generation of tryptic variable region peptides from antibody sequence data. The tools described here are customized to process data units available from the Observed Antibody Space (<https://opig.stats.ox.ac.uk/webapps/oas/>).[^1][^2]

1. Extraction and checking of amino acid sequence and optional sampling of N entries

```
read-vrseq-vOAS.pl OAS-data-unit.csv.gz seq-check sample=1000
```

*Note:*

- `seq-check` ensures consistency between the full variable region amino acid sequence and the individual CDR and FWR region sequence.
- `sample=N` reads a sample of N entries, giving preference to entries with *Redundancy* greater than 1.
- The output will be a tab-separated file containing amino acid sequence data with the extension `*.OAS.N1000.txt`.


2. Compute tryptic peptides, region assignments, and mutations

```
extract-vrseq-vOAS.pl OAS-data-unit.OAS.N1000.txt
```

*Note:* Following output files ware written:

- OAS-data-unit.OAS.OAS.N1000.peptides.txt: table of tryptic peptides 
- OAS-data-unit.OAS.OAS.N1000-PEP.txt: antibody variable region assignment to tryptic peptides
- OAS-data-unit.OAS.OAS.N1000-MUT.txt: table of all mutation (defined as deviation to the germline sequence), linked to the region and tryptic peptide


[^1]: Olsen, T.H., Boyles, F., and Deane C.M. (2021). Protein Science.
[^2]: Kovaltsuk, A., Leem, J. et al (2018). J. Immunol.

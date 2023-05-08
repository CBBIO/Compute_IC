# Compute_IC
Compute Information Content from a list of GO IDs

## Install
+ Create a new virtual environment ```python3 -m venv venv_name```
+ Activate new enviroment ```source venv_name/bin/activate```
+ Install requirements ```pip install -r requirements.txt```

## Input Format
```
#Accession\tGOid
lcl|CADEPI010000048.1_prot_CAB3369837.1_10820	GO:0000922
lcl|CADEPI010000048.1_prot_CAB3369837.1_10820	GO:0005876
```

## Usage
+ As a CLI app
  + Not providing GPAD annotations ```python compute_IC.py path/to/input/file```
  + Providing GPAD annotations ```python compute_IC.py path/to/input/file --annotation path/to/gpad/file```
+ As a module
  + Simply ```import compute_IC``` and use the functions in your code

## Outputs
With ```python compute_IC.py...``` it generates some output folders (```outputs``` and ```plots```)
+ outputs/filename.tsv has the following format[^1]
  + Seq_ID: Sequence IDs as per input file entry
  + Category: biological_process, molecular_function or cellular_component
  + GOid: GO:XXXXXXX identifier
  + Descentant_count: Number of offpring for GO:XXXXXXX
  + Depth_level: Distance from root 
  + Information content: Information content for GO:XXXXXXX
  + Description: description for GO:XXXXXXX
```
# Seq_ID\tCategory\tGOid\tDescendant_count\tDepth_level\tInformation_content\tDescription
lcl|CADEPI010000384.1_prot_CAB3384984.1_25967	molecular_function	GO:0005525	0	8	7.173966949456064	GTP binding
lcl|CADEPI010000015.1_prot_CAB3364279.1_5262	molecular_function	GO:0005509	4	5	6.926897219781138	calcium ion binding
lcl|CADEPI010000037.1_prot_CAB3368322.1_9305	molecular_function	GO:0003899	9	7	9.69997493502881	DNA-directed 5'-3' RNA polymerase activity
lcl|CADEPI010000066.1_prot_CAB3371947.1_12930	molecular_function	GO:0003697	8	5	9.39086641174421	single-stranded DNA binding
lcl|CADEPI010000037.1_prot_CAB3368322.1_9305	biological_process	GO:0006351	29	9	7.605963657147811	DNA-templated transcription
lcl|CADEPI010000066.1_prot_CAB3371947.1_12930	molecular_function	GO:0005524	0	8	5.127132583684966	ATP binding
lcl|CADEPI010000002.1_prot_CAB3359968.1_951	cellular_component	GO:0019773	2	2	11.524807227649527	proteasome core complex, alpha-subunit complex
lcl|CADEPI010000066.1_prot_CAB3371947.1_12930	biological_process	GO:0006281	46	7	6.919242348331503	DNA repair
```
+ plots/Information_content_```filename```.jpeg represents Information Content (X axis) against probability density (Y axis) for BP, CC and MF.
  + Higher Information Content (right on the X axis) implies more specific GO term and a smaller (left on the X axis) value a more general GO term[^2].
  + Information content is calculated using the formula[^2]: $IC=-log{_2}{\(p\(t\)\)}$
    + Where $p(t)= 1 - \frac{'Offspring\ Count'}{'Offspring\ Count'\ +\ 'Ancestors\ Count'}$

![InformationContent_GCA_HMMER](https://user-images.githubusercontent.com/84094170/236842144-e9f0d29e-0267-4212-b25a-fab8e85d316b.jpeg)

[^1]: Klopfenstein, D.V., Zhang, L., Pedersen, B.S. et al. GOATOOLS: A Python library for Gene Ontology analyses. Sci Rep 8, 10872 (2018). https://doi.org/10.1038/s41598-018-28948-z
[^2]: Louie B, Bergen S, Higdon R, Kolker E. Quantifying protein function specificity in the gene ontology. Stand Genomic Sci. 2010 Mar 30;2(2):238-44. doi: [10.4056/sigs.561626](https://doi.org/10.4056%2Fsigs.561626). PMID: 21304708; PMCID: PMC3035283.

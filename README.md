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
+ outputs/filename.tsv has the following format
```
# Seq_ID\tCategory\tGOid\tDescendant_count\tDepth_level\tInformation_content\tDescription
lcl|CADEPI010000384.1_prot_CAB3384984.1_25967	molecular_function	GO:0005525	0	8	4.972614964445703	GTP binding
lcl|CADEPI010000015.1_prot_CAB3364279.1_5262	molecular_function	GO:0005509	4	5	4.801359277919819	calcium ion binding
lcl|CADEPI010000037.1_prot_CAB3368322.1_9305	molecular_function	GO:0003899	9	7	6.723510277717358	DNA-directed 5'-3' RNA polymerase activity
lcl|CADEPI010000066.1_prot_CAB3371947.1_12930	molecular_function	GO:0003697	8	5	6.509252576315588	single-stranded DNA binding
lcl|CADEPI010000037.1_prot_CAB3368322.1_9305	molecular_function	GO:0003677	130	4	3.142995268340097	DNA binding
lcl|CADEPI010000002.1_prot_CAB3359968.1_951	cellular_component	GO:0019773	2	2	7.988387636342149	proteasome core complex, alpha-subunit complex
lcl|CADEPI010000066.1_prot_CAB3371947.1_12930	biological_process	GO:0006281	46	7	4.796053325356956	DNA repair
lcl|CADEPI010000066.1_prot_CAB3371947.1_12930	molecular_function	GO:0005524	0	8	3.5538574947382617	ATP binding
```
+ plots/Information_content_```filename```.jpeg represents Information Content (X axis) against probability density (Y axis) for BP, CC and MF.
  + Higher Information Content (right on the X axis) implies more specific GO term and a smaller (left on the X axis) value a more general GO term.

![InformationContent_GCA_HMMER](https://user-images.githubusercontent.com/84094170/236820920-b35e2c43-9590-4b15-bfd8-7b57f9c5444e.jpeg)

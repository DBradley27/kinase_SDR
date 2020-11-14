# kinase_SDR
Code and data for the kinase SDR paper (2020)

The kinase PWMs and kinase sequence alignments are given in the following files:

```
PWMs_type_good.rds
report_kinases_al.fasta
```
The following script will take these two files as an input and generate SDR predictions from them:

```
SDR_pipeline_report.R
```
Running this R script successfully will require local installation of MAFFT, trimAl, and SPEER

The following Python script generates a mapping between the full-length multiple sequence alignment (MSA)
and the trimmed version of the alignment

```
corresponding_numbering.py
```

The following groups of files contains intramolecular contact data for CDK2 and PKA structures.
This is needed for the MultiRelief-3D method. The Python script maps this contact data to positions
in the trimmed alignment

```
CDK2_neighbours.txt
PKA_neighbours.txt
WHATIF_interchain_parser_neighbours_almap_notrim.py
```
For the SPEER, GroupSim, and MultiRelief-3D methods, the following Python scripts map the
results to positions in the untrimmed alignment and also to the human PKA sequence

```
mrelief_mapping_new_ensembl_notrim.py
SPEER_mapping_final_notrim.py
groupsim_mapping_notrim.py
```
Finally, I have also uploaded code written by Dr. Omar Wagih for the construction of kinase
PWMs from kinase substrate data

```
match-tm.r
```

# Furna: a structural database for FUnctions of RNA #

### 1. Download raw data ###
```bash
script/download_pdb.pl
script/download_rfam.pl
script/download_rnacentral.pl
script/download_goa.pl
script/download_sifts.pl
```

### 2. Curate data ###
```bash
script/curate_pdb.pl
script/curate_rfam.pl
script/curate_rnacentral.pl
script/curate_goa.pl
script/combine_chain.pl
```

### 3. Download data for each rna entry ###
```bash
script/download_pubmed.pl
script/download_assembly.pl
script/combine_interaction.pl
script/download_ligand.pl
```

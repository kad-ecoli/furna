#!/bin/bash
### 1. Download raw data ###
script/download_pdb.pl
script/download_rfam.pl
script/download_rnacentral.pl
script/download_goa.pl
script/download_sifts.pl

### 2. Curate data ###
script/curate_pdb.pl
script/curate_rfam.pl
script/curate_rnacentral.pl
script/curate_goa.pl
script/combine_chain.pl

### 3. Download data for each rna entry ###
script/download_pubmed.pl
script/download_assembly.pl
script/combine_interaction.pl
script/download_ligand.pl
script/combine_protein.pl
script/download_taxon.pl
script/fasta2taxon.pl
script/plot_goa.pl
script/combine_dna.pl
script/curate_blastdb.pl
script/download_uniprot.pl
script/download_attract.pl
script/make_summary.pl

# Documentation for all tab-separated tables available for download #

## rna.tsv.gz ##
All RNA entries in the FURNA database with the following columns:
1. 4-letter PDB ID
2. Chain ID (author assigned) in the asymmetric unit
3. Length (number of nucleotides with experimentally determined coordinates)
4. Resolution ("NA" for NMR structures without resolution)
5. Comma-separated list of Rfam familues
6. RNAcentral ID
7. Comma-separated list of PubMed ID
8. EC number
9. Comma-separated list of GO terms (Molecular Function), excluding parent terms
10. Comma-separated list of GO terms (Biological Process), excluding parent terms
11. Comma-separated list of GO terms (Cellular Component), excluding parent terms
12. NCBI Taxonomy ID
13. Sequence for nucleotides with experimentally determined coordinates
14. Dot-bracket format secondary structure assigned by CSSR
15. Dot-bracket format secondary structure assigned by DSSR
16. Title of the full PDB structure
17. Name of the PDB chain

## interaction.tsv.gz ##
All ligand-RNA interactions
1. 4-letter PDB ID
2. Chain ID (author assigned) of the "receptor" RNA in the asymmetric unit
3. Biological assembly (i.e., biounit) for which the ligand-RNA interaction is identified
4. Ligand ID. For small molecules, the ligand ID is the 3-letter Chemical Component Dictionary (CCD) ID used by the PDB database. For macromolecules, the  ligand ID is one of "protein", "dna", or "rna".
5. Chain ID of the ligand in the biological assembly.
6. Ligand index. For macromolecules, this value is always 0. For small molecules, this values is 1 or greater. For example, the 5-th MG in the biological assembly that interact with the given "receptor" RNA chain will have a ligand index of 5.
7. The original residue numbers for nucleotides on the "receptor" RNA that interact with the ligand.
8. The residue numbers for nucleotides on the "receptor" RNA that interact with the ligand, renumbered such that the first nucleotide in the "receptor" RNA starts from 1.
9. The residue sequence number of the macromolecular ligand. If the ligand consist of multiple residues, the first and last residue is separated by a "~".
10. The amino acid/nucleotide sequence of the ligand. Empty for small molecule ligands.


## dna.tsv.gz ##
All DNAs that interact with at least one RNA
1. 4-letter PDB ID
2. Chain ID (author assigned) in the asymmetric unit
3. Length (number of nucleotides with experimentally determined coordinates)
4. NCBI Taxonomy ID
5. Sequence for amino acid residues with experimentally determined coordinates
6. Name of the PDB chain

## protein.tsv.gz ##
All proteins that interact with at least one RNA
1. 4-letter PDB ID
2. Chain ID (author assigned) in the asymmetric unit
3. Length (number of amino acid residues with experimentally determined coordinates)
4. UniProt accession ID
5. EC number
6. Comma-separated list of GO terms (Molecular Function), excluding parent terms
7. Comma-separated list of GO terms (Biological Process), excluding parent terms
8. Comma-separated list of GO terms (Cellular Component), excluding parent terms
9. NCBI Taxonomy ID
10. Sequence for amino acid residues with experimentally determined coordinates
11. Name of the PDB chain

## ligand.tsv.gz ##
All small molecules (regular ligands and metal ions) that interact with at least one RNA in FURNA
1. 3-letter Chemical Component Dictionary (CCD) ID used by the PDB database to uniquely identify a small molecule
2. Chemical formula of the small molecule
3. International Chemical Identifier (InChI)
4. InChIKey (hashed InChI)
5. Semicolon separated list of SMILES string
6. Name of the small molecule
7. ChEMBL ID
8. DrugBank ID
9. ZINC ID

## attract_fimo.tsv.gz ##
Match to ATtRACT motifs for protein binding
1. 4-letter PDB ID
2. Chain ID (author assigned) in the asymmetric unit
3. First nucleotide for the matched motif. The nucleotide number is the position in the sequence for nucleotides with experimentally determined coordinates. This may be different from the residue sequence number in the mmCIF or PDB file.
4. Last nucleotide for the matched motif.
5. q-value of the sequence-motif match
6. The local region of the nucleotide sequence that matches the motif. Nucleotide U is converted to T.
7. ATtRACT motif ID.
8. Gene name of the RNA-binding protein for the motif.
9. Gene ID of the RNA-binding protein for the motif.
10. Semicolon separated list of PubMed ID for the citation of the motif. "PubMed ID is not available" if no citation is found.


## parent.tsv.gz ##
Full GO terms for each RNA, including the parent terms
1. 4-letter PDB ID
2. Chain ID (author assigned) in the asymmetric unit
3. Comma-separated list of GO terms (Molecular Function), including parent terms
4. Comma-separated list of GO terms (Biological Process), including parent terms
5. Comma-separated list of GO terms (Cellular Component), including parent terms

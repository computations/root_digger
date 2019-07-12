# Study

- Title: Phylogeny with introgression in Habronattus jumping spiders (Araneae:
  Salticidae)
- Dryad: https://www.datadryad.org/resource/doi:10.5061/dryad.5kg33

Original data files are placed in the `original` directory.

# Modifications

The trees were extracted from the file
`Leduc-Robert&Maddison_Habronattus-Trees.nex` and placed into the following
files:

- `ml_tree_from_primary_concat_1877_loci -> ml_primary.tree`
- `ml_tree_from_noncoding_loci_remnant_matrix -> ml_noncoding.tree`
- `ml_tree_from_missing_species_remnant_matrix -> ml_missing_species.tree`
- `ml_tree_from_concat._mitochondrial_loci -> ml_mitochondrial.tree`

And I renamed and converted the alignments:

- `Leduc-Robert&Maddison_Habronattus-MitochondrialTranscriptomes.nex -> mitochondrial.fasta`
- `Leduc-Robert&Maddison_Habronattus-MS-LociWithMissingSpecies1019Nuclear.nex -> missing_species.fasta`
- `Leduc-Robert&Maddison_Habronattus-NL-NoncodingLoci236Nuclear.nex -> noncoding.fasta`
- `Leduc-Robert&Maddison_Habronattus-Primary-1877NuclearGenes.nex -> primary.fasta`

Finally, the pairing between tree and alignment is

+---------------------------+-------------------------+
| tree file name            | alignment file name     |
+===========================+=========================+
| `ml_primary.tree`         | `primary.fasta`         |
+---------------------------+-------------------------+
| `ml_missing_species.tree` | `missing_species.fasta` |
+---------------------------+-------------------------+
| `ml_mitochondrial.tree`   | `mitochondrial.fasta`   |
+---------------------------+-------------------------+
| `ml_noncoding.tree`       | `noncoding.fasta`       |
+---------------------------+-------------------------+

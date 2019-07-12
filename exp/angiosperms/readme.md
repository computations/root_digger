# Study

- Title: Phylogenomics resolves the deep phylogeny of seed plants and indicates
  partial convergent or homoplastic evolution between Gnetales and angiosperms
- Dryad: https://www.datadryad.org/resource/doi:10.5061/dryad.f7f57

# Modifications

The following taxa were identfied as an outgroup:

- `Lja`
- `Paq`
- `Cga`

These taxa were used to reroot the tree in `Dendroscope`, and then removed.
They were also removed from the corrispoding alignment file. For convience the
files have been renamed:

- `CDS12_FcC_supermatrix.fas -> cds12_outgroup_removed.fasta`
- `CDS_FcC_supermatrix.fas -> cds_outgroup_removed.fasta`
- `RAxML_bipartitions.CDS12_partition -> cds12_outgroup_removed.tree`
- `RAxML_bipartitions.CDS_FcC_partition -> cds_outgroup_removed.tree`

And the data files are paired like in the following table

+-------------------------------+--------------------------------+
| tree file name                | alignment file name            |
+===============================+================================+
| `cds_outgroup_removed.tree`   | `cds_outgroup_removed.fasta`   |
+-------------------------------+--------------------------------+
| `cds12_outgroup_removed.tree` | `cds12_outgroup_removed.fasta` |
+-------------------------------+--------------------------------+

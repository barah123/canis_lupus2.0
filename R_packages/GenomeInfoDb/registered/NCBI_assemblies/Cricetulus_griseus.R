ORGANISM <- "Cricetulus griseus"

### List of assemblies by date.
ASSEMBLIES <- list(
    ## 109,152 sequences.
    list(assembly="CriGri_1.0",
         assembly_level="Scaffold",
         date="2011/08/23",
         ## GCF_000223135.1 and GCA_000223135.1 are identical.
         assembly_accession="GCF_000223135.1",
         circ_seqs="MT"),

    ## 52,710 sequences.
    list(assembly="C_griseus_v1.0",
         assembly_level="Scaffold",
         date="2013/07/12",
         ## GCF_000419365.1 and GCA_000419365.1 are identical.
         assembly_accession="GCA_000419365.1",   # criGri1
         circ_seqs=character(0)),

    ## 647 sequences.
    list(assembly="CriGri-PICRH-1.0",
         assembly_level="Chromosome",
         date="2020/06/11",
         ## GCF_003668045.3 and GCA_003668045.2 are identical.
         assembly_accession="GCF_003668045.3",
         circ_seqs=character(0))
)


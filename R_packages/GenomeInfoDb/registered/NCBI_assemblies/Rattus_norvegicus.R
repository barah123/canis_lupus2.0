ORGANISM <- "Rattus norvegicus"

### List of assemblies by date.
ASSEMBLIES <- list(
    ## 471 sequences.
    list(assembly="RGSC_v3.4",  # rn4
         assembly_level="Chromosome",
         date="2004/12/13",
         extra_info=c(strain="BN/SsNHsdMCW"),
         assembly_accession="GCF_000001895.3",
         circ_seqs="MT"),

    ## 2739 sequences.
    list(assembly="Rnor_5.0",  # rn5
         assembly_level="Chromosome",
         date="2012/03/16",
         extra_info=c(strain="BN/SsNHsdMCW"),
         assembly_accession="GCA_000001895.3",  # same as GCF_000001895.4
         circ_seqs="MT"),

    ## 955 sequences.
    list(assembly="Rnor_6.0",  # rn6
         assembly_level="Chromosome",
         date="2014/07/01",
         extra_info=c(strain="mixed"),
         assembly_accession="GCA_000001895.4",  # same as GCF_000001895.5
         circ_seqs="MT"),

    ## 176 sequences.
    list(assembly="mRatBN7.2",  # rn7
         assembly_level="Chromosome",
         date="2020/11/10",
         extra_info=c(strain="BN/NHsdMcwi"),
         assembly_accession="GCF_015227675.2",  # same as GCA_015227675.2
         circ_seqs="MT"),

    ## 77 sequences.
    list(assembly="GRCr8",  # will this become rn8?
         assembly_level="Chromosome",
         date="2024/01/31",
         extra_info=c(strain="BN/NHsdMcwi"),
         assembly_accession="GCF_036323735.1",  # same as GCA_036323735.1
         circ_seqs="chrM")
)


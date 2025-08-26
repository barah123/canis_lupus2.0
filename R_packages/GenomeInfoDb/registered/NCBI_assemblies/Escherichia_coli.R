ORGANISM <- "Escherichia coli"

### List of assemblies first by serotype, then by strain, then by substrain,
### then by date.
ASSEMBLIES <- list(

    ## --- strain: 536 ---

    ## Current scientific name: Escherichia coli 536
    ## NCBI Taxonomy ID: 362663

    list(assembly="ASM1330v1",
         assembly_level="Complete Genome",
         date="2006/07/13",
         extra_info=c(strain="536"),
         assembly_accession="GCF_000013305.1",
         circ_seqs="CP000247.1"),

    ## --- strain: APEC O1 ---

    ## Current scientific name: Escherichia coli APEC O1
    ## NCBI Taxonomy ID: 405955

    list(assembly="ASM1484v1",
         assembly_level="Complete Genome",
         date="2006/10/04",
         extra_info=c(strain="APEC O1"),
         assembly_accession="GCF_000014845.1",
         circ_seqs=c("CP000468.1", "pAPEC-O1-ColBM", "pAPEC-O1-R")),

    list(assembly="APEC O1",
         assembly_level="Contig",
         date="2020/08/04",
         extra_info=c(strain="APEC O1"),
         assembly_accession="GCF_902880315.1",
         circ_seqs=character(0)),

    ## --- strain: ATCC 8739 ---

    ## Current scientific name: Escherichia coli ATCC 8739
    ## NCBI Taxonomy ID: 481805

    list(assembly="ASM1938v1",
         assembly_level="Complete Genome",
         date="2008/03/11",
         extra_info=c(strain="ATCC 8739"),
         assembly_accession="GCF_000019385.1",
         circ_seqs="CP000946.1"),

    list(assembly="ASM359159v1",
         assembly_level="Complete Genome",
         date="2018/09/24",
         extra_info=c(strain="ATCC 8739"),
         assembly_accession="GCF_003591595.1",
         circ_seqs="CP022959.1"),

    list(assembly="ASM1686447v1",
         assembly_level="Complete Genome",
         date="2021/02/11",
         extra_info=c(strain="ATCC 8739"),
         assembly_accession="GCF_016864475.1",
         circ_seqs=c("CP043852.1", "unnamed")),

    list(assembly="ASM2601678v1",
         assembly_level="Complete Genome",
         date="2022/11/07",
         extra_info=c(strain="ATCC 8739"),
         assembly_accession="GCF_026016785.1",
         circ_seqs="CP033020.1"),

    ## --- strain: CFT073 ---

    ## Current scientific name: Escherichia coli CFT073
    ## NCBI Taxonomy ID: 199310

    list(assembly="ASM744v1",
         assembly_level="Complete Genome",
         date="2002/12/06",
         extra_info=c(strain="CFT073"),
         assembly_accession="GCF_000007445.1",
         circ_seqs="AE014075.1"),

    list(assembly="ASM1426294v1",
         assembly_level="Complete Genome",
         date="2020/08/21",
         extra_info=c(strain="CFT073"),
         assembly_accession="GCF_014262945.1",
         circ_seqs="GNE_CFT073_1"),

    ## --- strain: HS ---

    ## Current scientific name: Escherichia coli HS
    ## NCBI Taxonomy ID: 331112

    list(assembly="ASM1776v1",
         assembly_level="Complete Genome",
         date="2007/09/10",
         extra_info=c(strain="HS"),
         assembly_accession="GCF_000017765.1",
         circ_seqs="CP000802.1"),

    ## --- strain: K-12 ---

    ## Substrain: DH10B
    ## Current scientific name: Escherichia coli str. K-12 substr. DH10B
    ## NCBI Taxonomy ID: 316385
    list(assembly="ASM1942v1",
         assembly_level="Complete Genome",
         date="2008/03/14",
         extra_info=c(strain="K-12", substrain="DH10B"),
         assembly_accession="GCF_000019425.1",
         circ_seqs="CP000948.1"),

    ## Substrain: MG1655
    ## Current scientific name: Escherichia coli str. K-12 substr. MG1655
    ## NCBI Taxonomy ID: 511145
    list(assembly="ASM584v2",
         assembly_level="Complete Genome",
         date="2013/09/26",
         extra_info=c(strain="K-12", substrain="MG1655"),
         assembly_accession="GCF_000005845.2",
         circ_seqs="U00096.3"),

    ## Substrain: W3110
    ## Current scientific name: Escherichia coli str. K-12 substr. W3110
    ## NCBI Taxonomy ID: 316407
    list(assembly="ASM1024v1",
         assembly_level="Complete Genome",
         date="2006/01/23",
         extra_info=c(strain="K-12", substrain="W3110"),
         assembly_accession="GCF_000010245.2",
         circ_seqs="AP009048.1"),

    ## --- strain: SMS-3-5 ---

    ## Current scientific name: Escherichia coli SMS-3-5 
    ## NCBI Taxonomy ID: 439855

    list(assembly="ASM1964v1",
         assembly_level="Complete Genome",
         date="2008/03/20",
         extra_info=c(strain="SMS-3-5"),
         assembly_accession="GCF_000019645.1",
         circ_seqs=c("CP000970.1", "pSMS35_130", "pSMS35_3",
                     "pSMS35_4", "pSMS35_8")),

    ## --- strain: UTI89 ---

    ## Current scientific name: Escherichia coli UTI89
    ## NCBI Taxonomy ID: 364106

    list(assembly="ASM1326v1",
         assembly_level="Complete Genome",
         date="2006/04/05",
         extra_info=c(strain="UTI89"),
         assembly_accession="GCF_000013265.1",
         circ_seqs=c("CP000243.1", "pUTI89")),

    list(assembly="ASM1493087v1",
         assembly_level="Complete Genome",
         date="2020/10/25",
         extra_info=c(strain="UTI89"),
         assembly_accession="GCF_014930875.1",
         circ_seqs="CP062985.1"),

    list(assembly="ASM1564476v1",
         assembly_level="Complete Genome",
         date="2020/11/22",
         extra_info=c(strain="UTI89"),
         assembly_accession="GCF_015644765.1",
         circ_seqs="CP064825.1"),

    ## ========================= serotype O139:H28 ========================

    ## --- strain: E24377A ---

    ## Current scientific name: Escherichia coli O139:H28 str. E24377A
    ## NCBI Taxonomy ID: 331111

    list(assembly="ASM1774v1",
         assembly_level="Complete Genome",
         date="2007/09/11",
         extra_info=c(serotype="O139:H28", strain="E24377A"),
         assembly_accession="GCF_000017745.1",
         circ_seqs=c("CP000800.1", "pETEC_5", "pETEC_6", "pETEC_35",
                     "pETEC_73", "pETEC_74", "pETEC_80")),

    list(assembly="ena-yuan-GCF_000017745.1",
         assembly_level="Contig",
         date="2024/03/04",
         extra_info=c(serotype="O139:H28", strain="E24377A"),
         assembly_accession="GCA_963860295.1",
         circ_seqs=character(0)),

    ## ========================= serotype O157:H7 =========================

    ## --- strain: EDL933 ---

    ## Current scientific name: Escherichia coli O157:H7 str. EDL933
    ## NCBI Taxonomy ID: 155864

    list(assembly="ASM666v1",
         assembly_level="Chromosome",
         date="2004/12/06",
         extra_info=c(serotype="O157:H7", strain="EDL933"),
         assembly_accession="GCF_000006665.1",
         circ_seqs=c("AE005174.2", "pO157")),

    ## --- strain: Sakai ---

    ## Substrain: RIMD 0509952
    ## Current scientific name: Escherichia coli O157:H7 str. Sakai
    ## NCBI Taxonomy ID: 386585
    list(assembly="ASM886v2",
         assembly_level="Complete Genome",
         date="2018/06/08",
         extra_info=c(serotype="O157:H7",
                      strain="Sakai", substrain="RIMD 0509952"),
         assembly_accession="GCF_000008865.2",
         circ_seqs=c("BA000007.3", "pOSAK1", "pO157"))
)


# =============================================================================
# CanisLupus 2.0 - Microbiome Analysis Dashboard
# Author: Philip Appiah
# =============================================================================

# ── Package Installation (run once) ──────────────────────────────────────────
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("phyloseq", "ggtree", "metacoder", "microbiome", "DESeq2"))
# install.packages(c("shiny","plotly","DT","vegan","ape","ggplot2","tidyverse",
#   "networkD3","heatmaply","patchwork","RColorBrewer","viridis","igraph",
#   "cowplot","ggdendro","dendextend","broom","cluster","reshape2","scales"))

# ── Libraries ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(phyloseq)
  library(microbiome)
  library(plotly)
  library(DT)
  library(vegan)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(tidyverse)
  library(metacoder)
  library(networkD3)
  library(heatmaply)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(igraph)
  library(cowplot)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(scales)
})

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Remove singletons (taxa present in only 1 read across all samples)
remove_singletons <- function(ps) {
  prune_taxa(taxa_sums(ps) > 1, ps)
}

#' CLR transformation (robust for compositional data)
clr_transform <- function(ps) {
  microbiome::transform(ps, "clr")
}

#' Bray-Curtis normalised object (compositional)
compositional <- function(ps) {
  microbiome::transform(ps, "compositional")
}

#' Safe tax_glom: skips if rank missing
safe_tax_glom <- function(ps, rank) {
  if (rank %in% rank_names(ps)) tax_glom(ps, rank) else ps
}

wolf_pal <- c("#2D2926","#D4AF37","#8B5A2B","#6D7B8D","#A0856C",
              "#E5E4E2","#5C4033","#B8A99A","#4A6741","#7B8F6A")

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      :root {
        --wolf-dark:  #2D2926;
        --wolf-gold:  #D4AF37;
        --wolf-warm:  #8B5A2B;
        --wolf-cool:  #6D7B8D;
        --wolf-light: #E5E4E2;
      }
      body { font-family: Arial, sans-serif; background-color: #F0EFED; }
      .main-container { background: rgba(255,255,255,.95); border-radius:10px;
                        padding:20px; margin-top:15px;
                        box-shadow:0 3px 8px rgba(0,0,0,.15); }
      .sidebar-panel  { background:#fff; border-radius:8px; padding:14px;
                        border-left:4px solid var(--wolf-warm);
                        box-shadow:0 2px 4px rgba(0,0,0,.08); }
      .main-panel     { background:#fff; border-radius:8px; padding:18px;
                        box-shadow:0 2px 4px rgba(0,0,0,.08); }
      .btn-primary    { background-color:var(--wolf-dark); border-color:var(--wolf-dark); color:#fff; }
      .btn-primary:hover { background-color:var(--wolf-warm); border-color:var(--wolf-warm); }
      .btn-success    { background-color:#4A6741; border-color:#4A6741; color:#fff; }
      h4 { color:var(--wolf-dark); border-bottom:2px solid var(--wolf-warm); padding-bottom:4px; }
      .tabbable > .nav > li > a { color:var(--wolf-cool); }
      .tabbable > .nav > li.active > a { color:var(--wolf-dark); font-weight:bold;
                                         border-bottom:2px solid var(--wolf-gold); }
      .info-box { background:#fffbf0; border:1px solid var(--wolf-gold); border-radius:6px;
                  padding:10px; margin-bottom:12px; font-size:13px; color:#555; }
      .footer { background:var(--wolf-dark); color:var(--wolf-light); padding:18px;
                margin-top:25px; text-align:center; font-size:13px;
                border-top:3px solid var(--wolf-gold); }
      .footer a { color:var(--wolf-gold); text-decoration:none; margin:0 12px; font-weight:bold; }
      .footer a:hover { color:var(--wolf-light); }
      hr { border-top:1px solid var(--wolf-warm); }
    "))
  ),
  
  div(class = "main-container",
      
      # ── Header ──────────────────────────────────────────────────────────────
      div(style = "display:flex; align-items:center; margin-bottom:18px;",
          img(src = "canis2.png", height = "72px",
              style = "margin-right:18px; border:3px solid var(--wolf-gold); border-radius:50%;",
              onerror = "this.style.display='none'"),
          div(
            h2("CanisLupus 2.0", style = "color:var(--wolf-dark); margin-bottom:2px;"),
            h4("Microbiome Analysis Dashboard", style = "color:var(--wolf-cool); font-style:italic; border:none; margin-top:6px;")
          )
      ),
      
      # ── Navigation ──────────────────────────────────────────────────────────
      navlistPanel(
        id = "mainTabs", widths = c(2, 10),
        
        # ────────────────────────────────────────────────────────────────────
        "Data",
        # ────────────────────────────────────────────────────────────────────
        
        tabPanel("Data Upload & QC",
                 sidebarLayout(
                   sidebarPanel(width = 3, class = "sidebar-panel",
                                h4("Upload Files"),
                                div(class="info-box","All files are optional; demo data loads automatically if none provided."),
                                fileInput("asvFile",  "ASV Table (CSV, taxa as rows)", accept=".csv"),
                                fileInput("taxFile",  "Taxonomy Table (CSV)",          accept=".csv"),
                                fileInput("metaFile", "Metadata Table (CSV)",          accept=".csv"),
                                fileInput("treeFile", "Phylogenetic Tree (.nwk/.tree)", accept=c(".tree",".tre",".nwk",".txt")),
                                hr(),
                                h4("Singleton Removal"),
                                div(class="info-box","Singletons = taxa with total read count = 1 across all samples. Removing them reduces noise and spurious diversity estimates."),
                                checkboxInput("removeSingletons", "Remove Singletons", value = TRUE),
                                hr(),
                                h4("Sample Filtering"),
                                numericInput("minReads", "Min reads per sample:", value = 1000, min = 0, step = 100),
                                hr(),
                                actionButton("update", "Load & Process Data", class="btn-primary", style="width:100%;"),
                                br(), br(),
                                downloadButton("downloadPS", "Download Filtered Data Summary", class="btn-success", style="width:100%;")
                   ),
                   mainPanel(width = 9, class="main-panel",
                             h4("Dataset Summary"),
                             verbatimTextOutput("dataSummary"),
                             hr(),
                             h4("QC: Read Counts per Sample"),
                             plotlyOutput("readCountPlot", height="350px"),
                             hr(),
                             h4("Sample Metadata Table"),
                             DTOutput("sampleTable")
                   )
                 )
        ),
        
        # ────────────────────────────────────────────────────────────────────
        "Diversity",
        # ────────────────────────────────────────────────────────────────────
        
        tabPanel("Rarefaction Curve",
                 sidebarLayout(
                   sidebarPanel(width = 3, class="sidebar-panel",
                                h4("Rarefaction Settings"),
                                div(class="info-box","Rarefaction curves show whether sequencing depth was sufficient to capture community diversity. A plateau indicates saturation."),
                                selectInput("regionRare", "Color samples by:", choices=NULL),
                                numericInput("rareStep", "Step size:", value=100, min=10, step=50),
                                numericInput("rareMax",  "Max depth (0 = auto):", value=0, min=0),
                                checkboxInput("showRareLine","Show rarefaction depth line", value=TRUE),
                                numericInput("rareDepth", "Rarefaction depth:", value=10000, min=100),
                                actionButton("updateRare","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("rarefactionCurve", height="550px"),
                             hr(),
                             verbatimTextOutput("rareSummary")
                   )
                 )
        ),
        
        tabPanel("Alpha Diversity",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Alpha Diversity"),
                                div(class="info-box","Alpha diversity measures within-sample species richness and evenness. Multiple indices capture different aspects."),
                                checkboxGroupInput("alphaMeasure","Diversity indices:",
                                                   choices = c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Fisher"),
                                                   selected = c("Shannon","Observed")),
                                selectInput("regionAlpha","Group by:", choices=NULL),
                                checkboxInput("alphaStats","Show Kruskal-Wallis p-value", value=TRUE),
                                actionButton("updateAlpha","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("alphaDiv",      height="450px"),
                             plotlyOutput("alphaDivBoxplot",height="450px"),
                             verbatimTextOutput("alphaStats")
                   )
                 )
        ),
        
        tabPanel("Beta Diversity",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Beta Diversity"),
                                div(class="info-box","Beta diversity measures between-sample dissimilarity. PERMANOVA tests whether group centroids differ significantly."),
                                selectInput("betaMethod","Ordination method:",
                                            choices=c("PCoA","NMDS","RDA","CCA"), selected="PCoA"),
                                selectInput("betaDistance","Distance metric:",
                                            choices=c("bray","jaccard","unifrac","wunifrac","euclidean"), selected="bray"),
                                selectInput("regionBeta","Color/group by:", choices=NULL),
                                checkboxInput("betaEllipse","Add 95% ellipses", value=TRUE),
                                checkboxInput("betaPermanova","Run PERMANOVA", value=TRUE),
                                actionButton("updateBeta","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("betaPlot",  height="500px"),
                             plotlyOutput("nmdsPlot",  height="500px"),
                             verbatimTextOutput("permanovaResult")
                   )
                 )
        ),
        
        # ────────────────────────────────────────────────────────────────────
        "Taxonomy",
        # ────────────────────────────────────────────────────────────────────
        
        tabPanel("Stacked Bar Plots",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Bar Plot Options"),
                                selectInput("taxLevelBar","Taxonomic level:",
                                            choices=c("Phylum","Class","Order","Family","Genus","Species"), selected="Phylum"),
                                selectInput("regionBar","Facet by:", choices=NULL),
                                sliderInput("topNBar","Top N taxa:", min=3, max=30, value=10),
                                selectInput("barTransform","Abundance type:",
                                            choices=c("Relative"="relative","Absolute"="absolute"), selected="relative"),
                                actionButton("updateBar","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("taxaBarplot",        height="500px"),
                             plotlyOutput("relativeAbundancePlot",height="450px")
                   )
                 )
        ),
        
        tabPanel("Interactive Pie Chart",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Pie Chart Options"),
                                selectInput("taxLevelPie","Taxonomic level:",
                                            choices=c("Phylum","Class","Order","Family","Genus","Species"), selected="Phylum"),
                                selectInput("regionPie","Group by:", choices=NULL),
                                actionButton("updatePie","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("pieChart", height="700px")
                   )
                 )
        ),
        
        tabPanel("Core Microbiome",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Core Microbiome"),
                                div(class="info-box","Core = taxa present above detection threshold in a defined fraction of samples. These are the most ecologically consistent members."),
                                selectInput("taxLevelCore","Taxonomic level:",
                                            choices=c("Phylum","Class","Order","Family","Genus","Species"), selected="Genus"),
                                selectInput("regionCore","Group by:", choices=NULL),
                                sliderInput("prevalenceCore","Prevalence threshold (%):", min=5,  max=100, value=50),
                                sliderInput("detectionCore", "Detection threshold (%):",  min=0.001, max=10, value=0.1, step=0.01),
                                actionButton("updateCore","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("coreHeatmap", height="550px"),
                             hr(),
                             verbatimTextOutput("coreSummary")
                   )
                 )
        ),
        
        # ────────────────────────────────────────────────────────────────────
        "Phylogenetics",
        # ────────────────────────────────────────────────────────────────────
        tabPanel("Phylogenetic Tree",
                 sidebarLayout(
                   sidebarPanel(width = 3, class = "sidebar-panel",
                                
                                h4("Tree Options"),
                                
                                div(class="info-box",
                                    "Displays the evolutionary relationships among ASVs/OTUs.
           A random tree is generated if none is uploaded."),
                                
                                selectInput("treeLayout","Layout:",
                                            choices=c("rectangular","circular","fan","radial"),
                                            selected="rectangular"),
                                
                                # ADD THIS BLOCK HERE
                                selectInput("treeColorBy","Color tips by:",
                                            choices=NULL),
                                
                                sliderInput("treeTipSize","Tip label size:",
                                            min=0, max=5, value=2, step=0.5),
                                
                                numericInput("treePruneN","Show top N taxa (0=all):",
                                             value=50, min=0),
                                
                                actionButton("updateTree","Update", class="btn-primary")
                   ),
                   
                   mainPanel(width = 9,
                             plotOutput("phylogeneticTree", height="750px")
                   )
                 )
        ),
        
        # ────────────────────────────────────────────────────────────────────
        "Network & Clustering",
        # ────────────────────────────────────────────────────────────────────
        
        tabPanel("Interactive Heatmap",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Heatmap Options"),
                                div(class="info-box","Clustered heatmap reveals sample and taxon groupings. CLR transformation stabilises variance for compositional data."),
                                selectInput("taxLevelHeatmap","Taxonomic level:",
                                            choices=c("Phylum","Class","Order","Family","Genus"), selected="Genus"),
                                selectInput("regionHeatmap","Annotate by:", choices=NULL),
                                selectInput("hmTransform","Transformation:",
                                            choices=c("CLR"="clr","Relative"="compositional","Log10p"="log10p"), selected="clr"),
                                sliderInput("hmTopN","Top N taxa:", min=5, max=50, value=20),
                                actionButton("updateHeatmap","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("interactiveHeatmap", height="750px")
                   )
                 )
        ),
        
        tabPanel("Dendrogram",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Dendrogram Options"),
                                div(class="info-box","Hierarchical clustering of samples based on community composition. Ward.D2 + Bray-Curtis is the standard microbiome approach."),
                                selectInput("distanceMethod","Distance method:",
                                            choices=c("bray","jaccard","euclidean","manhattan","canberra"), selected="bray"),
                                selectInput("clusterMethod","Clustering algorithm:",
                                            choices=c("ward.D2","complete","average","single","mcquitty"), selected="ward.D2"),
                                selectInput("dendroColorBy","Color labels by:", choices=NULL),
                                checkboxInput("dendroShowBar","Show abundance bar alongside", value=FALSE),
                                actionButton("updateDendro","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("dendrogram", height="650px"),
                             conditionalPanel("input.dendroShowBar",
                                              plotlyOutput("dendroBar", height="400px"))
                   )
                 )
        ),
        
        tabPanel("Correlation Network",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Network Options"),
                                div(class="info-box","Spearman correlations between taxa abundance profiles. Edges represent statistically meaningful co-occurrence or mutual-exclusion relationships."),
                                selectInput("taxLevelNetwork","Taxonomic level:",
                                            choices=c("Phylum","Class","Order","Family","Genus"), selected="Genus"),
                                sliderInput("corThreshold","|Correlation| threshold:", min=0.3, max=0.99, value=0.6, step=0.05),
                                selectInput("corMethod","Correlation method:",
                                            choices=c("spearman","pearson"), selected="spearman"),
                                numericInput("netTopN","Top N taxa by abundance:", value=30, min=5, max=80),
                                checkboxInput("netNegEdges","Show negative correlations (red)", value=TRUE),
                                actionButton("updateNetwork","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             forceNetworkOutput("correlationNetwork", height="700px"),
                             hr(),
                             verbatimTextOutput("networkSummary")
                   )
                 )
        ),
        
        # ────────────────────────────────────────────────────────────────────
        "Statistics",
        # ────────────────────────────────────────────────────────────────────
        
        tabPanel("Differential Abundance",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Differential Abundance"),
                                div(class="info-box","Kruskal-Wallis + Benjamini-Hochberg FDR correction identifies taxa significantly different across groups."),
                                selectInput("taxLevelDA","Taxonomic level:",
                                            choices=c("Phylum","Class","Order","Family","Genus"), selected="Genus"),
                                selectInput("regionDA","Grouping variable:", choices=NULL),
                                numericInput("daFDR","FDR threshold:", value=0.05, min=0.001, max=0.2, step=0.01),
                                actionButton("updateDA","Run Analysis", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("daVolcano",   height="500px"),
                             hr(),
                             DTOutput("daTable")
                   )
                 )
        ),
        
        tabPanel("Transformation Explorer",
                 sidebarLayout(
                   sidebarPanel(width=3, class="sidebar-panel",
                                h4("Transformation"),
                                div(class="info-box","Visualise the effect of different normalisations on your data distribution before downstream analysis."),
                                selectInput("transMethod","Transformation:",
                                            choices=c("Raw"="raw","Relative"="compositional","CLR"="clr",
                                                      "Hellinger"="hellinger","Log10p"="log10p"), selected="clr"),
                                selectInput("transColorBy","Color by:", choices=NULL),
                                actionButton("updateTrans","Update", class="btn-primary")
                   ),
                   mainPanel(width=9,
                             plotlyOutput("transPCoA",  height="500px"),
                             plotlyOutput("transHist",  height="350px")
                   )
                 )
        )
        
      ) # end navlistPanel
  ), # end main-container
  
  # ── Footer ──────────────────────────────────────────────────────────────────
  div(class="footer",
      h4("Contact", style="color:var(--wolf-gold);"),
      a(href="mailto:pyappiah561@gmail.com","✉ Email"),
      a(href="https://www.linkedin.com/in/philip-appiah","in LinkedIn"),
      a(href="https://github.com/barah123","⌥ GitHub"),
      p(style="color:var(--wolf-gold); margin-top:10px;",
        "© 2025 CanisLupus 2.0 | Canis lupus Microbiome Project | Philip Appiah | Tel: +1 202 500 8302")
  )
)

# =============================================================================
# SERVER
# =============================================================================

server <- function(input, output, session) {
  
  # ── Reactive: processed phyloseq object ─────────────────────────────────────
  ps <- reactiveVal(NULL)
  
  # ── Update all group selectors when ps changes ────────────────────────────
  observe({
    req(ps())
    vars <- colnames(sample_data(ps()))
    for (id in c("regionBar","regionPie","regionRare","regionAlpha","regionBeta",
                 "regionCore","regionHeatmap","dendroColorBy","regionDA",
                 "transColorBy","treeColorBy")) {
      updateSelectInput(session, id, choices = c("None" = "", vars), selected = "")
    }
  })
  
  # ── Load Data ────────────────────────────────────────────────────────────────
  observeEvent(input$update, {
    withProgress(message = "Loading and processing data…", value = 0, {
      
      tryCatch({
        incProgress(0.1, detail = "Reading ASV table")
        asv <- if (!is.null(input$asvFile))
          read.csv(input$asvFile$datapath,  row.names=1, check.names=FALSE)
        else
          read.csv("demo_asv.csv",           row.names=1, check.names=FALSE)
        
        incProgress(0.2, detail = "Reading taxonomy")
        tax_mat <- if (!is.null(input$taxFile))
          as.matrix(read.csv(input$taxFile$datapath,  row.names=1))
        else
          as.matrix(read.csv("demo_taxonomy.csv",      row.names=1))
        
        incProgress(0.2, detail = "Reading metadata")
        meta <- if (!is.null(input$metaFile))
          read.csv(input$metaFile$datapath, row.names=1)
        else
          read.csv("demo_metadata.csv",     row.names=1)
        
        incProgress(0.1, detail = "Building phyloseq object")
        otu  <- otu_table(as.matrix(asv), taxa_are_rows=TRUE)
        tax  <- tax_table(tax_mat)
        sam  <- sample_data(meta)
        ps_obj <- phyloseq(otu, tax, sam)
        
        # Load or generate tree
        if (!is.null(input$treeFile)) {
          tree <- read.tree(input$treeFile$datapath)
          ps_obj <- merge_phyloseq(ps_obj, tree)
        } else {
          # Generate random tree so tree-based analyses always work
          set.seed(42)
          rand_tree <- ape::rtree(ntaxa(ps_obj), rooted=TRUE,
                                  tip.label=taxa_names(ps_obj))
          ps_obj <- merge_phyloseq(ps_obj, rand_tree)
        }
        
        incProgress(0.2, detail = "Applying filters")
        # Minimum reads filter
        ps_obj <- prune_samples(sample_sums(ps_obj) >= input$minReads, ps_obj)
        
        # Singleton removal
        n_before <- ntaxa(ps_obj)
        if (input$removeSingletons) {
          ps_obj <- remove_singletons(ps_obj)
        }
        n_after <- ntaxa(ps_obj)
        
        ps(ps_obj)
        
        incProgress(0.2, detail = "Done")
        showNotification(
          paste0("Data loaded! ", nsamples(ps_obj), " samples | ",
                 n_after, " taxa",
                 if (input$removeSingletons)
                   paste0(" (", n_before - n_after, " singletons removed)")
                 else ""),
          type="message", duration=6
        )
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type="error", duration=10)
      })
    })
  })
  
  # ── Download summary ─────────────────────────────────────────────────────
  output$downloadPS <- downloadHandler(
    filename = function() paste0("canis_lupus_summary_", Sys.Date(), ".txt"),
    content  = function(file) {
      req(ps())
      sink(file)
      cat("=== CanisLupus 2.0 - Filtered Dataset Summary ===\n\n")
      print(ps())
      cat("\nSample read counts:\n")
      print(sort(sample_sums(ps())))
      cat("\nTax ranks:\n")
      print(rank_names(ps()))
      sink()
    }
  )
  
  # ── Data Summary ─────────────────────────────────────────────────────────
  output$dataSummary <- renderPrint({
    req(ps())
    cat("=== Microbiome Dataset ===\n\n")
    cat("Samples    :", nsamples(ps()), "\n")
    cat("Taxa       :", ntaxa(ps()),    "\n")
    cat("Tax ranks  :", paste(rank_names(ps()), collapse=" > "), "\n")
    cat("Tree tips  :", ifelse(!is.null(phy_tree(ps(), errorIfNULL=FALSE)),
                               length(phy_tree(ps())$tip.label), "none"), "\n")
    cat("\nRead depth (min/median/max):",
        min(sample_sums(ps())), "/",
        round(median(sample_sums(ps()))), "/",
        max(sample_sums(ps())), "\n")
    cat("\nSample variables:\n")
    print(colnames(sample_data(ps())))
  })
  
  # ── QC Read Count Plot ───────────────────────────────────────────────────
  output$readCountPlot <- renderPlotly({
    req(ps())
    df <- data.frame(
      Sample = sample_names(ps()),
      Reads  = sample_sums(ps())
    ) %>% arrange(Reads)
    df$Sample <- factor(df$Sample, levels=df$Sample)
    
    plot_ly(df, x=~Sample, y=~Reads, type="bar",
            marker=list(color=wolf_pal[2])) %>%
      layout(title  = "Read Counts per Sample",
             xaxis  = list(tickangle=-45, title=""),
             yaxis  = list(title="Total Reads"),
             shapes = list(list(type="line", x0=0, x1=1, xref="paper",
                                y0=input$minReads, y1=input$minReads,
                                line=list(color="red", dash="dash"))))
  })
  
  output$sampleTable <- renderDT({
    req(ps())
    datatable(as.data.frame(sample_data(ps())),
              options=list(scrollX=TRUE, pageLength=10),
              class="cell-border stripe")
  })
  
  # ============================================================
  # RAREFACTION CURVE (improved: proper per-sample curves)
  # ============================================================
  output$rarefactionCurve <- renderPlotly({
    req(ps())
    input$updateRare
    isolate({
      withProgress(message="Computing rarefaction curves…", {
        otu_mat <- t(as(otu_table(ps()), "matrix"))
        maxd     <- if (input$rareMax > 0) input$rareMax else max(rowSums(otu_mat))
        step_sz  <- max(input$rareStep, 1)
        
        rc <- vegan::rarecurve(otu_mat, step=step_sz, tmax=maxd,
                               sample=input$rareDepth, label=FALSE)
        
        samp_data <- as.data.frame(sample_data(ps()))
        grp_var   <- if (!is.null(input$regionRare) && nchar(input$regionRare)>0 &&
                         input$regionRare %in% colnames(samp_data))
          input$regionRare else NULL
        
        plot_list <- lapply(seq_along(rc), function(i) {
          x <- attr(rc[[i]], "Subsample")
          y <- rc[[i]]
          grp <- if (!is.null(grp_var)) as.character(samp_data[[grp_var]][i]) else "All"
          data.frame(Sample=rownames(samp_data)[i], Reads=x, OTUs=y, Group=grp)
        })
        plot_df <- bind_rows(plot_list)
        
        p <- ggplot(plot_df, aes(x=Reads, y=OTUs, group=Sample, color=Group)) +
          geom_line(alpha=0.7, linewidth=0.8) +
          theme_minimal(base_size=13) +
          scale_color_brewer(palette="Set1") +
          labs(title="Rarefaction Curves", x="Sequencing Depth (Reads)", y="Observed OTUs",
               color=if(!is.null(grp_var)) grp_var else "")
        
        if (input$showRareLine) {
          p <- p + geom_vline(xintercept=input$rareDepth, linetype="dashed", color="red", linewidth=0.8)
        }
        ggplotly(p)
      })
    })
  })
  
  output$rareSummary <- renderPrint({
    req(ps())
    input$updateRare
    isolate({
      ss <- sort(sample_sums(ps()))
      cat("=== Read-depth summary ===\n")
      cat("Min:", min(ss), "\n")
      cat("1st Qu:", quantile(ss, 0.25), "\n")
      cat("Median:", median(ss), "\n")
      cat("Mean:", round(mean(ss)), "\n")
      cat("3rd Qu:", quantile(ss, 0.75), "\n")
      cat("Max:", max(ss), "\n")
      cat("\nSamples below rarefaction depth (", input$rareDepth, "):",
          sum(ss < input$rareDepth), "\n")
    })
  })
  
  # ============================================================
  # ALPHA DIVERSITY
  # ============================================================
  output$alphaDiv <- renderPlotly({
    req(ps())
    input$updateAlpha
    isolate({
      measures <- input$alphaMeasure
      if (length(measures)==0) measures <- "Shannon"
      p <- plot_richness(ps(), measures=measures) +
        theme_minimal(base_size=12) +
        scale_color_brewer(palette="Set1") +
        labs(title="Alpha Diversity Indices")
      if (!is.null(input$regionAlpha) && nchar(input$regionAlpha)>0 &&
          input$regionAlpha %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color=input$regionAlpha), size=3)
      } else {
        p <- p + geom_point(size=3)
      }
      ggplotly(p)
    })
  })
  
  output$alphaDivBoxplot <- renderPlotly({
    req(ps())
    input$updateAlpha
    isolate({
      if (is.null(input$regionAlpha) || nchar(input$regionAlpha)==0 ||
          !input$regionAlpha %in% colnames(sample_data(ps()))) return(NULL)
      measures <- input$alphaMeasure
      if (length(measures)==0) measures <- "Shannon"
      p <- plot_richness(ps(), x=input$regionAlpha, measures=measures) +
        geom_boxplot(aes_string(fill=input$regionAlpha), alpha=0.6, outlier.shape=NA) +
        geom_jitter(aes_string(color=input$regionAlpha), width=0.15, size=2) +
        theme_minimal(base_size=12) +
        scale_fill_brewer(palette="Set1") +
        scale_color_brewer(palette="Set1") +
        labs(title=paste("Alpha Diversity by", input$regionAlpha))
      ggplotly(p)
    })
  })
  
  output$alphaStats <- renderPrint({
    req(ps())
    input$updateAlpha
    isolate({
      if (is.null(input$regionAlpha) || nchar(input$regionAlpha)==0 ||
          !input$alphaStats) return(invisible(NULL))
      grp_var <- input$regionAlpha
      if (!grp_var %in% colnames(sample_data(ps()))) return(invisible(NULL))
      measures <- input$alphaMeasure
      if (length(measures)==0) measures <- "Shannon"
      rich <- estimate_richness(ps(), measures=measures)
      rich$Group <- as.factor(sample_data(ps())[[grp_var]])
      cat("=== Kruskal-Wallis Tests ===\n\n")
      for (m in measures) {
        if (m %in% colnames(rich)) {
          kw <- kruskal.test(rich[[m]] ~ rich$Group)
          cat(m, ": p =", round(kw$p.value, 4), "\n")
        }
      }
    })
  })
  
  # ============================================================
  # BETA DIVERSITY
  # ============================================================
  output$betaPlot <- renderPlotly({
    req(ps())
    input$updateBeta
    isolate({
      ord <- ordinate(ps(), method=input$betaMethod, distance=input$betaDistance)
      p   <- plot_ordination(ps(), ord, type="samples") +
        theme_minimal(base_size=12) +
        labs(title=paste(input$betaMethod, "(", input$betaDistance, ")"))
      if (!is.null(input$regionBeta) && nchar(input$regionBeta)>0 &&
          input$regionBeta %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color=input$regionBeta), size=3)
        if (input$betaEllipse)
          p <- p + stat_ellipse(aes_string(color=input$regionBeta), level=0.95)
      } else {
        p <- p + geom_point(color=wolf_pal[2], size=3)
      }
      ggplotly(p)
    })
  })
  
  output$nmdsPlot <- renderPlotly({
    req(ps())
    input$updateBeta
    isolate({
      suppressMessages({
        ord <- ordinate(ps(), method="NMDS", distance=input$betaDistance)
      })
      p <- plot_ordination(ps(), ord, type="samples") +
        theme_minimal(base_size=12) +
        labs(title=paste("NMDS (", input$betaDistance, ")"))
      if (!is.null(input$regionBeta) && nchar(input$regionBeta)>0 &&
          input$regionBeta %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color=input$regionBeta), size=3)
        if (input$betaEllipse)
          p <- p + stat_ellipse(aes_string(color=input$regionBeta), level=0.95)
      } else {
        p <- p + geom_point(color=wolf_pal[3], size=3)
      }
      ggplotly(p)
    })
  })
  
  output$permanovaResult <- renderPrint({
    req(ps())
    input$updateBeta
    isolate({
      if (!input$betaPermanova) return(invisible(NULL))
      grp_var <- input$regionBeta
      if (is.null(grp_var) || nchar(grp_var)==0 ||
          !grp_var %in% colnames(sample_data(ps()))) {
        cat("Select a grouping variable to run PERMANOVA.\n")
        return(invisible(NULL))
      }
      dist_mat <- phyloseq::distance(ps(), method=input$betaDistance)
      meta_df  <- as.data.frame(sample_data(ps()))
      set.seed(42)
      perm <- vegan::adonis2(dist_mat ~ meta_df[[grp_var]], permutations=999)
      cat("=== PERMANOVA (adonis2) ===\n")
      cat("Response variable:", grp_var, "\n")
      cat("Distance metric:  ", input$betaDistance, "\n\n")
      print(perm)
      
      # Homogeneity of dispersion
      bd <- vegan::betadisper(dist_mat, meta_df[[grp_var]])
      cat("\n=== Homogeneity of dispersion (betadisper) ===\n")
      print(anova(bd))
    })
  })
  
  # ============================================================
  # STACKED BAR PLOTS
  # ============================================================
  output$taxaBarplot <- renderPlotly({
    req(ps())
    input$updateBar
    isolate({
      ps_use <- if (input$barTransform=="relative")
        transform_sample_counts(ps(), function(x) x/sum(x)) else ps()
      top_t <- names(sort(taxa_sums(ps_use), decreasing=TRUE)[1:input$topNBar])
      ps_top <- prune_taxa(top_t, ps_use)
      p <- plot_bar(ps_top, fill=input$taxLevelBar) +
        theme_minimal(base_size=12) +
        labs(title=paste("Top", input$topNBar, input$taxLevelBar)) +
        scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(input$topNBar))
      if (!is.null(input$regionBar) && nchar(input$regionBar)>0 &&
          input$regionBar %in% colnames(sample_data(ps())))
        p <- p + facet_wrap(as.formula(paste("~", input$regionBar)), scales="free_x")
      ggplotly(p)
    })
  })
  
  output$relativeAbundancePlot <- renderPlotly({
    req(ps())
    input$updateBar
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelBar)
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      p <- plot_bar(ps_rel, fill=input$taxLevelBar) +
        theme_minimal(base_size=12) +
        labs(title=paste(input$taxLevelBar, "Relative Abundance")) +
        scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Set3"))(ntaxa(ps_rel)))
      if (!is.null(input$regionBar) && nchar(input$regionBar)>0 &&
          input$regionBar %in% colnames(sample_data(ps())))
        p <- p + facet_wrap(as.formula(paste("~", input$regionBar)), scales="free_x")
      ggplotly(p)
    })
  })
  
  # ============================================================
  # PIE CHART
  # ============================================================
  output$pieChart <- renderPlotly({
    req(ps())
    input$updatePie
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelPie)
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      melted  <- psmelt(ps_rel)
      grp_var <- input$regionPie
      lvl_col <- input$taxLevelPie
      
      if (!is.null(grp_var) && nchar(grp_var)>0 && grp_var %in% colnames(melted)) {
        agg <- melted %>%
          group_by(across(all_of(c(lvl_col, grp_var)))) %>%
          summarise(Abundance=mean(Abundance), .groups="drop")
        plot_ly(agg, labels=~get(lvl_col), values=~Abundance,
                color=~get(grp_var), type="pie",
                textposition="inside", textinfo="label+percent") %>%
          layout(title=paste(lvl_col, "by", grp_var))
      } else {
        agg <- melted %>% group_by(across(all_of(lvl_col))) %>%
          summarise(Abundance=mean(Abundance), .groups="drop")
        plot_ly(agg, labels=~get(lvl_col), values=~Abundance, type="pie",
                textposition="inside", textinfo="label+percent") %>%
          layout(title=paste(lvl_col, "Composition"))
      }
    })
  })
  
  # ============================================================
  # CORE MICROBIOME
  # ============================================================
  output$coreHeatmap <- renderPlotly({
    req(ps())
    input$updateCore
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelCore)
      core_ids <- microbiome::core_members(ps_glom,
                                           detection=input$detectionCore/100,
                                           prevalence=input$prevalenceCore/100)
      if (length(core_ids)==0) {
        return(plot_ly() %>%
                 add_annotations(text="No core taxa found – try lower thresholds",
                                 x=0.5, y=0.5, showarrow=FALSE) %>%
                 layout(title="Core Microbiome"))
      }
      ps_core <- prune_taxa(core_ids, ps_glom)
      ps_core <- transform_sample_counts(ps_core, function(x) x/sum(x))
      mat <- as(otu_table(ps_core), "matrix")
      rownames(mat) <- as.character(tax_table(ps_core)[, input$taxLevelCore])
      heatmaply(t(mat), colors=viridis(256),
                main=paste("Core Microbiome –", input$taxLevelCore))
    })
  })
  
  output$coreSummary <- renderPrint({
    req(ps())
    input$updateCore
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelCore)
      core_ids <- microbiome::core_members(ps_glom,
                                           detection=input$detectionCore/100,
                                           prevalence=input$prevalenceCore/100)
      cat("=== Core Microbiome ===\n")
      cat("Taxonomic level:", input$taxLevelCore, "\n")
      cat("Detection threshold:", input$detectionCore, "% relative abundance\n")
      cat("Prevalence threshold:", input$prevalenceCore, "% samples\n")
      cat("Core taxa found:", length(core_ids), "\n\n")
      if (length(core_ids) > 0) {
        cat("Core taxa names:\n")
        print(as.character(tax_table(ps_glom)[core_ids, input$taxLevelCore]))
      }
    })
  })
  
  # ============================================================
  # PHYLOGENETIC TREE (improved with ggtree options)
  # ============================================================
  output$phylogeneticTree <- renderPlot({
    
    req(ps())
    input$updateTree
    
    isolate({
      
      tree_ps <- ps()
      
      # prune taxa if user requests
      if (input$treePruneN > 0 && ntaxa(tree_ps) > input$treePruneN) {
        
        top_t <- names(sort(taxa_sums(tree_ps), decreasing = TRUE))[1:input$treePruneN]
        
        tree_ps <- prune_taxa(top_t, tree_ps)
        
        # IMPORTANT: prune tree to match taxa
        phy_tree(tree_ps) <- ape::keep.tip(
          phy_tree(tree_ps),
          taxa_names(tree_ps)
        )
      }
      
      tree <- phy_tree(tree_ps)
      
      p <- ggtree::ggtree(tree, layout = input$treeLayout) +
        ggtree::geom_tiplab(size = input$treeTipSize) +
        ggtree::theme_tree2() +
        ggplot2::ggtitle("Phylogenetic Tree")
      
      
      print(p)
      
    })
    
  })
  
  # ============================================================
  # HEAT TREE (metacoder – improved)
  # ============================================================
  output$heatTree <- renderPlot({
    req(ps())
    input$updateHeatTree
    isolate({
      tryCatch({
        obj <- metacoder::parse_phyloseq(ps())
        obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "otu_table")
        obj$data$tax_prop  <- obj$data$tax_abund
        obj$data$tax_prop[,-1] <- obj$data$tax_abund[,-1] /
          colSums(obj$data$tax_abund[,-1], na.rm=TRUE)
        
        metacoder::heat_tree(obj,
                             node_label  = taxon_names,
                             node_size   = n_obs,
                             node_color  = n_obs,
                             node_size_axis_label = "OTU count",
                             node_color_axis_label= "OTU count",
                             initial_layout       = "reingold-tilford",
                             layout               = "davidson-harel",
                             title = paste("Taxonomic Heat Tree –", input$taxLevelHeatTree, "level"))
      }, error = function(e) {
        plot.new()
        title(paste("Heat Tree error:", e$message))
      })
    })
  })
  
  # ============================================================
  # INTERACTIVE HEATMAP (CLR-transformed, improved)
  # ============================================================
  output$interactiveHeatmap <- renderPlotly({
    req(ps())
    input$updateHeatmap
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelHeatmap)
      
      # Select top N taxa
      top_t <- names(sort(taxa_sums(ps_glom), decreasing=TRUE))[1:min(input$hmTopN, ntaxa(ps_glom))]
      ps_top <- prune_taxa(top_t, ps_glom)
      
      # Apply transformation
      ps_trans <- tryCatch(
        microbiome::transform(ps_top, input$hmTransform),
        error = function(e) transform_sample_counts(ps_top, function(x) x/sum(x))
      )
      
      mat <- as(otu_table(ps_trans), "matrix")
      rownames(mat) <- as.character(tax_table(ps_trans)[, input$taxLevelHeatmap])
      
      # Column annotation if group var selected
      col_ann <- NULL
      grp_var <- input$regionHeatmap
      if (!is.null(grp_var) && nchar(grp_var)>0 &&
          grp_var %in% colnames(sample_data(ps()))) {
        grps <- as.character(sample_data(ps_trans)[[grp_var]])
        col_ann <- data.frame(Group=grps, row.names=colnames(mat))
      }
      
      heatmaply(t(mat),
                col_side_colors = col_ann,
                colors = viridis(256),
                main   = paste(input$hmTransform, "Heatmap –", input$taxLevelHeatmap),
                xlab   = "Taxa", ylab = "Samples",
                margins = c(80, 80, 40, 40))
    })
  })
  
  # ============================================================
  # DENDROGRAM (Bray-Curtis + Ward, standardised microbiome approach)
  # ============================================================
  output$dendrogram <- renderPlotly({
    req(ps())
    input$updateDendro
    isolate({
      otu_mat <- t(as(otu_table(ps()), "matrix"))
      # Normalise before distance calculation
      otu_norm <- vegan::decostand(otu_mat, "total")
      
      dist_mat <- if (input$distanceMethod %in% c("bray","jaccard"))
        vegan::vegdist(otu_norm, method=input$distanceMethod)
      else
        dist(otu_norm, method=input$distanceMethod)
      
      hc   <- hclust(dist_mat, method=input$clusterMethod)
      dend <- as.dendrogram(hc)
      
      # Color labels by group if selected
      grp_var <- input$dendroColorBy
      if (!is.null(grp_var) && nchar(grp_var)>0 &&
          grp_var %in% colnames(sample_data(ps()))) {
        meta_df  <- as.data.frame(sample_data(ps()))
        grp_fac  <- as.factor(meta_df[[grp_var]])
        pal      <- setNames(brewer.pal(min(nlevels(grp_fac),8),"Set1"), levels(grp_fac))
        label_colors <- pal[grp_fac[hc$order]]
        dend <- dendextend::color_labels(dend, col=label_colors)
      }
      
      # Convert to ggplot via ggdendro
      dend_data <- ggdendro::dendro_data(hc, type="rectangle")
      p <- ggplot() +
        geom_segment(data=dend_data$segments,
                     aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.5) +
        geom_text(data=dend_data$labels,
                  aes(x=x, y=y, label=label), hjust=1.1, size=3, angle=90) +
        scale_y_reverse(expand=c(0.3,0)) +
        theme_minimal() +
        labs(title=paste("Sample Dendrogram |",input$distanceMethod,"distance |",
                         input$clusterMethod,"clustering"),
             x="", y="Height")
      ggplotly(p)
    })
  })
  
  output$dendroBar <- renderPlotly({
    req(ps())
    input$updateDendro
    isolate({
      if (!input$dendroShowBar) return(NULL)
      ps_glom <- safe_tax_glom(ps(), "Phylum")
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      p <- plot_bar(ps_rel, fill="Phylum") +
        theme_minimal(base_size=11) +
        labs(title="Phylum Composition (ordered as dendrogram)")
      ggplotly(p)
    })
  })
  
  # ============================================================
  # CORRELATION NETWORK (Spearman, real correlations)
  # ============================================================
  output$correlationNetwork <- renderForceNetwork({
    req(ps())
    input$updateNetwork
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelNetwork)
      
      # Top N taxa
      n_taxa <- min(input$netTopN, ntaxa(ps_glom))
      top_t  <- names(sort(taxa_sums(ps_glom), decreasing=TRUE))[1:n_taxa]
      ps_top <- prune_taxa(top_t, ps_glom)
      
      # CLR-transform before correlation
      ps_clr <- tryCatch(clr_transform(ps_top),
                         error=function(e) compositional(ps_top))
      
      otu_mat <- as(otu_table(ps_clr), "matrix")   # taxa x samples
      # Spearman/Pearson correlation across taxa
      cor_mat <- cor(t(otu_mat), method=input$corMethod)
      
      # Build edge list
      tax_labels <- as.character(tax_table(ps_top)[, input$taxLevelNetwork])
      diag(cor_mat) <- 0
      idx <- which(abs(cor_mat) >= input$corThreshold & upper.tri(cor_mat), arr.ind=TRUE)
      
      if (nrow(idx)==0) {
        # Return empty network
        nodes <- data.frame(name=tax_labels[1], group=1, stringsAsFactors=FALSE)
        links <- data.frame(source=integer(0), target=integer(0), value=numeric(0))
        return(forceNetwork(Links=links, Nodes=nodes,
                            Source="source", Target="target",
                            NodeID="name", Group="group",
                            opacity=0.8, zoom=TRUE))
      }
      
      cor_vals <- cor_mat[idx]
      links <- data.frame(
        source = idx[,1] - 1,
        target = idx[,2] - 1,
        value  = abs(cor_vals),
        sign   = ifelse(cor_vals > 0, 1, 2)
      )
      
      if (!input$netNegEdges) links <- links[links$sign==1,]
      
      # Node abundance as group bin
      abund   <- rowSums(as.matrix(otu_table(ps_top)))
      grp_bin <- as.integer(cut(abund, breaks=5, labels=FALSE))
      nodes   <- data.frame(name=tax_labels, group=grp_bin, stringsAsFactors=FALSE)
      
      forceNetwork(Links=links, Nodes=nodes,
                   Source="source", Target="target",
                   NodeID="name", Group="group",
                   Value="value",
                   opacity=0.85, zoom=TRUE,
                   linkDistance=120, charge=-50,
                   colourScale=JS('d3.scaleOrdinal(d3.schemeCategory10)'))
    })
  })
  
  output$networkSummary <- renderPrint({
    req(ps())
    input$updateNetwork
    isolate({
      ps_glom <- safe_tax_glom(ps(), input$taxLevelNetwork)
      n_taxa  <- min(input$netTopN, ntaxa(ps_glom))
      top_t   <- names(sort(taxa_sums(ps_glom), decreasing=TRUE))[1:n_taxa]
      ps_top  <- prune_taxa(top_t, ps_glom)
      ps_clr  <- tryCatch(clr_transform(ps_top), error=function(e) compositional(ps_top))
      otu_mat <- as(otu_table(ps_clr), "matrix")
      cor_mat <- cor(t(otu_mat), method=input$corMethod)
      diag(cor_mat) <- 0
      n_pos <- sum(cor_mat >=  input$corThreshold) / 2
      n_neg <- sum(cor_mat <= -input$corThreshold) / 2
      cat("=== Correlation Network Summary ===\n")
      cat("Taxa analysed      :", n_taxa, "\n")
      cat("Correlation method :", input$corMethod, "\n")
      cat("Threshold          :", input$corThreshold, "\n")
      cat("Positive edges     :", n_pos, "\n")
      cat("Negative edges     :", n_neg, "\n")
      cat("Total edges        :", n_pos + n_neg, "\n")
    })
  })
  
  # ============================================================
  # DIFFERENTIAL ABUNDANCE (Kruskal-Wallis + BH)
  # ============================================================
  output$daVolcano <- renderPlotly({
    req(ps())
    input$updateDA
    isolate({
      grp_var <- input$regionDA
      if (is.null(grp_var) || nchar(grp_var)==0 ||
          !grp_var %in% colnames(sample_data(ps()))) {
        return(plot_ly() %>%
                 add_annotations(text="Select a grouping variable", x=0.5, y=0.5, showarrow=FALSE))
      }
      
      ps_glom <- safe_tax_glom(ps(), input$taxLevelDA)
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      melted  <- psmelt(ps_rel)
      tax_col <- input$taxLevelDA
      
      grps <- as.character(melted[[grp_var]])
      results <- melted %>%
        group_by(across(all_of(tax_col))) %>%
        summarise(
          mean_abund = mean(Abundance),
          kw_p = tryCatch(
            kruskal.test(Abundance ~ factor(.data[[grp_var]]))$p.value,
            error=function(e) NA_real_
          ),
          .groups="drop"
        ) %>%
        mutate(
          fdr    = p.adjust(kw_p, method="BH"),
          neg_log_p = -log10(pmax(fdr, 1e-10)),
          log2_mean = log2(mean_abund + 1e-10),
          sig = fdr < input$daFDR
        )
      
      plot_ly(results,
              x=~log2_mean, y=~neg_log_p,
              color=~sig,
              colors=c("grey70", wolf_pal[2]),
              text=~get(tax_col),
              type="scatter", mode="markers",
              marker=list(size=8, opacity=0.8)) %>%
        layout(
          title  = paste("Differential Abundance –", input$taxLevelDA),
          xaxis  = list(title="log2(Mean Relative Abundance)"),
          yaxis  = list(title="-log10(FDR)"),
          shapes = list(list(type="line", x0=min(results$log2_mean),
                             x1=max(results$log2_mean),
                             y0=-log10(input$daFDR), y1=-log10(input$daFDR),
                             line=list(color="red", dash="dash")))
        )
    })
  })
  
  output$daTable <- renderDT({
    req(ps())
    input$updateDA
    isolate({
      grp_var <- input$regionDA
      if (is.null(grp_var) || nchar(grp_var)==0 ||
          !grp_var %in% colnames(sample_data(ps()))) return(NULL)
      
      ps_glom <- safe_tax_glom(ps(), input$taxLevelDA)
      ps_rel  <- transform_sample_counts(ps_glom, function(x) x/sum(x))
      melted  <- psmelt(ps_rel)
      tax_col <- input$taxLevelDA
      
      results <- melted %>%
        group_by(across(all_of(tax_col))) %>%
        summarise(
          Mean_Abundance = round(mean(Abundance), 6),
          KW_pvalue      = round(tryCatch(
            kruskal.test(Abundance ~ factor(.data[[grp_var]]))$p.value,
            error=function(e) NA_real_), 4),
          .groups="drop"
        ) %>%
        mutate(FDR = round(p.adjust(KW_pvalue, method="BH"), 4),
               Significant = FDR < input$daFDR) %>%
        arrange(FDR)
      
      datatable(results,
                options=list(scrollX=TRUE, pageLength=15),
                class="cell-border stripe") %>%
        formatStyle("Significant",
                    backgroundColor=styleEqual(TRUE, "#fff3cd"))
    })
  })
  
  # ============================================================
  # TRANSFORMATION EXPLORER
  # ============================================================
  output$transPCoA <- renderPlotly({
    req(ps())
    input$updateTrans
    isolate({
      ps_use <- if (input$transMethod == "raw") ps() else
        tryCatch(microbiome::transform(ps(), input$transMethod),
                 error=function(e) ps())
      suppressMessages({
        ord <- ordinate(ps_use, method="PCoA", distance="euclidean")
      })
      p <- plot_ordination(ps_use, ord, type="samples") +
        theme_minimal(base_size=12) +
        labs(title=paste("PCoA –", input$transMethod, "transformation"))
      grp_var <- input$transColorBy
      if (!is.null(grp_var) && nchar(grp_var)>0 &&
          grp_var %in% colnames(sample_data(ps())))
        p <- p + geom_point(aes_string(color=grp_var), size=3) +
        scale_color_brewer(palette="Set1")
      else
        p <- p + geom_point(color=wolf_pal[2], size=3)
      ggplotly(p)
    })
  })
  
  output$transHist <- renderPlotly({
    req(ps())
    input$updateTrans
    isolate({
      ps_use <- if (input$transMethod == "raw") ps() else
        tryCatch(microbiome::transform(ps(), input$transMethod),
                 error=function(e) ps())
      vals <- as.vector(as(otu_table(ps_use), "matrix"))
      vals <- vals[vals > 0]
      plot_ly(x=vals, type="histogram",
              marker=list(color=wolf_pal[2], opacity=0.7),
              nbinsx=60) %>%
        layout(title=paste("Abundance Distribution –", input$transMethod),
               xaxis=list(title="Abundance value"),
               yaxis=list(title="Frequency"))
    })
  })
  
}

# =============================================================================
# RUN
# =============================================================================
shinyApp(ui = ui, server = server)



# Load the required libraries
library(shiny, lib.loc = custom_lib_path)
library(phyloseq, lib.loc = custom_lib_path)
library(plotly, lib.loc = custom_lib_path)
library(DT, lib.loc = custom_lib_path)
library(vegan, lib.loc = custom_lib_path)
library(ape, lib.loc = custom_lib_path)
library(ggtree, lib.loc = custom_lib_path)
library(ggplot2, lib.loc = custom_lib_path)
library(tidyverse, lib.loc = custom_lib_path)
library(metacoder, lib.loc = custom_lib_path)
library(networkD3, lib.loc = custom_lib_path)
library(heatmaply, lib.loc = custom_lib_path)
library(patchwork, lib.loc = custom_lib_path)



ui <- fluidPage(
  # Custom CSS for styling with wolf-inspired color scheme
  tags$head(
    tags$style(HTML("
      /* Wolf-inspired color scheme */
      :root {
        --color-wolf-dark: #2D2926;
        --color-wolf-light: #E5E4E2;
        --color-wolf-gold: #D4AF37;
        --color-wolf-warm: #8B5A2B;
        --color-wolf-cool: #6D7B8D;
      }

      body {
        font-family: 'Arial', sans-serif;
        background-color: #F5F5F5;
        background-image: url('canis1.png');
        background-size: 300px;
        background-position: right bottom;
        background-repeat: no-repeat;a
        background-attachment: fixed;
        opacity: 0.95;
      }

      .main-container {
        background-color: rgba(255, 255, 255, 0.9);
        border-radius: 8px;
        padding: 20px;
        margin-top: 20px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }

      .navbar {
        background-color: var(--color-wolf-dark) !important;
      }

      .footer {
        background-color: var(--color-wolf-dark);
        color: var(--color-wolf-light);
        padding: 20px;
        margin-top: 30px;
        text-align: center;
        font-size: 14px;
        border-top: 3px solid var(--color-wolf-gold);
      }

      .contact-links a {
        margin: 0 15px;
        color: var(--color-wolf-gold);
        text-decoration: none;
        font-weight: bold;
      }

      .contact-links a:hover {
        color: var(--color-wolf-light);
      }

      .main-title {
        margin-bottom: 5px;
        color: var(--color-wolf-dark);
      }

      .subtitle {
        color: var(--color-wolf-cool);
        font-style: italic;
        margin-bottom: 20px;
      }

      .sidebar-panel {
        background-color: white;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        padding: 15px;
        border-left: 4px solid var(--color-wolf-warm);
      }

      .main-panel {
        background-color: white;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        padding: 20px;
      }

      .btn-primary {
        background-color: var(--color-wolf-dark);
        border-color: var(--color-wolf-dark);
        color: white;
      }

      .btn-primary:hover {
        background-color: var(--color-wolf-warm);
        border-color: var(--color-wolf-warm);
        color: white;
      }

      .tabbable > .nav > li > a {
        color: var(--color-wolf-cool);
      }

      .tabbable > .nav > li.active > a {
        color: var(--color-wolf-dark);
        font-weight: bold;
        border-bottom: 2px solid var(--color-wolf-gold);
      }

      h4 {
        color: var(--color-wolf-dark);
        border-bottom: 2px solid var(--color-wolf-warm);
        padding-bottom: 5px;
      }

      hr {
        border-top: 1px solid var(--color-wolf-warm);
      }
    "))
  ),
  
  # Main title with logo and subtitle
  div(class = "main-container",
      div(
        style = "display: flex; align-items: center; margin-bottom: 20px;",
        img(src = "canis2.png", height = "80px", style = "margin-right: 20px; border: 3px solid var(--color-wolf-gold); border-radius: 50%;", alt = "Canis lupus logo"),
        div(
          h2("CanisLupus2.0", class = "main-title", style = "margin-bottom: 2;"),
          h4("Microbiome Analysis Dashboard", class = "subtitle", style = "margin-top: 10px;")
        )
      ),
      
      # Main content
      navlistPanel(
        id = "mainTabs",
        widths = c(2, 10),
        "Home",
        tabPanel("Data Upload",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     h4("Data Input", style = "color: var(--color-wolf-dark);"),
                     fileInput("asvFile", "Upload ASV Table (CSV)", accept = ".csv", width = "100%"),
                     fileInput("taxFile", "Upload Taxonomy Table (CSV)", accept = ".csv", width = "100%"),
                     fileInput("metaFile", "Upload Metadata Table (CSV)", accept = ".csv", width = "100%"),
                     fileInput("treeFile", "Upload Phylogenetic Tree", accept = c(".tree", ".tre", ".nwk", ".txt"), width = "100%"),
                     actionButton("update", "Load Data", class = "btn-primary", style = "width: 100%;")
                   ),
                   mainPanel(
                     width = 9,
                     class = "main-panel",
                     h4("Data Summary", style = "color: var(--color-wolf-dark);"),
                     verbatimTextOutput("dataSummary"),
                     hr(),
                     h4("Sample Table", style = "color: var(--color-wolf-dark);"),
                     DTOutput("sampleTable")
                   )
                 )
        ),
        "Visual Exploration",
        tabPanel("Stacked Bar Plots",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("taxLevelBar", "Taxonomic Level:",
                                 choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                 selected = "Phylum"),
                     selectInput("regionBar", "Region/Sample Group:",
                                 choices = NULL),
                     sliderInput("topNBar", "Number of Top Taxa to Show:",
                                 min = 5, max = 20, value = 10),
                     actionButton("updateBar", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("taxaBarplot", height = "600px"),
                     plotlyOutput("relativeAbundancePlot", height = "600px")
                   )
                 )
        ),
        tabPanel("Interactive Pie Chart",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("taxLevelPie", "Taxonomic Level:",
                                 choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                 selected = "Phylum"),
                     selectInput("regionPie", "Region/Sample Group:",
                                 choices = NULL),
                     actionButton("updatePie", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("pieChart", height = "800px")
                   )
                 )
        ),
        tabPanel("Rarefaction Curve",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("regionRare", "Region/Sample Group:",
                                 choices = NULL),
                     actionButton("updateRare", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("rarefactionCurve", height = "600px")
                   )
                 )
        ),
        tabPanel("Phylogenetic Tree",
                 plotOutput("phylogeneticTree", height = "800px")
        ),
        tabPanel("Heat Tree",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("taxLevelHeatTree", "Taxonomic Level:",
                                 choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                 selected = "Phylum"),
                     actionButton("updateHeatTree", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotOutput("heatTree", height = "800px")
                   )
                 )
        ),
        "Community Profile",
        tabPanel("Alpha Diversity",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("alphaMeasure", "Diversity Measure:",
                                 choices = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"),
                                 selected = "Shannon"),
                     selectInput("regionAlpha", "Group by:",
                                 choices = NULL),
                     actionButton("updateAlpha", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("alphaDiv", height = "600px"),
                     plotlyOutput("alphaDivBoxplot", height = "600px")
                   )
                 )
        ),
        tabPanel("Beta Diversity",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("betaMethod", "Method:",
                                 choices = c("PCoA", "NMDS", "CCA", "RDA"),
                                 selected = "PCoA"),
                     selectInput("betaDistance", "Distance Metric:",
                                 choices = c("bray", "jaccard", "unifrac", "wunifrac", "euclidean"),
                                 selected = "bray"),
                     selectInput("regionBeta", "Color by:",
                                 choices = NULL),
                     actionButton("updateBeta", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("betaPlot", height = "600px"),
                     plotlyOutput("nmdsPlot", height = "600px")
                   )
                 )
        ),
        tabPanel("Core Microbiome",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("taxLevelCore", "Taxonomic Level:",
                                 choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                 selected = "Genus"),
                     selectInput("regionCore", "Region/Sample Group:",
                                 choices = NULL),
                     sliderInput("prevalenceCore", "Prevalence Threshold (%):",
                                 min = 5, max = 100, value = 50),
                     sliderInput("detectionCore", "Detection Threshold (%):",
                                 min = 0.001, max = 10, value = 0.1, step = 0.01),
                     actionButton("updateCore", "Update Analysis", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("coreHeatmap", height = "600px"),
                     verbatimTextOutput("coreSummary")
                   )
                 )
        ),
        "Network Analysis",
        tabPanel("Interactive Heatmap",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("taxLevelHeatmap", "Taxonomic Level:",
                                 choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                 selected = "Genus"),
                     selectInput("regionHeatmap", "Region/Sample Group:",
                                 choices = NULL),
                     actionButton("updateHeatmap", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("interactiveHeatmap", height = "800px")
                   )
                 )
        ),
        tabPanel("Dendrogram",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("distanceMethod", "Distance Method:",
                                 choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                                 selected = "euclidean"),
                     selectInput("clusterMethod", "Clustering Method:",
                                 choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"),
                                 selected = "complete"),
                     actionButton("updateDendro", "Update Plot", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     plotlyOutput("dendrogram", height = "800px")
                   )
                 )
        ),
        tabPanel("Correlation Network",
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     class = "sidebar-panel",
                     selectInput("taxLevelNetwork", "Taxonomic Level:",
                                 choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                                 selected = "Genus"),
                     sliderInput("corThreshold", "Correlation Threshold:",
                                 min = 0.5, max = 0.99, value = 0.7, step = 0.05),
                     actionButton("updateNetwork", "Update Network", class = "btn-primary")
                   ),
                   mainPanel(
                     width = 9,
                     forceNetworkOutput("correlationNetwork", height = "800px")
                   )
                 )
        )
      )
  ),
  
  # Footer with contact information
  div(class = "footer",
      h4("Contact Information", style = "color: var(--color-wolf-gold);"),
      div(class = "contact-links",
          a(href = "mailto:pyappiah561@gmail.com", icon("envelope", style = "margin-right: 5px;"), "Email"),
          a(href = "https://www.linkedin.com/in/philip-appiah", icon("linkedin", style = "margin-right: 5px;"), "LinkedIn"),
          a(href = "https://github.com/barah123", icon("github", style = "margin-right: 5px;"), "GitHub")
      ),
      p(style = "color: var(--color-wolf-gold); margin-top: 15px;",
        "© 2025 Canis Lupus Microbiome Project.|Designed by Philip Appiah|Tel:+1 202 500 8302|All rights reserved.")
  )
)

server <- function(input, output, session) {
  ps <- reactiveVal(NULL)
  
  # Update region choices based on metadata
  observe({
    req(ps())
    meta_vars <- colnames(sample_data(ps()))
    updateSelectInput(session, "regionBar", choices = meta_vars)
    updateSelectInput(session, "regionPie", choices = meta_vars)
    updateSelectInput(session, "regionRare", choices = meta_vars)
    updateSelectInput(session, "regionAlpha", choices = meta_vars)
    updateSelectInput(session, "regionBeta", choices = meta_vars)
    updateSelectInput(session, "regionCore", choices = meta_vars)
    updateSelectInput(session, "regionHeatmap", choices = meta_vars)
  })
  
  # Load and process data
  observeEvent(input$update, {
    tryCatch({
      showNotification("Loading and processing data...", type = "message", duration = NULL)
      
      # Load ASV table
      asv <- if (!is.null(input$asvFile)) {
        read.csv(input$asvFile$datapath, row.names = 1, check.names = FALSE)
      } else {
        read.csv("demo_asv.csv", row.names = 1, check.names = FALSE)
      }
      
      # Load taxonomy table
      tax <- if (!is.null(input$taxFile)) {
        tax_table(as.matrix(read.csv(input$taxFile$datapath, row.names = 1)))
      } else {
        tax_table(as.matrix(read.csv("demo_taxonomy.csv", row.names = 1)))
      }
      
      # Load metadata
      meta <- if (!is.null(input$metaFile)) {
        sample_data(read.csv(input$metaFile$datapath, row.names = 1))
      } else {
        sample_data(read.csv("demo_metadata.csv", row.names = 1))
      }
      
      # Create phyloseq object
      otu <- otu_table(as.matrix(asv), taxa_are_rows = TRUE)
      ps_obj <- phyloseq(otu, tax, meta)
      
      # Load phylogenetic tree if provided
      if (!is.null(input$treeFile)) {
        tree <- read.tree(input$treeFile$datapath)
        phy_tree(ps_obj) <- tree
      }
      
      ps(ps_obj)
      removeNotification()
      showNotification("Data loaded successfully!", type = "message", duration = 3)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Data summary
  output$dataSummary <- renderPrint({
    req(ps())
    cat("=== Microbiome Data Summary ===\n\n")
    cat("Number of samples:", nsamples(ps()), "\n")
    cat("Number of features:", ntaxa(ps()), "\n")
    cat("\nSample variables:\n")
    print(colnames(sample_data(ps())))
    cat("\nTaxonomic ranks:\n")
    print(rank_names(ps()))
    if (!is.null(phy_tree(ps(), error = FALSE))) {
      cat("\nPhylogenetic tree tips:", length(phy_tree(ps())$tip.label), "\n")
    }
  })
  
  # Sample table
  output$sampleTable <- renderDT({
    req(ps())
    datatable(as.data.frame(sample_data(ps())),
              options = list(scrollX = TRUE, pageLength = 10),
              class = 'cell-border stripe',
              style = "bootstrap")
  })
  
  # Stacked Bar Plots
  output$taxaBarplot <- renderPlotly({
    req(ps())
    input$updateBar
    
    isolate({
      ps_rel <- transform_sample_counts(ps(), function(x) x / sum(x))
      top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE)[1:input$topNBar])
      ps_top <- prune_taxa(top_taxa, ps_rel)
      
      p <- plot_bar(ps_top, fill = input$taxLevelBar) +
        theme_minimal() +
        labs(title = paste("Top", input$topNBar, input$taxLevelBar, "Abundance")) +
        scale_fill_brewer(palette = "Set2")
      
      if (!is.null(input$regionBar) && input$regionBar %in% colnames(sample_data(ps()))) {
        p <- p + facet_wrap(as.formula(paste("~", input$regionBar)), scales = "free_x")
      }
      
      ggplotly(p)
    })
  })
  
  # Relative Abundance Plot
  output$relativeAbundancePlot <- renderPlotly({
    req(ps())
    input$updateBar
    
    isolate({
      ps_glom <- tax_glom(ps(), input$taxLevelBar)
      ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
      
      p <- plot_bar(ps_rel, fill = input$taxLevelBar) +
        theme_minimal() +
        labs(title = paste(input$taxLevelBar, "-level Relative Abundance")) +
        scale_fill_brewer(palette = "Set3")
      
      if (!is.null(input$regionBar) && input$regionBar %in% colnames(sample_data(ps()))) {
        p <- p + facet_wrap(as.formula(paste("~", input$regionBar)), scales = "free_x")
      }
      
      ggplotly(p)
    })
  })
  
  # Interactive Pie Chart
  output$pieChart <- renderPlotly({
    req(ps())
    input$updatePie
    
    isolate({
      ps_glom <- tax_glom(ps(), input$taxLevelPie)
      ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
      
      if (!is.null(input$regionPie) && input$regionPie %in% colnames(sample_data(ps()))) {
        melted <- psmelt(ps_rel)
        agg <- melted %>%
          group_by(!!sym(input$taxLevelPie), !!sym(input$regionPie)) %>%
          summarize(Abundance = mean(Abundance))
        
        plot_ly(agg, labels = ~get(input$taxLevelPie), values = ~Abundance,
                color = ~get(input$regionPie), type = 'pie',
                textposition = 'inside', textinfo = 'label+percent',
                marker = list(colors = RColorBrewer::brewer.pal(8, "Set3"))) %>%
          layout(title = paste(input$taxLevelPie, 'Composition by', input$regionPie))
      } else {
        melted <- psmelt(ps_rel)
        agg <- melted %>% group_by(!!sym(input$taxLevelPie)) %>% summarize(Abundance = mean(Abundance))
        
        plot_ly(agg, labels = ~get(input$taxLevelPie), values = ~Abundance, type = 'pie',
                textposition = 'inside', textinfo = 'label+percent',
                marker = list(colors = RColorBrewer::brewer.pal(8, "Set3"))) %>%
          layout(title = paste(input$taxLevelPie, '-level Composition'))
      }
    })
  })
  
  # Rarefaction Curve
  output$rarefactionCurve <- renderPlotly({
    req(ps())
    input$updateRare
    
    isolate({
      if (!is.null(input$regionRare) && input$regionRare %in% colnames(sample_data(ps()))) {
        # Color by selected region
        sample_data <- as(sample_data(ps()), "data.frame")
        colors <- factor(sample_data[[input$regionRare]])
        
        rare_curve <- rarecurve(t(otu_table(ps())), step = 50, label = FALSE,
                                col = colors)
        
        # Create a custom plotly version
        samples <- row.names(sample_data)
        n_samples <- length(samples)
        
        # Create a data frame for plotting
        plot_data <- map_df(1:n_samples, function(i) {
          curve <- rare_curve[[i]]
          data.frame(
            Sample = samples[i],
            Reads = attr(curve, "Subsample"),
            OTUs = curve,
            Group = sample_data[[input$regionRare]][i]
          )
        })
        
        ggplotly(
          ggplot(plot_data, aes(x = Reads, y = OTUs, group = Sample, color = Group)) +
            geom_line() +
            theme_minimal() +
            labs(title = "Rarefaction Curves",
                 x = "Number of Reads",
                 y = "Number of OTUs")
        )
      } else {
        rarecurve(t(otu_table(ps())), step = 50, label = FALSE) %>%
          ggplotly() %>%
          layout(title = "Rarefaction Curves")
      }
    })
  })
  
  # Phylogenetic Tree
  output$phylogeneticTree <- renderPlot({
    req(ps())
    if (!is.null(phy_tree(ps(), error = FALSE))) {
      ggtree(phy_tree(ps())) +
        geom_tiplab(size = 3) +
        ggtitle("Phylogenetic Tree") +
        theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # Heat Tree
  output$heatTree <- renderPlot({
    req(ps())
    input$updateHeatTree
    
    isolate({
      # Convert to metacoder object
      obj <- parse_phyloseq(ps())
      
      # Aggregate at selected taxonomic level
      obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table")
      
      heat_tree(obj,
                node_label = taxon_names,
                node_size = n_obs,
                node_color = n_obs,
                initial_layout = "reingold-tilford",
                title = paste("Taxonomic Heat Tree at", input$taxLevelHeatTree, "Level"))
    })
  })
  
  # Alpha Diversity
  output$alphaDiv <- renderPlotly({
    req(ps())
    input$updateAlpha
    
    isolate({
      measures <- input$alphaMeasure
      if (length(measures) == 0) measures <- "Shannon"
      
      p <- plot_richness(ps(), measures = measures) +
        theme_minimal() +
        labs(title = paste("Alpha Diversity (", measures, ")"))
      
      if (!is.null(input$regionAlpha) && input$regionAlpha %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color = input$regionAlpha), size = 3)
      }
      
      ggplotly(p)
    })
  })
  
  output$alphaDivBoxplot <- renderPlotly({
    req(ps())
    input$updateAlpha
    
    isolate({
      if (!is.null(input$regionAlpha) && input$regionAlpha %in% colnames(sample_data(ps()))) {
        p <- plot_richness(ps(), x = input$regionAlpha, measures = input$alphaMeasure) +
          geom_boxplot(fill = "lightblue") +
          theme_minimal() +
          labs(title = paste(input$alphaMeasure, "Diversity by", input$regionAlpha))
        
        ggplotly(p)
      }
    })
  })
  
  # Beta Diversity
  output$betaPlot <- renderPlotly({
    req(ps())
    input$updateBeta
    
    isolate({
      ord <- ordinate(ps(), method = input$betaMethod, distance = input$betaDistance)
      
      p <- plot_ordination(ps(), ord) +
        theme_minimal() +
        labs(title = paste(input$betaMethod, "Plot (", input$betaDistance, "Distance)"))
      
      if (!is.null(input$regionBeta) && input$regionBeta %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color = input$regionBeta), size = 3)
      } else {
        p <- p + geom_point(color = "blue", size = 3)
      }
      
      ggplotly(p)
    })
  })
  
  output$nmdsPlot <- renderPlotly({
    req(ps())
    input$updateBeta
    
    isolate({
      ord <- ordinate(ps(), method = "NMDS", distance = input$betaDistance)
      
      p <- plot_ordination(ps(), ord) +
        theme_minimal() +
        labs(title = paste("NMDS Plot (", input$betaDistance, "Distance)"))
      
      if (!is.null(input$regionBeta) && input$regionBeta %in% colnames(sample_data(ps()))) {
        p <- p + geom_point(aes_string(color = input$regionBeta), size = 3)
      } else {
        p <- p + geom_point(color = "blue", size = 3)
      }
      
      ggplotly(p)
    })
  })
  
  # Core Microbiome
  output$coreHeatmap <- renderPlotly({
    req(ps())
    input$updateCore
    
    isolate({
      # First aggregate at selected taxonomic level
      ps_glom <- tax_glom(ps(), input$taxLevelCore)
      
      # Calculate core microbiome
      core <- core_members(ps_glom,
                           detection = input$detectionCore/100,
                           prevalence = input$prevalenceCore/100)
      
      if (length(core) > 0) {
        ps_core <- prune_taxa(core, ps_glom)
        
        # If a region is selected, group samples by region
        if (!is.null(input$regionCore) && input$regionCore %in% colnames(sample_data(ps()))) {
          sample_data(ps_core)$Group <- sample_data(ps_core)[[input$regionCore]]
          ps_core <- merge_samples(ps_core, "Group")
          ps_core <- transform_sample_counts(ps_core, function(x) x / sum(x))
        }
        
        heatmaply(otu_table(ps_core),
                  colors = viridis::viridis(256),
                  main = paste("Core Microbiome at", input$taxLevelCore, "Level"))
      } else {
        plot_ly() %>%
          add_annotations(text = "No core taxa found with current thresholds",
                          x = 0.5, y = 0.5, showarrow = FALSE) %>%
          layout(title = paste("Core Microbiome at", input$taxLevelCore, "Level"))
      }
    })
  })
  
  output$coreSummary <- renderPrint({
    req(ps())
    input$updateCore
    
    isolate({
      ps_glom <- tax_glom(ps(), input$taxLevelCore)
      core <- core_members(ps_glom,
                           detection = input$detectionCore/100,
                           prevalence = input$prevalenceCore/100)
      
      cat("=== Core Microbiome Summary ===\n\n")
      cat("Taxonomic level:", input$taxLevelCore, "\n")
      cat("Detection threshold:", input$detectionCore, "%\n")
      cat("Prevalence threshold:", input$prevalenceCore, "%\n")
      cat("\nCore microbiome features:", length(core), "\n")
      
      if (length(core) > 0) {
        cat("\nCore taxa:\n")
        print(tax_table(ps_glom)[core, input$taxLevelCore])
      } else {
        cat("\nNo core taxa found with current thresholds\n")
      }
    })
  })
  
  # Interactive Heatmap
  output$interactiveHeatmap <- renderPlotly({
    req(ps())
    input$updateHeatmap
    
    isolate({
      # Aggregate at selected taxonomic level
      ps_glom <- tax_glom(ps(), input$taxLevelHeatmap)
      
      # If a region is selected, group samples by region
      if (!is.null(input$regionHeatmap) && input$regionHeatmap %in% colnames(sample_data(ps()))) {
        sample_data(ps_glom)$Group <- sample_data(ps_glom)[[input$regionHeatmap]]
        ps_glom <- merge_samples(ps_glom, "Group")
        ps_glom <- transform_sample_counts(ps_glom, function(x) x / sum(x))
      }
      
      heatmaply(otu_table(ps_glom),
                colors = viridis::viridis(256),
                main = paste("Interactive Heatmap at", input$taxLevelHeatmap, "Level"))
    })
  })
  
  # Dendrogram
  output$dendrogram <- renderPlotly({
    req(ps())
    input$updateDendro
    
    isolate({
      # Calculate distance matrix
      dist_matrix <- dist(t(otu_table(ps())), method = input$distanceMethod)
      
      # Perform hierarchical clustering
      hc <- hclust(dist_matrix, method = input$clusterMethod)
      
      # Convert to dendrogram
      dend <- as.dendrogram(hc)
      
      # Create plotly dendrogram
      plotly::plot_ly() %>%
        add_segments(
          x = dendextend::get_nodes_xy(dend)[,1],
          y = dendextend::get_nodes_xy(dend)[,2],
          xend = dendextend::get_nodes_xy(dend)[,3],
          yend = dendextend::get_nodes_xy(dend)[,4]
        ) %>%
        layout(
          title = paste("Sample Dendrogram (",
                        input$distanceMethod, "distance,",
                        input$clusterMethod, "clustering)"),
          xaxis = list(showticklabels = FALSE),
          yaxis = list(title = "Height")
        )
    })
  })
  
  # Correlation Network
  output$correlationNetwork <- renderForceNetwork({
    req(ps())
    input$updateNetwork
    
    isolate({
      # Aggregate at selected taxonomic level
      ps_glom <- tax_glom(ps(), input$taxLevelNetwork)
      otu_mat <- as(otu_table(ps_glom), "matrix")
      
      # Calculate correlations (simplified example)
      # In a real app, you'd want to calculate proper correlations
      n_taxa <- nrow(otu_mat)
      if (n_taxa > 20) {
        top_taxa <- names(sort(rowSums(otu_mat), decreasing = TRUE))[1:20]
        otu_mat <- otu_mat[top_taxa, ]
        n_taxa <- 20
      }
      
      # Create a simple network for demonstration
      # In reality, you'd calculate correlations between taxa
      src <- sample(0:(n_taxa-1), min(30, n_taxa*2), replace = TRUE)
      target <- sample(0:(n_taxa-1), min(30, n_taxa*2), replace = TRUE)
      
      # Filter based on correlation threshold (simulated)
      keep <- runif(length(src)) > (1 - input$corThreshold)
      src <- src[keep]
      target <- target[keep]
      
      # Remove self-links
      no_self <- src != target
      src <- src[no_self]
      target <- target[no_self]
      
      # Create nodes and links
      tax_names <- rownames(otu_mat)
      nodes <- data.frame(name = tax_names, group = 1)
      links <- data.frame(source = src, target = target, value = 1)
      
      forceNetwork(Links = links, Nodes = nodes,
                   Source = "source", Target = "target",
                   NodeID = "name", Group = "group",
                   opacity = 0.8, zoom = TRUE,
                   linkDistance = 100,
                   charge = -30)
    })
  })
}


shinyApp(ui = ui, server = server)

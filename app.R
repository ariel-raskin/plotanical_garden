# GeneTrack Studio — Phase 1 Shiny App
# A GUI for building genomic track figures with plotgardener

suppressPackageStartupMessages({
  library(shiny)
  library(plotgardener)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(TxDb.Hsapiens.UCSC.hg38.refGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(grid)
})

source("R/utils.R")
source("R/render_figure.R")

# Available TxDb assemblies (extend later)
ASSEMBLIES <- list(
  "hg38 (Human)" = "TxDb.Hsapiens.UCSC.hg38.refGene"
)

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background: #f8f9fa; font-family: 'Segoe UI', Tahoma, sans-serif; }
    .sidebar-panel { background: white; border-radius: 8px; padding: 20px;
                     box-shadow: 0 1px 3px rgba(0,0,0,0.12); }
    .main-panel { background: white; border-radius: 8px; padding: 20px;
                  box-shadow: 0 1px 3px rgba(0,0,0,0.12); }
    h3 { color: #2c3e50; margin-top: 0; }
    .track-panel { border: 1px solid #dee2e6; border-radius: 6px; padding: 12px;
                   margin-bottom: 10px; background: #fdfdfd; }
    .track-header { font-weight: 600; color: #495057; margin-bottom: 8px; }
    .btn-add-track { margin-bottom: 15px; }
    #render_btn { font-size: 16px; padding: 10px 30px; }
  "))),

  titlePanel(div(
    span("\U0001F9EC", style="font-size:28px;"),
    "GeneTrack Studio",
    style = "color: #2c3e50;"
  )),

  sidebarLayout(
    sidebarPanel(
      class = "sidebar-panel", width = 4,

      h3("\U0001F4CD Region"),
      fluidRow(
        column(4, textInput("chrom", "Chrom", value = "chr1")),
        column(4, numericInput("start", "Start", value = 109758170)),
        column(4, numericInput("end", "End", value = 109945373))
      ),
      selectInput("assembly", "Assembly", choices = names(ASSEMBLIES),
                  selected = names(ASSEMBLIES)[1]),

      hr(),
      h3("\U0001F3B5 Tracks"),
      helpText("Add data tracks to your figure. They render top-to-bottom."),

      # Dynamic track list
      div(id = "track_list", uiOutput("track_ui")),

      fluidRow(
        column(6, actionButton("add_signal", "\u2795 Signal (bigWig)",
                               class = "btn btn-outline-primary btn-sm btn-add-track")),
        column(6, actionButton("add_hic", "\u2795 Hi-C (.hic)",
                               class = "btn btn-outline-primary btn-sm btn-add-track"))
      ),
      fluidRow(
        column(6, actionButton("add_ranges", "\u2795 Ranges (BED)",
                               class = "btn btn-outline-primary btn-sm btn-add-track")),
        column(6, actionButton("add_genes", "\u2795 Gene Track",
                               class = "btn btn-outline-success btn-sm btn-add-track"))
      ),

      hr(),
      h3("\U0001F4D0 Page Settings"),
      fluidRow(
        column(4, numericInput("page_w", "Width (in)", value = 9, min = 4, max = 20, step = 0.5)),
        column(4, numericInput("page_h", "Height (in)", value = 6, min = 3, max = 20, step = 0.5)),
        column(4, numericInput("dpi", "DPI", value = 600, min = 72, max = 1200, step = 50))
      ),

      hr(),
      actionButton("render_btn", "\U0001F3A8 Render Figure",
                    class = "btn btn-primary btn-lg", width = "100%"),
      br(), br(),
      downloadButton("download_png", "\U0001F4E5 Download PNG", class = "btn btn-success", width = "100%")
    ),

    mainPanel(
      class = "main-panel", width = 8,
      h3("\U0001F5BC Preview"),
      div(
        style = "text-align: center; min-height: 400px; border: 2px dashed #dee2e6; border-radius: 8px; padding: 20px;",
        conditionalPanel(
          condition = "!output.figure_ready",
          div(style = "color: #adb5bd; padding-top: 150px; font-size: 18px;",
              "\u2190 Add tracks and click Render")
        ),
        imageOutput("figure_plot", height = "auto", width = "100%")
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # Reactive: track list as a list of lists
  tracks_rv <- reactiveVal(list())
  track_counter <- reactiveVal(0)

  # Add track helpers
  add_track <- function(type, defaults = list()) {
    n <- track_counter() + 1
    track_counter(n)
    id <- paste0("track_", n)
    tk <- c(list(id = id, type = type), defaults)
    tracks_rv(c(tracks_rv(), list(tk)))
  }

  observeEvent(input$add_signal, {
    add_track("signal", list(label = paste0("Signal ", track_counter() + 1),
                             color = "#1f77b4", height = 0.75, file = NULL))
  })
  observeEvent(input$add_hic, {
    add_track("hic", list(label = "Hi-C", color = NA, height = 1.60,
                          norm = "KR", file = NULL))
  })
  observeEvent(input$add_ranges, {
    add_track("ranges", list(label = "Ranges", color = "#2ca02c", height = 0.30, file = NULL))
  })
  observeEvent(input$add_genes, {
    add_track("genes", list(label = "Genes", color = "#0B3D91", height = 0.90))
  })

  # Render track UI panels
  output$track_ui <- renderUI({
    tks <- tracks_rv()
    if (length(tks) == 0) return(helpText("No tracks added yet."))

    lapply(seq_along(tks), function(i) {
      tk <- tks[[i]]
      id <- tk$id

      # Common elements
      remove_btn <- actionButton(paste0("rm_", id), "\u274C", class = "btn btn-sm btn-outline-danger",
                                  style = "float:right; margin-top:-5px;")

      panel_contents <- list(
        div(class = "track-header",
            paste0(toupper(tk$type), ": "), remove_btn),
        textInput(paste0("label_", id), "Label", value = tk$label),
        numericInput(paste0("height_", id), "Height (in)", value = tk$height,
                     min = 0.1, max = 5, step = 0.05)
      )

      if (tk$type %in% c("signal", "hic", "ranges")) {
        panel_contents <- c(panel_contents, list(
          textInput(paste0("file_", id), "File path", value = tk$file %||% "",
                    placeholder = "Full path to file...")
        ))
      }

      if (tk$type == "signal") {
        panel_contents <- c(panel_contents, list(
          textInput(paste0("color_", id), "Color", value = tk$color %||% "#1f77b4")
        ))
      }

      if (tk$type == "hic") {
        panel_contents <- c(panel_contents, list(
          selectInput(paste0("norm_", id), "Normalization",
                      choices = c("KR", "VC", "VC_SQRT", "NONE"),
                      selected = tk$norm %||% "KR")
        ))
      }

      if (tk$type == "genes") {
        panel_contents <- c(panel_contents, list(
          textInput(paste0("color_", id), "Exon color", value = tk$color %||% "#0B3D91"),
          numericInput(paste0("maxrows_", id), "Max rows", value = 5, min = 1, max = 20)
        ))
      }

      div(class = "track-panel", panel_contents)
    })
  })

  # Remove track observers (dynamic)
  observe({
    tks <- tracks_rv()
    lapply(tks, function(tk) {
      local({
        tid <- tk$id
        observeEvent(input[[paste0("rm_", tid)]], {
          current <- tracks_rv()
          tracks_rv(Filter(function(x) x$id != tid, current))
        }, ignoreInit = TRUE, once = TRUE)
      })
    })
  })

  # Collect current track configs from inputs
  get_track_configs <- function() {
    tks <- tracks_rv()
    lapply(tks, function(tk) {
      id <- tk$id
      cfg <- list(
        type   = tk$type,
        label  = input[[paste0("label_", id)]] %||% tk$label,
        height = input[[paste0("height_", id)]] %||% tk$height,
        file   = input[[paste0("file_", id)]] %||% tk$file,
        color  = input[[paste0("color_", id)]] %||% tk$color
      )
      if (tk$type == "hic") cfg$norm <- input[[paste0("norm_", id)]] %||% "KR"
      if (tk$type == "genes") cfg$max_rows <- input[[paste0("maxrows_", id)]] %||% 5
      cfg
    })
  }

  # TxDb (reactive, could expand for multiple assemblies)
  get_txdb <- reactive({
    TxDb.Hsapiens.UCSC.hg38.refGene
  })

  # Rendered figure path
  rendered_path <- reactiveVal(NULL)

  output$figure_ready <- reactive({ !is.null(rendered_path()) })
  outputOptions(output, "figure_ready", suspendWhenHidden = FALSE)

  # RENDER
  observeEvent(input$render_btn, {
    req(input$chrom, input$start, input$end)

    track_cfgs <- get_track_configs()
    if (length(track_cfgs) == 0) {
      showNotification("Add at least one track first!", type = "warning")
      return()
    }

    region <- list(chrom = input$chrom, start = input$start, end = input$end)
    page <- list(width = input$page_w, height = input$page_h, dpi = input$dpi)

    # Compute shared signal range for signal tracks
    signal_files <- Filter(function(tk) tk$type == "signal" && !is.null(tk$file) && nchar(tk$file) > 0,
                           track_cfgs)
    shared_range <- NULL
    if (length(signal_files) > 0) {
      tryCatch({
        shared_range <- plotgardener::calcSignalRange(
          data = lapply(signal_files, function(x) x$file),
          chrom = region$chrom, chromstart = region$start, chromend = region$end,
          assembly = "hg38")
      }, error = function(e) {
        showNotification(paste("Signal range error:", e$message), type = "error")
      })
    }
    # Inject shared range
    track_cfgs <- lapply(track_cfgs, function(tk) {
      if (tk$type == "signal" && !is.null(shared_range)) tk$range <- shared_range
      tk
    })

    withProgress(message = "Rendering figure...", {
      tryCatch({
        path <- render_figure(track_cfgs, region, page, get_txdb())
        rendered_path(path)
        showNotification("Figure rendered!", type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("Render error:", e$message), type = "error", duration = 10)
      })
    })
  })

  # Display rendered figure
  output$figure_plot <- renderImage({
    req(rendered_path())
    list(src = rendered_path(), contentType = "image/png", width = "100%",
         alt = "Rendered gene track figure")
  }, deleteFile = FALSE)

  # Download
  output$download_png <- downloadHandler(
    filename = function() {
      paste0("genetrack_", input$chrom, "_", input$start, "-", input$end, ".png")
    },
    content = function(file) {
      req(rendered_path())
      file.copy(rendered_path(), file)
    }
  )
}

shinyApp(ui, server)

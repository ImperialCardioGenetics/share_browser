library(shiny)
library(tidyverse)
library(plotly)
library(reactable)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(jsonlite)

options(warn = -1)
options(plotly.event_data.warn = FALSE)
enableBookmarking("url")

app_dir <- getwd()
load(file.path(app_dir, "share_info_2026Q1.RData"))

required_objects <- c(
  "share",
  "gene_info_filtered",
  "gene_info_all",
  "hcm_genes",
  "multi_tx_genes",
  "multi_tx_exclude",
  "share_by_gene",
  "share_display",
  "share_display_by_gene",
  "share_by_variant",
  "variant_gene_lookup",
  "gene_structures",
  "gene_regions",
  "empty_exons",
  "empty_domains",
  "empty_meta",
  "gene_clinvar_counts",
  "genomic_classifications",
  "size_map"
)
missing_objects <- setdiff(required_objects, ls())
if (length(missing_objects) > 0) {
  stop(
    "share_info_2026Q1.RData is missing prepared objects: ",
    paste(missing_objects, collapse = ", "),
    ". Run ingest_share_data.R to rebuild the data file."
  )
}

# Keep the UI constants local to the app.
extract_exons <- function(df) {
  df %>%
    select(-contains(c("pfam", "cds"))) %>%
    distinct(exon_stable_id, .keep_all = TRUE)
}

extract_domains <- function(df) {
  df %>%
    select(contains(c("gene_name", "pfam"))) %>%
    distinct(pfam_name, .keep_all = TRUE) %>%
    arrange(pfam_start_g)
}

extract_meta <- function(df) {
  df %>%
    select(gene_start_bp, gene_end_bp, strand) %>%
    distinct() %>%
    slice_head(n = 1)
}

# Classification helpers used by the UI.

add_pill <- function(value) {
  value <- as.character(value)
  color <- switch(
    value,
    "P/LP" = "#FF6666",
    "VUS"  = "#FFFF66",
    "B/LB" = "#B2FF66",
    "conflicting" = "#C0C0C0",
    "NA"   = "#C0C0C0",
    "#C0C0C0"
  )
  span(value,
       style = paste("display: inline-block; padding: 0.25em 0.85em 0.50em; border-radius: 999px;",
                     "color: black; font-weight: bold; background-color:",
                     color
       ))
}

clinvar_colors <- c(
  "P/LP" = "#FF6666",
  "VUS"  = "#FFFF66",
  "B/LB" = "#B2FF66",
  "conflicting" = "#C0C0C0",
  "NA"   = "#C0C0C0"
)

# remove empty gene
hcm_genes <- setdiff(hcm_genes, "MT-TI")

# ---------------- UI ----------------
ui <- navbarPage(
  
  title = "SHaRe Genomic Data Browser V0.2.5",
  id = "navbar",
  theme = shinytheme("flatly"),
  header = tagList(
    tags$head(
      tags$style(HTML("
        .share-github-link {
          position: fixed;
          top: 8px;
          right: 18px;
          z-index: 3000;
          color: #ffffff;
          font-size: 2.1rem;
          line-height: 44px;
          width: 90px;
          height: 90px;
          text-align: center;
          text-decoration: none;
        }
        .share-github-link:hover,
        .share-github-link:focus {
          color: #dce4ec;
          text-decoration: none;
          outline: none;
        }
        @media (max-width: 800px) {
          .share-github-link {
            right: 58px;
          }
        }
      ")),
      tags$script(HTML("
      (function() {
        var lastNotifiedPath = null;

        function currentSharePath() {
          return window.location.hash || window.location.pathname || '';
        }

        function notifySharePath() {
          if (!window.Shiny) return;
          var path = currentSharePath();
          if (path === lastNotifiedPath) return;
          lastNotifiedPath = path;
          Shiny.setInputValue('deep_link_path', {path: path, nonce: Date.now()}, {priority: 'event'});
        }

        Shiny.addCustomMessageHandler('shareSetPath', function(msg) {
          if (!msg || !msg.path) return;
          var method = (msg.mode === 'push') ? 'pushState' : 'replaceState';
          var base = window.location.pathname;
          if (base.length === 0) base = '/';
          var path = msg.path;
          if (path.charAt(0) !== '/') path = '/' + path;
          window.history[method]({}, '', base + '#' + path);
          lastNotifiedPath = currentSharePath();
        });
        window.addEventListener('popstate', notifySharePath);
        window.addEventListener('hashchange', notifySharePath);
      })();
      "))
    ),
    tags$a(
      id = "share-github-link",
      class = "share-github-link",
      href = "https://github.com/ImperialCardioGenetics/share_browser",
      target = "_blank",
      rel = "noopener noreferrer",
      title = "Open GitHub repository",
      `aria-label` = "Open GitHub repository",
      icon("github")
    )
  ),
  tabPanel("Home",
        fluidPage(
                  useShinyjs(),
                  tags$style(HTML("
                                  /* Label text color */
                                  .pretty input[value='P/LP'] ~ .state span { color: #FF6666; font-weight: bold; }
                                  .pretty input[value='VUS']  ~ .state span { color: #FFFF66; font-weight: bold; }
                                  .pretty input[value='B/LB'] ~ .state span { color: #B2FF66; font-weight: bold; }
                            
                                  /* Checkbox background color */
                                  .pretty input[value='P/LP'] ~ .state label:after,
                                  .pretty input[value='P/LP'] ~ .state label:before {
                                    background-color: #FF6666;
                                  }
                                  .pretty input[value='VUS'] ~ .state label:after,
                                  .pretty input[value='VUS'] ~ .state label:before {
                                    background-color: #FFFF66;
                                  }
                                  .pretty input[value='B/LB'] ~ .state label:after,
                                  .pretty input[value='B/LB'] ~ .state label:before {
                                    background-color: #B2FF66;
                                  }
                                  .home-page {
                                    min-height: 80vh;
                                    display: flex;
                                    flex-direction: column;
                                  }
                                  .home-logos {
                                    margin-top: auto;
                                    width: 100%;
                                    display: flex;
                                    justify-content: center;
                                    align-items: flex-start;
                                    gap: 2rem;
                                    flex-wrap: wrap;
                                    padding: 1.5rem 0 2rem;
                                  }
                                  .home-logos a {
                                    display: flex;
                                    width: 220px;
                                    height: 120px;
                                    justify-content: center;
                                    align-items: center;
                                  }
                                  .home-logos img {
                                    max-width: 100%;
                                    max-height: 100%;
                                    width: auto;
                                    height: auto;
                                    object-fit: contain;
                                  }
                                  ")),
           div(class = "home-page",
           br(),
           br(),
           h1("SHaRe Genomic Data Browser", align = "center"),
           br(),br(),
             column(width = 4, offset = 4, align = "center",
                    wellPanel(h4("Select an HCM gene and view it in the Gene View tab."),
                              selectizeInput("gene",
                                             selected = "",
                                             width = 300,
                                             label = "", 
                                             choices = c("Search" = "",hcm_genes),
                                             multiple = FALSE,
                                             options = list(
                                               onChange = I("function(value) {
                                                 if (value) {
                                                   var $tab = $('a[data-value=\"Gene View\"]');
                                                   if ($tab.length) { $tab.tab('show'); }
                                                   Shiny.setInputValue('start', value, {priority: 'event'});
                                                 }
                                               }")
                                             )
                                             ),
                              
                              ),
                              style =  "background-color: #ffffff;
                                        border-bottom-color: #333333;
                                        border-left-color: #333333;
                                        border-right-color: #333333;
                                        margin-top: 0px"),
             div(class = "home-logos",
                 tags$a(href='https://www.imperial.ac.uk/', target="_blank",
                        tags$img(src='https://upload.wikimedia.org/wikipedia/commons/thumb/0/06/Imperial_College_London_new_logo.png/500px-Imperial_College_London_new_logo.png',height='50',width='100' )),
                 tags$a(href='https://www.theshareregistry.org/', target="_blank",
                        tags$img(src='https://www.theshareregistry.org/wp-content/uploads/2024/05/SHARE_Main_Logo_215px.png',height='300',width='300' )),
                 tags$a(href='https://lms.mrc.ac.uk/', target="_blank",
                        tags$img(src='https://lms.mrc.ac.uk/wp-content/uploads/logo.svg',height='150',width='400' )),
                 tags$a(href='https://www.garvan.org.au/', target="_blank",
                        tags$img(src='https://images.contentstack.io/v3/assets/blt324fd0a04af716e6/blt3f0048229394c515/6405de96205f2b7a60b745d6/gimr-logo.png',height='150',width='200' ))
             )
           
           )),
        div(
          style = "text-align:center; margin: 0.5rem 0 1.5rem 0; font-size: 1rem;",
          span("All rights reserved. Please report any issues to the web administrator "),
          tags$a(href = "mailto:p.theotokis@imperial.ac.uk", "here")
        )
        ),
        


  tabPanel("Gene View",
           fluidPage(
             shinycssloaders::withSpinner(uiOutput("gene_view_ui"), type = 4, color = "#2c3e50")
           )
  ),
  tabPanel("Variant View",
           fluidPage(
             shinycssloaders::withSpinner(uiOutput("variant_view_ui"), type = 4, color = "#2c3e50")
           )
  )
)

# ---------------- SERVER ----------------
server <- function(input, output, session) {

  create_spaced_buttons <- function(..., spacing = "8px") {
    btns <- list(...)
    tags$div(
      style = paste0("display:flex; align-items:center; gap:", spacing, "; flex-wrap:wrap;"),
      btns
    )
  }

  plain_action_button <- function(label, onclick, class = "btn btn-default btn-xs") {
    tags$button(
      type = "button",
      class = class,
      onclick = onclick,
      label
    )
  }

  modal_action_js <- function(action, variant_id) {
    payload <- jsonlite::toJSON(list(action = action, variant = variant_id), auto_unbox = TRUE)
    sprintf("var payload = %s; payload.nonce = Date.now(); Shiny.setInputValue('modal_action', payload, {priority: 'event'});", payload)
  }

  add_size_legend <- function(plot_obj) {
    legend_sizes <- c(
      #"SHaRe AC 0" = 2,
      "SHaRe AC 1-5" = 4,
      "SHaRe AC 6-10" = 6,
      "SHaRe AC 11-20" = 9,
      "SHaRe AC 21-50" = 12,
      "SHaRe AC >50" = 20
    )
    legend_labels <- names(legend_sizes)
    for (i in seq_along(legend_sizes)) {
      plot_obj <- plot_obj %>%
        add_trace(
          x = 0,
          y = 0,
          type = "scatter",
          mode = "markers",
            # marker = list(symbol = "diamond-tall", size = legend_sizes[[i]], color = "white", width = 1.5, line = "black"),
          marker = list(
            symbol = "diamond-tall",
            color = "rgba(0,0,0,0)",
            line = list(color = "black", width = 1.5),
            size = legend_sizes[[legend_labels[[i]]]]
          ),
          name = legend_labels[[i]],
          showlegend = TRUE,
          hoverinfo = "skip",
          visible = "legendonly",
          inherit = FALSE
        )
    }
    plot_obj %>%
      layout(legend = list(itemsizing = "trace"))
  }

  initial_load <- reactiveVal(TRUE)
  suppress_url_update <- reactiveVal(FALSE)

  selected_row <- reactive({
    getReactableState("variant_table", "selected")
  })

  selected_structures <- reactive({
    req(has_active_gene())
    res <- gene_structures[[active_gene()]]
    if (is.null(res)) {
      list(exons = empty_exons, domains = empty_domains, meta = empty_meta)
    } else {
      if (is.null(res$meta) || nrow(res$meta) == 0) {
        res$meta <- empty_meta
      }
      if (is.null(res$exons)) res$exons <- empty_exons
      if (is.null(res$domains)) res$domains <- empty_domains
      res
    }
  }) %>% bindCache(active_gene())
  
  
  sel_var <- reactiveVal(NULL)
  gene_selected <- reactiveVal(NULL)
  gene_loading <- reactiveVal(FALSE)
  initial_url_gene <- reactiveVal(NULL)
  pending_variant <- reactiveVal(NULL)
  pending_tab <- reactiveVal(NULL)
  invalid_variant_request <- reactiveVal(NULL)
  invalid_gene_request <- reactiveVal(NULL)

  get_variant_row <- function(variant_id) {
    if (is.null(variant_id) || !nzchar(variant_id)) return(share[0, ])
    if (!(variant_id %in% names(share_by_variant))) return(share[0, ])
    row <- share_by_variant[[variant_id]]
    if (is.null(row)) share[0, ] else row
  }

  get_variant_gene <- function(variant_id) {
    if (is.null(variant_id) || !nzchar(variant_id)) return(NULL)
    if (!(variant_id %in% names(variant_gene_lookup))) return(NULL)
    gene <- unname(variant_gene_lookup[variant_id])
    if (is.null(gene) || is.na(gene) || !nzchar(gene)) NULL else gene
  }

  select_variant <- function(variant_id, navigate = TRUE) {
    if (is.null(variant_id) || !nzchar(variant_id)) {
      pending_variant(NULL)
      return(FALSE)
    }
    row <- get_variant_row(variant_id)
    if (nrow(row) == 0) {
      pending_variant(NULL)
      return(FALSE)
    }
    target_gene <- get_variant_gene(variant_id)
    if (is.null(target_gene) || !nzchar(target_gene) || !(target_gene %in% hcm_genes)) return(FALSE)

    if (!identical(input$gene, target_gene)) {
      pending_variant(variant_id)
      updateSelectizeInput(session, "gene", selected = target_gene)
    } else {
      pending_variant(NULL)
      current_sel <- isolate(sel_var())
      if (identical(current_sel, variant_id)) {
        sel_var(NULL)
      }
      sel_var(variant_id)
    }

    if (isTRUE(navigate)) {
      updateTabsetPanel(session, inputId = "navbar", selected = "Variant View")
    }
    TRUE
  }


  open_gene_view <- function(variant_id) {
    if (is.null(variant_id) || !nzchar(variant_id)) return(FALSE)
    row <- get_variant_row(variant_id)
    if (nrow(row) == 0) return(FALSE)
    target_gene <- get_variant_gene(variant_id)
    if (is.null(target_gene) || !nzchar(target_gene)) return(FALSE)

    if (!identical(input$gene, target_gene)) {
      pending_variant(variant_id)
      updateSelectizeInput(session, "gene", selected = target_gene)
    } else {
      pending_variant(NULL)
      sel_var(NULL)
      sel_var(variant_id)
    }

    updateTabsetPanel(session, inputId = "navbar", selected = "Gene View")
    TRUE
  }

  selected_variant <- reactive({
    variant_id <- sel_var()
    if (is.null(variant_id) || !nzchar(variant_id)) return(NULL)
    row <- get_variant_row(variant_id)
    if (nrow(row) == 0) return(NULL)
    row
  }) %>% bindCache(sel_var())

  active_gene <- reactive({
    g <- gene_selected()
    if (is.null(g) || !nzchar(g)) g <- input$gene
    if (is.null(g) || !nzchar(g)) g <- initial_url_gene()
    g
  })

  has_active_gene <- function() {
    g <- active_gene()
    !is.null(g) && nzchar(g)
  }

  infer_gene_from_variant_id <- function(variant_id) {
    parts <- strsplit(variant_id, "-", fixed = TRUE)[[1]]
    if (length(parts) < 2) return(NULL)

    pos <- suppressWarnings(as.numeric(parts[2]))
    if (is.na(pos)) return(NULL)

    chr <- parts[1]
    chr <- if (startsWith(chr, "chr")) chr else paste0("chr", chr)

    hit <- gene_regions %>%
      filter(hg38_chr == chr, gene_start_bp <= pos, gene_end_bp >= pos) %>%
      mutate(span = gene_end_bp - gene_start_bp) %>%
      arrange(span)

    if (nrow(hit) == 0) return(NULL)
    hit$gene_name[1]
  }

  parse_initial_location <- function(pathname = "") {
    pathname <- pathname %||% ""
    pathname <- sub("^#/?", "/", pathname)
    parts <- strsplit(pathname, "/", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]

    tab <- gene <- variant <- NULL
    idx <- which(parts %in% c("gene", "variant"))[1]
    if (!is.na(idx)) {
      tab <- parts[idx]
      if (idx < length(parts)) {
        value <- URLdecode(paste(parts[(idx + 1):length(parts)], collapse = "/"))
        if (nzchar(value)) {
          if (identical(tab, "gene")) {
            gene <- value
          } else if (identical(tab, "variant")) {
            variant <- value
          }
        }
      }
    }

    list(tab = tab, gene = gene, variant = variant)
  }

  apply_initial_state <- function(tab = NULL, gene = NULL, variant = NULL) {
    state_variant <- variant
    state_gene <- gene
    invalid_gene <- FALSE

    if (!is.null(state_variant) && nzchar(state_variant)) {
      variant_row <- get_variant_row(state_variant)
      if (nrow(variant_row) == 0) {
        invalid_variant_request(list(
          variant = state_variant,
          gene = infer_gene_from_variant_id(state_variant),
          nonce = as.numeric(Sys.time())
        ))
        state_gene <- NULL
        state_variant <- NULL
      } else {
        state_gene <- get_variant_gene(state_variant)
      }
    }

    if (!is.null(state_gene) && nzchar(state_gene) && !(state_gene %in% hcm_genes)) {
      invalid_gene <- TRUE
      invalid_gene_request(list(
        gene = state_gene,
        tab = tab,
        nonce = as.numeric(Sys.time())
      ))
      initial_url_gene(NULL)
      state_gene <- NULL
      state_variant <- NULL
    }

    target_tab <- NULL
    if (!is.null(tab) && nzchar(tab) && !is.null(state_variant)) {
      target_tab <- slug_to_tab(tab)
    } else if (!is.null(tab) && nzchar(tab) && !is.null(state_gene)) {
      target_tab <- "Gene View"
    } else if (!is.null(state_variant)) {
      target_tab <- 'Variant View'
    } else if (!is.null(state_gene)) {
      target_tab <- 'Gene View'
    } else if (invalid_gene && identical(tab, "gene")) {
      target_tab <- 'Gene View'
    }

    if (!is.null(target_tab)) {
      pending_tab(target_tab)
    }
    if (!is.null(state_variant)) {
      pending_variant(state_variant)
    }
    if (!is.null(state_gene) && nzchar(state_gene) && state_gene %in% hcm_genes) {
      gene_loading(TRUE)
      initial_url_gene(state_gene)
      gene_selected(state_gene)
      updateSelectizeInput(session, 'gene', selected = state_gene)
    } else if (!is.null(tab) && identical(tab, "gene")) {
      initial_url_gene(NULL)
    }
    if (!is.null(state_variant)) {
      sel_var(state_variant)
    }
    if (!is.null(target_tab)) {
      updateTabsetPanel(session, inputId = 'navbar', selected = target_tab)
    }

    current_gene <- isolate(input$gene)
    if (!is.null(state_variant) && identical(current_gene, state_gene)) {
      sel_var(state_variant)
      pending_variant(NULL)
    }

    list(tab = target_tab, gene = state_gene, variant = state_variant)
  }

  tab_to_slug <- function(tab) {
    switch(tab,
           "Home" = "home",
           "Gene View" = "gene",
           "Variant View" = "variant",
           tolower(gsub("\\s+", "-", as.character(tab))))
  }

  slug_to_tab <- function(slug) {
    switch(tolower(as.character(slug)),
           "home" = "Home",
           "gene" = "Gene View",
           "variant" = "Variant View",
           "Home")
  }

  build_deep_link_path <- function(tab = "Home", gene = NULL, variant = NULL) {
    tab <- tab %||% "Home"
    if (identical(tab, "Gene View")) {
      if (is.null(gene) || !nzchar(gene)) return("/")
      return(paste0("/gene/", utils::URLencode(gene, reserved = TRUE)))
    }
    if (identical(tab, "Variant View")) {
      if (is.null(variant) || !nzchar(variant)) return("/")
      return(paste0("/variant/", utils::URLencode(variant, reserved = TRUE)))
    }
    "/"
  }

  get_initial_deep_link <- function() {
    loc_hash <- isolate(session$clientData$url_hash) %||% ""
    loc_path <- isolate(session$clientData$url_pathname) %||% ""
    if (nzchar(loc_hash)) loc_hash else loc_path
  }

  update_url_path <- function(mode = "replace") {
    if (isTRUE(initial_load())) return(invisible(NULL))
    if (isTRUE(suppress_url_update())) return(invisible(NULL))

    current_tab <- input$navbar
    if (is.null(current_tab) || !nzchar(current_tab)) current_tab <- "Home"

    path <- build_deep_link_path(
      tab = current_tab,
      gene = active_gene(),
      variant = sel_var()
    )

    session$sendCustomMessage("shareSetPath", list(path = path, mode = mode))
    invisible(NULL)
  }



  
  observeEvent(input$gene, {
    if (isTRUE(initial_load())) return()
    if (!is.null(pending_tab())) return()
    if (!is.null(pending_variant())) return()
    if (!is.null(input$gene) && nzchar(input$gene)) {
      gene_selected(input$gene)
      pending <- pending_variant()
      if (!is.null(pending)) {
        df_gene <- share_by_gene[[input$gene]]
        if (!is.null(df_gene) && pending %in% df_gene$VariantID) {
          pending_variant(NULL)
          if (!identical(isolate(sel_var()), pending)) {
            sel_var(pending)
          }
        }
      }
      if (!identical(input$navbar, "Gene View")) {
        updateTabsetPanel(session, inputId = "navbar", selected = "Gene View")
      }
      update_url_path(mode = "push")
    } else {
      if (!identical(input$navbar, "Home")) return()
      gene_selected(NULL)
      pending_variant(NULL)
      if (!is.null(sel_var())) sel_var(NULL)
      update_url_path()
    }
  }, ignoreNULL = FALSE)

  observeEvent(pending_variant(), {
    variant_id <- pending_variant()
    if (is.null(variant_id) || !nzchar(variant_id)) return()
    row <- get_variant_row(variant_id)
    if (nrow(row) == 0) {
      pending_variant(NULL)
      return()
    }
    target_gene <- get_variant_gene(variant_id)
    if (is.null(target_gene) || !nzchar(target_gene)) {
      pending_variant(NULL)
      return()
    }
    if (!identical(input$gene, target_gene)) return()

    current_sel <- isolate(sel_var())
    if (!identical(current_sel, variant_id)) {
      sel_var(variant_id)
    }
    pending_variant(NULL)
    if (identical(pending_tab(), "Variant View")) {
      updateTabsetPanel(session, inputId = "navbar", selected = "Variant View")
    }
  }, ignoreNULL = TRUE)

  observeEvent(input$navbar, {
    if (isTRUE(initial_load())) return()
    if (identical(input$navbar, "Home")) {
      pending_variant(NULL)
      sel_var(NULL)
    }
    if (!is.null(pending_tab()) && identical(input$navbar, pending_tab())) {
      pending_tab(NULL)
    }
    update_url_path()
  }, ignoreNULL = TRUE)

  observeEvent(input$start, {
    g <- input$start
    gene_loading(TRUE)
    if (!identical(input$navbar, "Gene View")) {
      updateTabsetPanel(session, inputId = "navbar", selected = "Gene View")
    }
    session$onFlushed(function() {
      if (!is.null(g) && nzchar(g) && !identical(isolate(gene_selected()), g)) {
        gene_selected(g)
      }
    }, once = TRUE)
  }, ignoreNULL = TRUE)

  observeEvent(gene_selected(), {
    g <- isolate(gene_selected())
    if (!is.null(g) && nzchar(g)) {
      gene_loading(FALSE)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$navbar, {
    if (identical(input$navbar, "Home")) {
      gene_loading(FALSE)
    }
  }, ignoreNULL = FALSE)

  observeEvent(invalid_variant_request(), {
    info <- invalid_variant_request()
    req(!is.null(info), !is.null(info$variant))

    has_gene <- !is.null(info$gene) && nzchar(info$gene)
    showModal(modalDialog(
      title = "Variant Not Found",
      p(sprintf("The requested variant '%s' was not found in this release.", info$variant)),
      if (has_gene) {
        p(sprintf("The requested position overlaps the %s gene region. Do you want to open that gene page instead?", info$gene))
      } else {
        p("No overlapping gene region could be identified from the requested variant position.")
      },
      easyClose = TRUE,
      footer = tagList(
        if (has_gene) {
          actionButton("invalid_variant_redirect_btn", sprintf("Open %s Gene Page", info$gene))
        },
        modalButton(if (has_gene) "Stay Here" else "Close")
      )
    ))
  }, ignoreNULL = TRUE)

  observeEvent(invalid_gene_request(), {
    info <- invalid_gene_request()
    req(!is.null(info), !is.null(info$gene), nzchar(info$gene))

    showModal(modalDialog(
      title = "Gene Not Available",
      p(sprintf("The requested gene '%s' is not available in the SHaRe browser for this release.", info$gene)),
      p("Use the Home tab to choose one of the supported HCM genes."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  }, ignoreNULL = TRUE)

  observeEvent(input$invalid_variant_redirect_btn, {
    info <- invalid_variant_request()
    req(!is.null(info), !is.null(info$gene), nzchar(info$gene))
    removeModal()
    invalid_variant_request(NULL)
    pending_variant(NULL)
    sel_var(NULL)
    gene_selected(info$gene)
    updateSelectizeInput(session, "gene", selected = info$gene)
    updateTabsetPanel(session, inputId = "navbar", selected = "Gene View")
    session$sendCustomMessage(
      "shareSetPath",
      list(path = build_deep_link_path(tab = "Gene View", gene = info$gene), mode = "push")
    )
  })

  observeEvent(input$gnomad_info_btn, {
    showModal(modalDialog(
      title = "gnomAD Info",
      HTML(paste0(
        "<p>This annotation can be used for filtering variants by allele frequency against a disease-specific threshold that can be set for each disease 
        (e.g. BA1 in the 2015 ACMG/AMP guidelines). 
        In this case the filtering allele frequency (FAF) is the maximum credible genetic ancestry group AF (e.g. the lower bound of the 95% confidence interval (CI)). 
        If the FAF is above the disease-specific threshold, then the observed AC is not compatible with pathogenicity. 
        See <a href='http://cardiodb.org/allelefrequencyapp/'>http://cardiodb.org/allelefrequencyapp/</a> and Whiffin et al. 2017 for additional information.</p>",
        
        "<p>Note that the GroupMax FAF contains filtering allele frequency information from the genetic ancestry group with the highest FAF, 
        not the filtering allele frequency information calculated on the genetic ancestry group with the highest AF.</p>",
        
        "<p>The filtering allele frequency calculation only includes non-bottlenecked genetic ancestry groups: 
        we did not calculate this metric on the Amish (ami), Ashkenazi Jewish (asj), European Finnish (fin), and 'Remaining Individuals' (rmi) groups. 
        Note that the exome FAF and joint (combined exome and genome) FAF calculations included the Middle Eastern (mid) group. 
        However, due to small group size, the genome FAF calculations did not include mid.</p>"
      )),
      easyClose = TRUE,
      size = "l",
      footer = modalButton("Close")
    ))
  })

  
  session$onFlushed(function() {
    if (!isTRUE(isolate(initial_load()))) return()
    loc <- parse_initial_location(get_initial_deep_link())

    apply_initial_state(tab = loc$tab, gene = loc$gene, variant = loc$variant)

    initial_load(FALSE)
  }, once = TRUE)

  observeEvent(input$deep_link_path, {
    if (isTRUE(initial_load())) return()

    msg <- input$deep_link_path
    path <- if (is.list(msg)) msg$path else msg
    path <- path %||% ""
    loc <- parse_initial_location(path)

    suppress_url_update(TRUE)

    if (is.null(loc$tab)) {
      pending_tab(NULL)
      pending_variant(NULL)
      initial_url_gene(NULL)
      invalid_variant_request(NULL)
      invalid_gene_request(NULL)
      gene_selected(NULL)
      sel_var(NULL)
      gene_loading(FALSE)
      updateSelectizeInput(session, "gene", selected = "")
      updateTabsetPanel(session, inputId = "navbar", selected = "Home")
    } else {
      apply_initial_state(tab = loc$tab, gene = loc$gene, variant = loc$variant)
    }

    session$onFlushed(function() {
      suppress_url_update(FALSE)
    }, once = TRUE)
  }, ignoreNULL = TRUE)

  
  observe({
    idx <- selected_row()
    df <- gene_df()
    if (!is.null(idx) && length(idx) == 1) {
      sv <- df$VariantID[idx]
      sel_var(sv)
    }
  })
  
  observe({
    df <- gene_df()
    sv <- sel_var()
    if (!is.null(sv) && !(sv %in% df$VariantID)) {
      sel_var(NULL)
    }
  })

  gene_df <- reactive({
    req(has_active_gene())
    classes <- input$select_class
    if (is.null(classes) || length(classes) == 0) classes <- genomic_classifications
    df <- share_display_by_gene[[active_gene()]]
    if (is.null(df)) df <- share_display[0, ]

    df %>%
      filter(VarClass %in% classes)
  }) %>% bindCache(active_gene(), input$select_class)

  build_multi_tx_gene_tracks <- function(plot_obj, tx_source, tx_ids, tx_y, pfam_colors = NULL, legend_seen = character(0)) {
    add_segment_trace <- function(plot_obj, data, x_col, xend_col, y_value, yend_value = y_value, line = list(), showlegend = FALSE, name = NULL, group_col = NULL) {
      if (is.null(data) || nrow(data) == 0) return(plot_obj)
      groups <- if (!is.null(group_col) && group_col %in% names(data)) split(data, data[[group_col]]) else list(default = data)
      for (grp_name in names(groups)) {
        grp <- groups[[grp_name]]
        if (nrow(grp) == 0) next
        show_this <- showlegend
        if (!is.null(group_col)) {
          show_this <- showlegend && !(grp_name %in% legend_seen)
          if (show_this) legend_seen <<- c(legend_seen, grp_name)
          if (!is.null(pfam_colors) && grp_name %in% names(pfam_colors)) {
            line$color <- pfam_colors[[grp_name]]
          }
        }
        x_vals <- as.vector(rbind(grp[[x_col]], grp[[xend_col]], NA))
        y_vals <- rep(c(y_value, yend_value, NA), nrow(grp))
        plot_obj <- plot_obj %>%
          add_trace(
            data = NULL,
            x = x_vals,
            y = y_vals,
            type = "scatter",
            mode = "lines",
            line = line,
            hoverinfo = "skip",
            showlegend = show_this,
            name = if (!is.null(group_col)) grp_name else name,
            legendgroup = if (!is.null(group_col)) grp_name else NULL,
            inherit = FALSE
          )
      }
      plot_obj
    }

    for (tx in tx_ids) {
      tx_df <- tx_source %>%
        filter(as.character(transcript_stable_id_version) == tx)

      tx_exons <- extract_exons(tx_df)
      if (nrow(tx_exons) > 0) {
        plot_obj <- add_segment_trace(
          plot_obj, tx_exons,
          x_col = "exon_region_start_bp", xend_col = "exon_region_end_bp",
          y_value = tx_y[[tx]],
          line = list(color = "steelblue", width = 26),
          showlegend = FALSE,
          name = tx
        )
      }

      tx_domains <- extract_domains(tx_df)
      if (nrow(tx_domains) > 0) {
        plot_obj <- add_segment_trace(
          plot_obj, tx_domains,
          x_col = "pfam_start_g", xend_col = "pfam_end_g",
          y_value = tx_y[[tx]] - 0.18,
          line = list(width = 10),
          showlegend = TRUE,
          group_col = "pfam_name"
        )
      }
    }
    plot_obj
  }
  
  
  output$gene_plot <- renderPlotly({
    suppressMessages(suppressWarnings({
      classes <- input$select_class
      if (is.null(classes) || length(classes) == 0) classes <- genomic_classifications
      req(has_active_gene())

      df <- gene_df()
    structures <- selected_structures()
    exon <- structures$exons
    dom <- structures$domains
    meta <- structures$meta
    multi_tx_gene <- active_gene() %in% multi_tx_genes
    plot_source <- NULL
    main_tx <- NULL
    extra_tx_ids <- character(0)
    extra_tx_y <- numeric(0)
    legend_seen <- character(0)
    pfam_colors <- NULL
    if (multi_tx_gene) {
      plot_source <- gene_info_all %>%
        filter(gene_name == active_gene()) %>%
        mutate(transcript_stable_id_version = as.character(transcript_stable_id_version))
      main_tx <- gene_info_filtered %>%
        filter(gene_name == active_gene()) %>%
        mutate(transcript_stable_id_version = as.character(transcript_stable_id_version)) %>%
        pull(transcript_stable_id_version) %>%
        first()
      tx_ids <- plot_source %>%
        filter(!is.na(transcript_stable_id_version), nzchar(transcript_stable_id_version)) %>%
        distinct(transcript_stable_id_version) %>%
        pull(transcript_stable_id_version)
      tx_ids <- unique(c(main_tx, setdiff(tx_ids, main_tx)))
      tx_ids <- tx_ids[!is.na(tx_ids) & nzchar(tx_ids)]
      tx_drop <- multi_tx_exclude[[active_gene()]]
      if (!is.null(tx_drop) && length(tx_drop) > 0) {
        tx_ids <- setdiff(tx_ids, tx_drop)
      }
      extra_tx_ids <- setdiff(tx_ids, main_tx)
      extra_tx_y <- setNames(seq(-1, by = -1, length.out = length(extra_tx_ids)), extra_tx_ids)
    }
    pfam_names <- unique(c(as.character(dom$pfam_name), if (!is.null(plot_source)) as.character(plot_source$pfam_name) else character(0)))
    pfam_names <- pfam_names[!is.na(pfam_names) & nzchar(pfam_names)]
    if (length(pfam_names) > 0) {
      pfam_colors <- setNames(grDevices::hcl.colors(length(pfam_names), palette = "Dark 2"), pfam_names)
    }

    add_direction_markers <- function(plot_obj, meta_df) {
      if (is.null(meta_df) || nrow(meta_df) == 0) {
        return(add_size_legend(plot_obj))
      }
      info <- meta_df[1, , drop = FALSE]
      if (any(is.na(info[c("gene_start_bp", "gene_end_bp")])) || is.na(info$strand[1])) {
        return(add_size_legend(plot_obj))
      }

      meta_df <- tibble::tibble(
        POS = c(info$gene_start_bp[1], info$gene_end_bp[1]),
        ypos = c(0, 0),
        symbol = switch(as.character(info$strand[1]),
                        "1" = "arrow-bar-right",
                        "-1" = "arrow-bar-left",
                        "triangle-up")
      )

      plot_obj <- plot_obj %>%
        add_trace(
          data = meta_df,
          x = ~POS,
          y = ~ypos,
          type = "scatter",
          mode = "markers",
          marker = list(symbol = ~symbol, color = "black", line = list(width = 1.5, color = "black"), size = 10),
          hoverinfo = "skip",
          name = "Gene direction",
          showlegend = FALSE,
          inherit = FALSE
        ) %>%
        add_trace(
          data = meta_df,
          x = ~POS,
          y = ~ypos,
          type = "scatter",
          mode = "lines",
          line = list(color = "black", width = 1.5),
          hoverinfo = "skip",
          name = "Gene span",
          showlegend = FALSE,
          inherit = FALSE
        )

      add_size_legend(plot_obj)
    }

    if (nrow(df) == 0) {
      p_empty <- plot_ly(source = "gene_plot")
      if (nrow(exon) > 0) {
        p_empty <- p_empty %>%
          add_trace(
            data = NULL,
            x = as.vector(rbind(exon$exon_region_start_bp, exon$exon_region_end_bp, NA)),
            y = rep(c(0, 0, NA), nrow(exon)),
            type = "scatter",
            mode = "lines",
            line = list(color = "steelblue", width = 30),
            hoverinfo = "skip",
            showlegend = FALSE,
            name = "Exons",
            inherit = FALSE
          )
      }
      if (nrow(dom) > 0) {
        for (pfam_name in unique(dom$pfam_name)) {
          dom_i <- dom %>% filter(pfam_name == !!pfam_name)
          legend_this <- !(pfam_name %in% legend_seen)
          if (legend_this) legend_seen <- c(legend_seen, pfam_name)
          p_empty <- p_empty %>%
            add_trace(
              data = NULL,
              x = as.vector(rbind(dom_i$pfam_start_g, dom_i$pfam_end_g, NA)),
              y = rep(c(-0.18, -0.18, NA), nrow(dom_i)),
              type = "scatter",
              mode = "lines",
              line = list(width = 10, color = if (!is.null(pfam_colors) && pfam_name %in% names(pfam_colors)) pfam_colors[[pfam_name]] else NULL),
              hoverinfo = "skip",
              showlegend = legend_this,
              legendgroup = pfam_name,
              name = pfam_name,
              inherit = FALSE
            )
        }
      }
      if (multi_tx_gene && length(extra_tx_ids) > 0) {
        p_empty <- build_multi_tx_gene_tracks(p_empty, plot_source, extra_tx_ids, extra_tx_y, pfam_colors = pfam_colors, legend_seen = legend_seen)
      }
      p_empty <- add_direction_markers(p_empty, meta)
      p_empty <- p_empty %>%
        layout(
          showlegend = TRUE,
          xaxis = list(title = "Genomic Position"),
          yaxis = if (multi_tx_gene && length(extra_tx_ids) > 0) {
            list(
              title = "SHaRe Classification",
              range = c(min(unname(extra_tx_y)) - 0.8, 4),
              tickvals = c(0, unname(extra_tx_y), 1, 2, 3),
              ticktext = c(main_tx, extra_tx_ids, genomic_classifications),
              fixedrange = TRUE
            )
          } else {
            list(title = "SHaRe Classification")
          }
        )
      return(event_register(p_empty, "plotly_click"))
    }

    p <- plot_ly(source = "gene_plot")
    if (nrow(exon) > 0) {
      p <- p %>%
        add_trace(
          data = NULL,
          x = as.vector(rbind(exon$exon_region_start_bp, exon$exon_region_end_bp, NA)),
          y = rep(c(0, 0, NA), nrow(exon)),
          type = "scatter",
          mode = "lines",
          line = list(color = "steelblue", width = 30),
          hoverinfo = "skip",
          showlegend = FALSE,
          name = "Exons",
          inherit = FALSE
        )
    }
    if (nrow(dom) > 0) {
      for (pfam_name in unique(dom$pfam_name)) {
        dom_i <- dom %>% filter(pfam_name == !!pfam_name)
        legend_this <- !(pfam_name %in% legend_seen)
        if (legend_this) legend_seen <- c(legend_seen, pfam_name)
        p <- p %>%
          add_trace(
            data = NULL,
            x = as.vector(rbind(dom_i$pfam_start_g, dom_i$pfam_end_g, NA)),
            y = rep(c(-0.18, -0.18, NA), nrow(dom_i)),
            type = "scatter",
            mode = "lines",
            line = list(width = 10, color = if (!is.null(pfam_colors) && pfam_name %in% names(pfam_colors)) pfam_colors[[pfam_name]] else NULL),
            hoverinfo = "skip",
            showlegend = legend_this,
            legendgroup = pfam_name,
            name = pfam_name,
            inherit = FALSE
          )
      }
    }
    if (multi_tx_gene && length(extra_tx_ids) > 0) {
      p <- build_multi_tx_gene_tracks(p, plot_source, extra_tx_ids, extra_tx_y, pfam_colors = pfam_colors, legend_seen = legend_seen)
    }

    for (class in intersect(genomic_classifications, levels(df$classification))) {
      df_class <- df[df$classification == class, ]

      p <- add_trace(
        p,
        data = df_class,
        x = ~POS,
        y = ~ypos,
        type = "scatter", mode = "markers",
        marker = list(
          symbol = "diamond-tall",
          color = clinvar_colors[[class]],
          line = list(color = "black", width = 1.5),
          size = ~size + 9
        ),
        text = ~paste0(
          "Variant: ", VariantID, "<br>",
          "Consequence: ", Consequence, "<br>",
          "SHaRe Classification: ", classification, "<br>",
          "ClinVar Classification: ", ClinVar_VarClass, "<br>"
        ),
        hoverinfo = "text",
        key = ~VariantID,
        name = class,
        showlegend = FALSE
      )
    }

    p <- add_direction_markers(p, meta)

    p <- p %>%
      layout(
        showlegend = TRUE,
        legend = list(
          itemsizing = "trace",
          itemclick = FALSE,
          itemdoubleclick = FALSE,
          font = list(color = "black")
        ),
        dragmode = "zoom",
        xaxis = list(
          title = "Genomic Position",
          showspikes = TRUE,
          spikemode = 'across',
          spikesnap = 'cursor',
          spikecolor = 'black',
          spikethickness = 1,
          tickformat = ",.0f",
          zeroline = FALSE,
          fixedrange = FALSE
        ),
        yaxis = list(
          range = if (multi_tx_gene && length(extra_tx_ids) > 0) {
            c(min(unname(extra_tx_y)) - 0.8, 4)
          } else {
            c(-0.8, 4)
          },
          title = "SHaRe Classification",
          tickvals = if (multi_tx_gene && length(extra_tx_ids) > 0) {
            c(0, unname(extra_tx_y), 1, 2, 3)
          } else {
            c(0, 1, 2, 3)
          },
          ticktext = if (multi_tx_gene && length(extra_tx_ids) > 0) {
            c(main_tx, extra_tx_ids, genomic_classifications)
          } else {
            c("Exons", genomic_classifications)
          },
          fixedrange = TRUE
        ),
        hovermode = "closest",
        clickmode = "event+select"
      ) %>%
      config(
        displayModeBar = TRUE,
        scrollZoom = TRUE,
        modeBarButtonsToRemove = list(
          "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d",
          "toggleSpikelines", 
          "zoom", "pan", "select", "zoomIn", "zoomOut", "autoScale"
        ),
        toImageButtonOptions = list(format = "png", width = 1400, height = 600, scale = 4),
        displaylogo = FALSE
      )

      event_register(p, "plotly_click")
      p
    }))
  })



  # Observe plot click -> scroll table and highlight
  observeEvent(event_data("plotly_click", source = "gene_plot"), {
    ed <- event_data("plotly_click", source = "gene_plot")
    if (!is.null(ed) && "key" %in% names(ed)) {
      sv <- as.character(ed$key[[1]])
      sel_var(sv)
      df <- gene_df()
      idx <- which(df$VariantID == sv)[1]
      if (length(idx) == 1 && !is.na(idx)) {
        updateReactable("variant_table", selected = idx)
        runjs(sprintf(
          "(function(){var t=document.getElementById('variant_table'); if(!t) return; 
          var rows=t.querySelectorAll('.rt-tr-group'); var i=%d-1; 
          if(rows[i]) rows[i].scrollIntoView({behavior:'smooth', block:'center'});})();",
          idx
        ))
      }
    }
  })

  
  
  output$variant_table <- renderReactable({
    df <- gene_df()
    if (nrow(df) == 0) return(NULL)

    status_cell <- function(value) {
      color <- switch(value,
                      "stop gained"  = "#dd2c00",
                      "splice donor" = "#dd2c00",
                      "splice acceptor" = "#dd2c00",
                      "frameshift" = "#dd2c00",

                      "missense" = "#ffa500",
                      "inframe insertion" = "#ffa500",
                      "inframe deletion" = "#ffa500",
                      "start lost" = "#ffa500",
                      "stop lost" = "#ffa500",

                      "splice region" = "#191919",
                      "splice donor region" = "#191919",
                      "splice donor 5th base" = "#191919",
                      "splice polypyrimidine tract" = "#191919",
                      "intron" = "#191919",
                      "5 prime UTR" = "#191919",
                      "3 prime UTR" = "#191919",

                      "synonymous" = "#2e7d32",
                      "#C0C0C0"
      )
      
      div(style = "display: flex; align-items: center;",
          div(style = paste0(
            "width: 12px; height: 12px; border-radius: 50%; margin-right: 8px;",
            "background-color:", color, ";",
            "border: 1px solid black;"
          )),
          span(value)
      )
    }

    df %>% 
      select("VariantID",
             "HGVS Coding Consequence" = cdot,
             "HGVS Protein Consequence" = pdot,
             "VEP_Consequence" = VEP_Consequence_display,
             "SHaRe AC" = AC,
             "SHaRe_VarClass" = VarClass ,
             "ClinVar_VarClass") %>% 
    reactable(
      selection = "single",
      highlight = TRUE,
      onClick = "select",
      pagination = FALSE,
      height = 600,
      searchable = TRUE,
      columns = list(
        VariantID= colDef(width = 250,
                          html = TRUE,
                          cell = function(value, index) {
                            payload <- jsonlite::toJSON(list(id = value, index = index), auto_unbox = TRUE)
                            as.character(tags$a(href = "#",onclick = sprintf("event.preventDefault(); event.stopPropagation(); Shiny.setInputValue('variant_navigate', %s, {priority: 'event'});", payload),value))
                            }),
        VEP_Consequence     = colDef(
          name = as.character(tags$span(
            "VEP Consequence ",
            tags$a(
              href = "https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html",
              target = "_blank",
              rel = "noopener noreferrer",
              title = "Open Ensembl consequence prediction guide",
              icon("circle-info", style = "font-size: 0.9em; vertical-align: middle; color: #2c3e50;")
            )
          )),
          html = TRUE,
          cell = status_cell,
          width = 320
        ),
        AC              = colDef(width = 100),
        SHaRe_VarClass  = colDef(cell = add_pill,name = "VarClass",maxWidth = 150,align = "center" ),
        ClinVar_VarClass= colDef(cell = add_pill,name = "ClinVar VarClass",maxWidth = 150,align = "center")
      ),
      details = colDef(
        cell = function(value, index) {
          as.character(
            tags$span(
              onclick = sprintf(
                "event.preventDefault(); event.stopPropagation(); Shiny.setInputValue('show_details', %d, {priority: 'event'}); return false;",
                index
              ),
              icon("info-circle", class = "text-primary", style = "font-size: 18px; cursor:pointer;")
            )
          )
        },
        html = TRUE,
        width = 80
      ),
      defaultColDef = colDef(vAlign = "center"),
      theme = reactableTheme(
        rowStyle = list(height = "50px"),
        rowSelectedStyle = list(backgroundColor = "#CEE1EE", boxShadow = "inset 5px 0 0 0 #ffa62d")
      )
    )
  })
  
  ###########
  # Modal UI
  ###########

  show_variant_modal <- function(row, index = NULL) {
    clinvar_modal_link <- if (!is.na(row$ClinVar[1]) && nzchar(row$ClinVar[1])) {
      HTML(sprintf("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/' target='_blank'>ClinVar (%s)</a>", row$ClinVar[1], row$ClinVar[1]))
    } else {
      HTML("ClinVar record not available")
    }
    clinvar_modal_class <- if (!is.na(row$ClinVar_VarClass[1]) && nzchar(as.character(row$ClinVar_VarClass[1]))) {
      row$ClinVar_VarClass[1]
    } else {
      "-"
    }

    showModal(modalDialog(
      title = div(
        style = "display: flex; justify-content: space-between; align-items: center;",
        h4("Variant Details", style = "margin: 0;"),
        create_spaced_buttons(
          plain_action_button(
            "Copy variant ID",
            onclick = modal_action_js("copy", row$VariantID[1])
          ),
          plain_action_button(
            "Open in Variant View",
            onclick = modal_action_js("open_variant", row$VariantID[1])
          ),
          plain_action_button(
            "Open in Gene View",
            onclick = modal_action_js("open_gene", row$VariantID[1])
          )
        )
      ),
      tags$head(
        tags$style(HTML(".modal-dialog {width: 1200px !important;max-width: 1200px;}.modal-body {min-height: 600px;}"))
      ),
      tags$script(HTML("Shiny.addCustomMessageHandler('copyText', function(txt) {navigator.clipboard.writeText(txt);});
    ")),
      size = "m",
      fluidRow(
        column(width = 6,  align = "center",offset = 1,
               h1(row$VariantID),
               reactableOutput("var_summary")
        )
      ),
      HTML("<br><br/>"),
      fluidRow(
        column(10, offset = 1,
               h2("SHaRe Classificaton"),
               reactableOutput("share_summary")
        )
      ),
      fluidRow(
        column(width = 10, offset = 1,
               h2("gnomAD version 4.1"),
               h4(HTML(paste0("<a href='", paste0("https://gnomad.broadinstitute.org/variant/",row$VariantID), "' target='_blank'>", row$VariantID, "</a>"))),
               reactableOutput("gnomad_summary"),
               h2("ClinVar"),
               h4(clinvar_modal_link),
               if (!is.na(row$ClinVar[1]) && nzchar(row$ClinVar[1])) {
                 tags$p(sprintf("ClinVar VarClass: %s", clinvar_modal_class), style = "margin-top: 0.25rem;")
               }
        )
      ),

      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    
    output$var_summary <- renderReactable({
      row %>%
        select(any_of(c("Gene", "Consequence","HGVSg", "HGVSc", "HGVSp"))) %>%
        distinct() %>%
        mutate(across(everything(), as.character)) %>%
        pivot_longer(cols = everything(),names_to = "column_name", values_to = "value") %>%
        reactable(borderless = TRUE,compact = T,
                  columns =  list(column_name = colDef(name = "",
                                                       align = "right",
                                                       maxWidth = 110,
                                                       style = function(value) {list(fontWeight = "bold")}),
                                  value = colDef(name = "",html = TRUE,align = "left")
                                  )
                  )
    })
    output$share_summary <- renderReactable({
      build_variant_share_summary_table(row)
    })

    output$gnomad_summary <- renderReactable({
      build_variant_gnomad_summary_table(row)
    })

    outputOptions(output, "share_summary", suspendWhenHidden = FALSE)
    outputOptions(output, "var_summary", suspendWhenHidden = FALSE)
    outputOptions(output, "gnomad_summary", suspendWhenHidden = FALSE)
  
  }
  
  
  
  observeEvent(input$variant_navigate, {
    data <- input$variant_navigate
    if (is.null(data)) return()
    if (is.list(data)) {
      variant_id <- data$id
    } else {
      variant_id <- data
    }
    if (is.null(variant_id) || !nzchar(variant_id)) return()
    select_variant(variant_id, navigate = TRUE)
  })

  
  
  #################
  # Variant Page UI
  #################
  
  
  output$variant_view_ui <- renderUI({
    row <- selected_variant()
    if (is.null(row) || nrow(row) == 0) {
      return(fluidRow(
        column(10, offset = 1,
               div(class = "alert alert-info", style = "margin-top: 2rem;",
                   p("Select a variant from the Gene View table or open a bookmarked variant URL to see full details.")
               )
        )
      ))
    }

    row <- row[1, , drop = FALSE]
    variant_id <- row$VariantID[1]
    chr <- row$CHR[1]
    pos <- row$POS[1]
    clinvar_id <- row$ClinVar[1]

    gnomad_link <- if (!is.na(variant_id) && nzchar(variant_id)) {
      HTML(sprintf("<a href='https://gnomad.broadinstitute.org/variant/%s' target='_blank'>%s</a>", variant_id, variant_id))
    } else {
      ""
    }
    clinvar_link <- if (!is.na(clinvar_id) && nzchar(clinvar_id)) {
      HTML(sprintf("<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/' target='_blank'>ClinVar (%s)</a>", clinvar_id, clinvar_id))
    } else {
      HTML("ClinVar record not available")
    }
    ensembl_link <- if (!is.na(chr) && !is.na(pos)) {
      region <- sprintf("%s:%s-%s", chr, pos, pos)
      HTML(sprintf("<a href='https://www.ensembl.org/Homo_sapiens/Location/View?r=%s' target='_blank'>Ensembl browser (region)</a>", region))
    } else {
      HTML("Ensembl browser (region)")
    }
    varsome_link <- if (!is.na(variant_id) && nzchar(variant_id)) {
      HTML(sprintf("<a href='https://varsome.com/variant/hg38/%s?annotation-mode=germline' target='_blank'>VarSome (%s)</a>", variant_id, variant_id))
    } else {
      HTML("VarSome")
    }

    #shinycssloaders::withSpinner(
    tagList(
      tags$script(HTML("if (typeof window.shareCopyHandler === 'undefined') {window.shareCopyHandler = true; Shiny.addCustomMessageHandler('copyText', function(txt) {navigator.clipboard.writeText(txt);});}")),

      
      
      fluidRow(
        column(width = 10, offset = 1,
               div(
                 style = "text-align:center;",
                 h1(variant_id, style = "margin: 0 0 0.75rem 0; width: 100%"),
                 div(
                   style = "display: flex; justify-content: center; gap: 0.5rem; margin-bottom: 0.75rem;",
                   create_spaced_buttons(
                     plain_action_button("Copy variant ID", "Shiny.setInputValue('variant_view_copy_btn', {nonce: Date.now()}, {priority: 'event'});"),
                     plain_action_button("Open in Gene View", "Shiny.setInputValue('variant_view_gene_btn', {nonce: Date.now()}, {priority: 'event'});")
                   )
                 ),
                 reactableOutput("variant_view_summary")
                 
               )
        )
      ),
      br(), br(),
      
      fluidRow(
        column(width = 5,offset = 1,
               h2("SHaRe Classification"),
               reactableOutput("variant_view_share_table")
        )
      ),
      br(), br(),
      
      
      fluidRow(
        column(width = 5, offset = 1,
               
               h2("GnomAD version 4.1"),
               h4(HTML(paste0("<a href='", paste0("https://gnomad.broadinstitute.org/variant/",row$VariantID), "' target='_blank'>", row$VariantID, "</a>"))),
               reactableOutput("variant_view_gnomad_joint")
        )
      ),
      br(), br(),
      
      
      fluidRow(
        column(
          width = 10,offset = 1,
          div(
            style = "max-width: 100%;",
            tags$div(
              style = "display: flex; align-items: center; gap: 0.35rem;",
              h2("GnomAD Genetic Ancestry Group Frequencies", style = "margin: 0;"),
              actionLink(
                "gnomad_info_btn",
                label = NULL,
                icon = icon("circle-info"),
                style = "color: #2c3e50; text-decoration: none; font-size: 0.9rem;"
              )
            ),
            reactableOutput("variant_view_gnomad_pop")
          )
        )
      ),
      br(), br(),
      
      
      fluidRow(
        column(
          width = 5,offset = 1,
          div(
            style = "max-width: 100%; margin-top: 1rem;",
            h2("Case - Control Odds Ratio (SHaRe vs. gnomAD)"),
            uiOutput("variant_view_or")
          )
        )
      ),
      br(), br(),

      fluidRow(
        column(width = 5,offset = 1,
               h2("SHaRe Classification Notes"),
               reactableOutput("variant_view_classification_notes")
        )
      ),
      br(), br(),
      
      fluidRow(
        column(width = 5,offset = 1,
               h2("ClinVar"),
               h4(clinvar_link),
               tags$ul(
                 tags$li(sprintf("Clinical Significance: %s", row$Clinvar_ClinicalSignificance)),
                 tags$li(sprintf("Last Evaluated: %s", row$Clinvar_LastEvaluated)),
                 tags$li(sprintf("Review Status: %s", row$Clinvar_ReviewStatus)),
                 tags$li(sprintf("Number of Submitters: %s", row$Clinvar_NumberSubmitters)),
                 tags$li(sprintf("Phenotype List: %s", row$Clinvar_PhenotypeList)),
                 tags$li(sprintf("ClinVar VarClass: %s", row$ClinVar_VarClass))
                 )
               )
        ),
      br(), br(),
      
      fluidRow(
        column(width = 5,offset = 1,
               h2("VarSome"),
               h4(varsome_link)
        )
      ),
      br(), br(),
      
      fluidRow(
        column(width = 5, offset = 1,
               h2("Ensembl Variant Effect Predictor"),
               h4(ensembl_link),
               tags$ul(
                 tags$li(sprintf("Consequence: %s", row$Consequence)),
                 tags$li(sprintf("Impact: %s", row$Impact)),
                 tags$li(sprintf("PTV: %s", row$Flag_PTV)),
                 tags$li(sprintf("PAV: %s", row$Flag_PAV))
               ),
               h2("In Silico Predictors"),
               tags$ul(
                 tags$li("LOFTEE",
                         tags$ul(
                           tags$li(sprintf("LoF: %s", row$LoF)),
                           tags$li(sprintf("LoF filters: %s", row$LoF_filter)),
                           tags$li(sprintf("LoF flags: %s", row$LoF_flags))
                         )
                 ),
                 tags$li("AlphaMissense",
                         tags$ul(
                           tags$li(sprintf("AM class: %s", row$am_class)),
                           tags$li(sprintf("AM pathogenicity: %s", row$am_pathogenicity))
                         )
                 ),
                 tags$li("REVEL",
                         tags$ul(
                           tags$li(sprintf("REVEL: %s", row$REVEL))
                         )
                 ),
                 tags$li("SpliceAI",
                         tags$ul(
                           tags$li(sprintf("SpliceAI gene: %s", row$SpliceAI_GENE)),
                           tags$li(sprintf("SpliceAI phred: %s", row$SpliceAI_pred)),
                           tags$li(sprintf("SpliceAI max: %s", row$SpliceAI_max))
                         )
                 ),
                 tags$li("CADD",
                         tags$ul(
                           tags$li(sprintf("CADD phred: %s", row$CADD_PHRED)),
                           tags$li(sprintf("CADD raw: %s", row$CADD_RAW))
                         )
                 )
               )
               )
      ),
      br(), br(),
      br(), br()
    )
    #)
  })

  build_variant_share_summary_table <- function(row) {
    row %>%
      select(
        "Allele Count" = AC,
        "Number of Carriers" = N_carriers,
        "Number of SHaRe Sites" = N_sites,
        VarClass,
        "VarClass Date" = VarClass_date
      ) %>%
      rename_with(~ str_replace_all(.x, "_", " ")) %>%
      reactable(
        columns = list(
          VarClass = colDef(name = "VarClass", cell = add_pill),
          VariantID = colDef(width = 250)
        ),
        defaultColDef = colDef(align = "left")
      )
  }

  fmt_gnomad_frequency <- function(x) {
    vapply(x, function(value) {
      if (length(value) == 0 || is.na(value) || !nzchar(as.character(value))) return("-")
      value <- suppressWarnings(as.numeric(value))
      if (is.na(value)) return("-")
      if (value == 0) return("0.0000")

      s <- format(value, scientific = FALSE, trim = TRUE, digits = 22)
      if (!str_detect(s, fixed("."))) {
        return(s)
      }

      parts <- str_split_fixed(s, fixed("."), 2)
      decimals <- str_replace(parts[, 2], "0+$", "")
      leading_zeros <- str_extract(decimals, "^0*")
      significant_decimals <- str_sub(decimals, nchar(leading_zeros) + 1L)

      if (!nzchar(significant_decimals)) {
        return(paste0(parts[, 1], ".", str_dup("0", 4)))
      }

      paste0(parts[, 1], ".", leading_zeros, str_sub(significant_decimals, 1L, 4L))
    }, character(1), USE.NAMES = FALSE)
  }

  build_variant_gnomad_summary_table <- function(row) {
    fmt_dash_num <- function(x) {
      if (is.na(x) || length(x) == 0 || !nzchar(as.character(x))) "-" else format(x, scientific = FALSE)
    }
    fmt_dash_chr <- function(x) {
      x <- as.character(x)
      if (is.na(x) || !nzchar(x)) "-" else x
    }
    gnomad_df <- row %>%
      select(
        "Allele Count" = gnomAD4_AC_joint,
        "Allele Number" = gnomAD4_AN_joint,
        "Number of Homozygotes" = gnomAD4_nhomalt_joint,
        "Allele Frequency" = gnomAD4_AF_joint,
        "Filtering AF (95% confidence) GroupMax" = gnomAD4_fafmax_faf95_max_joint,
        "Filtering AF (95% confidence) GroupMax Ancestry" = gnomAD4_fafmax_faf95_max_gen_anc_joint
      ) %>%
      distinct() %>%
      mutate(across(everything(), ~ {
        if (is.numeric(.x)) {
          fmt_dash_num(.x)
        } else {
          fmt_dash_chr(.x)
        }
      })) %>%
      mutate(
        `Allele Frequency` = fmt_gnomad_frequency(`Allele Frequency`),
        `Filtering AF (95% confidence) GroupMax` = fmt_gnomad_frequency(`Filtering AF (95% confidence) GroupMax`)
      ) %>%
      mutate(`Filtering AF (95% confidence) GroupMax Ancestry` = if_else(
        `Filtering AF (95% confidence) GroupMax Ancestry` == "-",
        "-",
        toupper(`Filtering AF (95% confidence) GroupMax Ancestry`)
      ))

    if (ncol(gnomad_df) == 0) return(NULL)

    reactable(gnomad_df, borderless = TRUE, compact = TRUE)
  }

  # variant view table summary
  output$variant_view_summary <- renderReactable({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    row %>%
      select(any_of(c("Gene", "Consequence", "Impact","SPDI","HGVSg", "HGVSc", "HGVSp", "EXON", "INTRON", "Amino_Acids","Codons"))) %>%
      distinct() %>%
      mutate(across(everything(), as.character)) %>%
      rename_with(~ gsub("_", " ", .x)) %>% 
      pivot_longer(cols = everything(), names_to = "column_name", values_to = "value") %>%
      reactable(borderless = TRUE, compact = TRUE,
                columns = list(
                  column_name = colDef(name = "", align = "right", #maxWidth = 140,
                                       style = function(value) {list(fontWeight = "bold")}),
                  value = colDef(name = "", html = TRUE, align = "left")
                ))
  })

  # share classification table
  output$variant_view_share_table <- renderReactable({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    build_variant_share_summary_table(row)
  })

  # share notes table
  output$variant_view_classification_notes <- renderReactable({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    row %>%
      select("VarClass Status" = VarClass_status,
             "ACMG Rules" = ACMG,
             "ACMG Score" = Score,
             "Classification with VUS subclass" = VUS_Class,
             "MYBPC3 SpliceStatus" = MYBPC3_SpliceStatus) %>% 
      mutate(across(everything(), ~ {
        x <- as.character(.x)
        x[is.na(x) | !nzchar(x)] <- "-"
        x
      })) %>%
      pivot_longer(cols = everything(), names_to = "column_name", values_to = "value") %>%
      reactable(
        borderless = TRUE,
        compact = TRUE,
        pagination = FALSE,
        columns = list(
          column_name = colDef(
            name = "",
            align = "right",
            maxWidth = 160,
            style = function(value) {list(fontWeight = "bold")}
          ),
          value = colDef(name = "", align = "left")
        )
      )
  })

  build_variant_view_gnomad_pop <- function(row) {
    pop_map <- tibble::tibble(
      pop = c("afr", "ami", "amr", "asj", "fin", "mid", "nfe", "sas", "remaining", "xx", "xy"),
      label = c(
        "African/African American",
        "Amish",
        "Admixed American",
        "Ashkenazi Jewish",
        "Finnish",
        "Middle Eastern",
        "Non-Finnish European",
        "South Asian",
        "Remaining",
        "XX",
        "XY"
      )
    )
    
    
    clean_decimals <- function(x) {
      fmt_gnomad_frequency(x)
    }

    fmt_count <- function(x) {
      vapply(x, function(value) {
        if (is.na(value) || length(value) == 0 || !nzchar(as.character(value))) {
          "-"
        } else {
          format(value, scientific = FALSE, trim = TRUE)
        }
      }, character(1), USE.NAMES = FALSE)
    }
    
    
    get_pop_val <- function(prefix, pop) {
      candidates <- unique(c(
        paste0(prefix, pop),
        paste0(prefix, toupper(pop))
      ))
      nm <- candidates[candidates %in% names(row)][1]
      if (is.na(nm) || !nzchar(nm)) NA else row[[nm]][1]
    }

    fmt_total_count <- function(x) {
      if (is.na(x) || length(x) == 0 || !nzchar(as.character(x))) {
        "-"
      } else {
        fmt_count(x)
      }
    }

    fmt_total_freq <- function(x) {
      clean_decimals(x)
    }

    pop_rows <- purrr::pmap_dfr(pop_map, function(pop, label) {
      tibble(
        Population = label,
        "Allele Count" = get_pop_val("gnomAD4_AC_joint_", pop),
        "Allele Number" = get_pop_val("gnomAD4_AN_joint_", pop),
        "Number of Homozygotes" = get_pop_val("gnomAD4_nhomalt_joint_", pop),
        "Allele Frequency" = get_pop_val("gnomAD4_AF_joint_", pop),
        "Filtering AF (95% confidence)" = get_pop_val("gnomAD4_faf95_joint_", pop)
      )
    }) %>%
      mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>% 
      mutate(
        `Allele Count` = fmt_count(`Allele Count`),
        `Allele Number` = fmt_count(`Allele Number`),
        `Number of Homozygotes` = fmt_count(`Number of Homozygotes`),
        `Allele Frequency` = clean_decimals(`Allele Frequency`),
        `Filtering AF (95% confidence)` = clean_decimals(`Filtering AF (95% confidence)`)
      ) %>%
      mutate(across(-Population, as.character)) %>% 
      mutate(.sex_block = Population %in% c("XX", "XY")) %>%
      arrange(.sex_block, desc(`Filtering AF (95% confidence)`)) %>%
      select(-.sex_block)

    total_row <- tibble(
      Population = "Total",
      "Allele Count" = fmt_total_count(row$gnomAD4_AC_joint[1]),
      "Allele Number" = fmt_total_count(row$gnomAD4_AN_joint[1]),
      "Number of Homozygotes" = "-",
      "Allele Frequency" = fmt_total_freq(row$gnomAD4_AF_joint[1]),
      "Filtering AF (95% confidence)" = fmt_total_freq(row$gnomAD4_fafmax_faf95_max_joint[1])
    )

    bind_rows(pop_rows, total_row)
  }

  
  build_variant_view_gnomad_joint <- function(row) {
    fmt_joint_num <- function(x) {
      if (is.na(x) || length(x) == 0 || !nzchar(as.character(x))) "-" else format(x, scientific = FALSE)
    }
    fmt_joint_chr <- function(x) {
      x <- as.character(x)
      if (is.na(x) || !nzchar(x)) "-" else x
    }

    tibble(
      Metric = c(
        "Allele Count",
        "Allele Number",
        "Allele Frequency",
        "Filtering AF (95% confidence) GroupMax",
        "Filtering AF (95% confidence) GroupMax Ancestry"
      ),
      Value = c(
        fmt_joint_num(row$gnomAD4_AC_joint[1]),
        fmt_joint_num(row$gnomAD4_AN_joint[1]),
        fmt_gnomad_frequency(row$gnomAD4_AF_joint[1]),
        fmt_gnomad_frequency(row$gnomAD4_fafmax_faf95_max_joint[1]),
        {
          anc <- fmt_joint_chr(row$gnomAD4_fafmax_faf95_max_gen_anc_joint[1])
          if (identical(anc, "-")) "-" else toupper(anc)
        }
      )
    ) %>%
      mutate(across(everything(), as.character))
  }

  build_variant_view_or <- function(row) {
    case_alt <- suppressWarnings(as.numeric(dplyr::coalesce(row$SHaRe_AC[1], row$AC[1])))
    if (is.na(case_alt)) case_alt <- 0
    case_an <- 13935 * 2
    case_ref <- max(case_an - case_alt, 0)
    case_af <- if (case_an > 0) case_alt / case_an else NA_real_

    ctrl_alt <- suppressWarnings(as.numeric(row$gnomAD4_AC_joint[1]))
    if (is.na(ctrl_alt)) ctrl_alt <- 0
    ctrl_an <- suppressWarnings(as.numeric(row$gnomAD4_AN_joint[1]))
    if (is.na(ctrl_an)) ctrl_an <- 0
    ctrl_ref <- max(ctrl_an - ctrl_alt, 0)
    ctrl_af <- if (ctrl_an > 0) ctrl_alt / ctrl_an else NA_real_

    tab <- matrix(
      c(case_alt, case_ref, ctrl_alt, ctrl_ref),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("Case", "Control"),
        c("Alt", "Ref")
      )
    )

    ft <- fisher.test(tab)
    or_est <- suppressWarnings(as.numeric(ft$estimate[[1]]))
    ci_low <- suppressWarnings(as.numeric(ft$conf.int[1]))
    ci_high <- suppressWarnings(as.numeric(ft$conf.int[2]))
    x_max <- if (is.finite(or_est) && or_est > 0) {
      max(or_est, 1)
    } else if (is.finite(ci_high) && ci_high > 0) {
      max(ci_high, 1)
    } else {
      1
    }
    ci_high_plot <- if (is.finite(ci_high)) min(ci_high, x_max) else x_max
    ci_low_plot <- if (is.finite(ci_low)) max(ci_low, 0) else 0

    list(
      table = tibble(
        Group = c("Case", "Control"),
        Alt = c(case_alt, ctrl_alt),
        Ref = c(case_ref, ctrl_ref),
        AF = fmt_gnomad_frequency(c(case_af, ctrl_af))
      ),
      summary = sprintf(
        "Fisher exact OR = %.2f (95%% CI %.2f to %.2f)",
        or_est,
        ci_low,
        ci_high
      ),
      or = or_est,
      ci_low = ci_low,
      ci_high = ci_high,
      x_max = x_max,
      ci_low_plot = ci_low_plot,
      ci_high_plot = ci_high_plot
    )
  }

  # gnomad summary freq table
  output$variant_view_gnomad_joint <- renderReactable({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    build_variant_gnomad_summary_table(row)
  })

  output$variant_view_or <- renderUI({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    or_data <- tryCatch(build_variant_view_or(row), error = function(e) NULL)
    if (is.null(or_data)) {
      return(tags$div(
        style = "margin: 0.5rem 0 0.75rem 0;",
      tags$em("Odds ratio could not be calculated for this variant.")
      ))
    }
    tagList(
      div(
        style = "text-align: left;",
        # tags$p(or_data$summary, style = "margin-top: 0.15rem; margin-bottom: 0.25rem; font-weight: bold;"),
        h4(or_data$summary),
        reactable(
          or_data$table,
          borderless = TRUE,
          compact = TRUE,
          pagination = FALSE,
          defaultColDef = colDef(align = "center"),
          columns = list(
            Group = colDef(name = "Group", align = "left", style = function(value) list(fontWeight = "bold")),
            Alt   = colDef(name = "Allele Count", align = "right"),
            Ref   = colDef(name = "Allele Number", align = "right"),
            AF    = colDef(name = "Allele Frequency",minWidth = 200, align = "left")
          )
        )# ,
        # tags$div(style = "height: 0.5rem;"),
        # plotOutput("variant_view_or_plot", height = "72px")
      )
    )
  })

  output$variant_view_or_plot <- renderPlot({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    or_data <- tryCatch(build_variant_view_or(row), error = function(e) NULL)
    req(!is.null(or_data))

    plot_df <- tibble(
      y = 1,
      or = or_data$or,
      ci_low = or_data$ci_low_plot,
      ci_high = or_data$ci_high_plot
    )

    ggplot(plot_df) +
      geom_segment(
        aes(x = ci_low, xend = ci_high, y = y, yend = y),
        colour = "black",
        linewidth = 0.35
      ) +
      geom_segment(
        aes(x = ci_low, xend = ci_low, y = y - 0.06, yend = y + 0.06),
        colour = "black",
        linewidth = 0.35
      ) +
      geom_segment(
        aes(x = ci_high, xend = ci_high, y = y - 0.06, yend = y + 0.06),
        colour = "black",
        linewidth = 0.35
      ) +
      #geom_vline(xintercept = 1, colour = "#d0d0d0", linewidth = 0.25) +
      geom_point(aes(x = or, y = y), shape = 15, size = 4, colour = "black") +
      scale_x_continuous(limits = c(0, or_data$x_max), expand = expansion(mult = c(0.08, 0.12))) +
      scale_y_continuous(limits = c(0.8, 1.2), breaks = NULL) +
      coord_cartesian(clip = "off") +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.x = element_text(size = 10, colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.ticks.x = element_line(colour = "#999999", linewidth = 0.2),
        plot.margin = margin(1, 1, 1, 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
        # panel.grid.major.x = element_line(colour = "#ededed", linewidth = 0.2)
      )
  })

  # gnomad population freq table
  output$variant_view_gnomad_pop <- renderReactable({
    row <- selected_variant()
    req(!is.null(row), nrow(row) > 0)
    reactable(
      build_variant_view_gnomad_pop(row),
      borderless = FALSE,
      compact = TRUE,
      pagination = FALSE,
      defaultColDef = colDef(align = "center"),
      columns = list(
        Population = colDef(minWidth = 200, align = "left"),
        "Allele Count" = colDef(minWidth = 80, align = "right"),
        "Allele Number" = colDef(minWidth = 100, align = "right"),
        "Number of Homozygotes" = colDef(minWidth = 120, align = "right"),
        "Allele Frequency" = colDef(minWidth = 130, align = "left"),
        "Filtering AF (95% confidence)" = colDef(minWidth = 150, align = "left")
      ),
      rowStyle = function(index) {
        rows <- build_variant_view_gnomad_pop(row)
        if (!length(index)) return(list())
        if (rows$Population[index[1]] == "XX") {
          list(borderTop = "3px solid #000")
        } else if (rows$Population[index[1]] == "Total") {
          list(borderTop = "3px solid #000", fontWeight = "bold")
        } else {
          list()
        }
      }
    )
  })

  
  
  
  ######## Listen for button clicks ######## 
  observeEvent(input$variant_view_copy_btn, {
    row <- selected_variant()
    if (is.null(row) || nrow(row) == 0) return()
    session$sendCustomMessage("copyText", row$VariantID[1])
  })

  observeEvent(input$variant_view_gene_btn, {
    row <- selected_variant()
    if (is.null(row) || nrow(row) == 0) return()
    open_gene_view(row$VariantID[1])
  })

  observeEvent(input$show_details, {
    # Extract the index from the button ID
    idx <- as.numeric(gsub("details_", "", input$show_details))
    
    df <- gene_df()
    if (idx > 0 && idx <= nrow(df)) {
      row <- df[idx, ]
      show_variant_modal(row)
    }
  })

  # Modal action dispatch (copy / open gene or variant)
  observeEvent(input$modal_action, {
    data <- input$modal_action
    req(is.list(data), !is.null(data$action), !is.null(data$variant))
    action <- data$action
    variant_id <- data$variant
    req(nzchar(variant_id))

    if (identical(action, "copy")) {
      session$sendCustomMessage("copyText", variant_id)
      return()
    }

    removeModal()
    if (identical(action, "open_variant")) {
      select_variant(variant_id, navigate = TRUE)
    } else if (identical(action, "open_gene")) {
      open_gene_view(variant_id)
    }
  }, ignoreNULL = TRUE)


  # Observe table click -> update plot and scroll
  observeEvent(input$variant_table_selected, {
    df <- gene_df()
    idx <- input$variant_table_selected
    if (!is.null(idx) && length(idx) >= 1) {
      i <- idx[[1]]
      sv <- df$VariantID[i]
      if (!identical(isolate(sel_var()), sv)) {
        sel_var(sv)
      }
      if (identical(input$navbar, "Gene View")) {
        runjs("var el=document.getElementById('gene_plot'); if(el) el.scrollIntoView({behavior:'smooth', block:'start'});")
      }
      updateReactable("variant_table", selected = i)
    }
  })
  observeEvent(sel_var(), {
    df <- tryCatch(gene_df(), error = function(e) NULL)
    sv <- sel_var()
    proxy <- plotlyProxy("gene_plot", session)

    valid_variant <- !is.null(sv) && !is.null(df) && nrow(df) > 0 && sv %in% df$VariantID

    if (!is.null(proxy)) {
      if (valid_variant) {
        pos <- df$POS[df$VariantID == sv][1]
        plotlyProxyInvoke(
          proxy,
          "relayout",
          list(shapes = list(list(
            type = "line",
            x0 = pos,
            x1 = pos,
            y0 = -1,
            y1 = 4,
            line = list(dash = "dash", width = 1, color = "black")
          )))
        )
      } else {
        plotlyProxyInvoke(proxy, "relayout", list(shapes = list()))
      }
    }

    if (valid_variant) {
      idx <- which(df$VariantID == sv)[1]
      current_sel <- isolate(input$variant_table_selected)
      current_idx <- if (!is.null(current_sel) && length(current_sel) >= 1) current_sel[[1]] else NULL
      if (is.null(current_idx) || current_idx != idx) {
        updateReactable("variant_table", selected = idx)
      }
      runjs(sprintf(
        "(function(){var t=document.getElementById('variant_table'); if(!t) return; var rows=t.querySelectorAll('.rt-tr-group'); var i=%d-1; if(rows[i]) rows[i].scrollIntoView({behavior:'smooth', block:'center'});})();",
        idx
      ))
      update_url_path(mode = "push")
    } else {
      updateReactable("variant_table", selected = NULL)
      update_url_path()
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$navbar, {
    if (identical(input$navbar, "Variant View")) {
      runjs("window.scrollTo({top: 0, behavior: 'smooth'});")
    }
  }, ignoreNULL = TRUE)


  
  
  
  
  
  #################
  # Gene Page UI
  #################
  
  output$gene_view_ui <- renderUI({
    on_gene_tab <- identical(input$navbar, "Gene View")
    body_style <- ""
    if (!on_gene_tab) return(NULL)
    if (isTRUE(gene_loading()) && !has_active_gene()) {
      return(fluidRow(
        column(
          10,
          offset = 1,
          div(
            style = "margin-top: 3rem; display: flex; justify-content: center;",
            shinycssloaders::withSpinner(div(style = "min-height: 120px;"), type = 4, color = "#2c3e50")
          )
        )
      ))
    }
    if (!has_active_gene()) {
      return(fluidRow(
        column(
          10,
          offset = 1,
          div(
            class = "alert alert-info",
            style = "margin-top: 2rem;",
            p("Select an HCM gene from the Home tab to load the overview.")
          )
        )
      ))
    }
    
    
    tagList(
      div(
        style = body_style,
        fluidRow(
          column(width = 10, offset = 1,
                 div(
                   style = "text-align:center;",
                   h1(textOutput("gene_title", inline = TRUE), style = "margin: 0 0 0.75rem 0; width: 100%;"),
                   reactableOutput("gene_summary")
                 )
          ),
        ),
        br(), br(),
        
        fluidRow(
          column(width = 5,offset = 1,
                 div(
                   style = "padding-left: 1rem; margin-top: 1rem;",
                   h2("SHaRe Variant Classification Summary", style = "margin-bottom: 0.5rem;"),
                   reactableOutput("gene_clinvar_summary")
                 )
          )
        ),
        #br(), br(),
        
        fluidRow(
          column(width = 10,offset = 1,
                 div(
                   style = "margin-top: 1rem;",
                   plotlyOutput("gene_plot", height = 300)
                 )
          )
        ),
        #br(), br(),
        
        fluidRow(
          column(width = 10,offset = 1,
                 prettyCheckboxGroup(
                   inputId = "select_class",
                   label = "Filter by Classification:",
                   choiceNames = c("P/LP", "VUS", "B/LB"),
                   choiceValues = c("P/LP", "VUS", "B/LB"),
                   selected = c("P/LP", "VUS", "B/LB"),
                   status = "default",
                   outline = TRUE,
                   shape = "round",
                   icon = icon("check"),
                   inline = TRUE,
                   bigger = TRUE,
                   thick = TRUE,
                   animation = "smooth"
                 )
          )
        ),
        # br(), br(),
        
        fluidRow(
          column(width = 10,offset = 1,
                 reactableOutput("variant_table")
          )
        ),
        br(), br(),
        br(), br()
      )
    )
  })

  output$gene_title <- renderText({
    req(has_active_gene())
    active_gene()
  })

  build_gene_summary_table <- function(gene_symbol) {
      df <- gene_info_all %>%
        filter(gene_name == gene_symbol) %>%
      select(-contains(c("5_","3_","exon","pfam","mrna"))) %>% 
      distinct(transcript_stable_id_version, .keep_all = TRUE) %>%
      mutate(
        transcript_rank = case_when(
          MANE_status == "MANE Select" ~ 1L,
          str_detect(MANE_status, "MANE Plus") ~ 2L,
          TRUE ~ 3L
        )
      ) %>%
      arrange(transcript_rank, transcript_stable_id_version)

    if (nrow(df) == 0) return(NULL)

    main <- df %>% slice(1)
    extra <- if (nrow(df) > 1) df %>% slice(2) else NULL

    format_transcript_value <- function(row) {
      pieces <- c(
        row$transcript_stable_id_version,
        row$refseq_match_transcript_mane_select
      )
      pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
      paste(pieces, collapse = " / ")
    }

    transcript_rows <- tibble(
      field = c(
        if (str_detect(as.character(main$MANE_status), "MANE Select")) "MANE Select" else "Transcript",
        if (!is.null(extra) && nrow(extra) == 1) {
          "MANE Plus Clinical"
        } else {
          NULL
        },
        if (!is.na(main$Ensembl_prot) && nzchar(as.character(main$Ensembl_prot))) {
          "Ensembl protein ID"
        } else {
          NULL
        }
      ),
      value = c(
        format_transcript_value(main),
        if (!is.null(extra) && nrow(extra) == 1) format_transcript_value(extra) else NULL,
        if (!is.na(main$Ensembl_prot) && nzchar(as.character(main$Ensembl_prot))) as.character(main$Ensembl_prot) else NULL
      )
    )

    gene_rows <- bind_rows(
      tibble(
        field = c(
          "Gene description",
          as.character(tags$span(
            "Disease ",
            tags$a(
              href = "https://imperialcardiogenetics.github.io/G2P-Cardiac-Panel/#hypertrophic-cardiomyopathy-hcm",
              target = "_blank",
              rel = "noopener noreferrer",
              title = "Open G2P HCM panel",
              icon("circle-info", style = "font-size: 0.9em; color: #2c3e50;")
            )
          )),
          "Evidence strength",
          "Inheritance",
          "Allelic requirement",
          "Genetic mechanism",
          "Ensembl gene ID"
        ),
        value = c(
          main$name,
          main$Disease,
          main$Evidence_strength,
          main$Inheritance,
          main$Allelic_requirement,
          main$Genetic_mechanism,
          paste0(
            "<a href='https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
            main$gene_stable_id_version,
            "' target='_blank'>",
            main$gene_stable_id_version,
            "</a>"
          )
        )
      ),
      transcript_rows,
      tibble(
        field = c(
          "UniProt ID",
          "Strand",
          "Genomic region"
        ),
        value = c(
          paste0(
            "<a href='https://www.uniprot.org/uniprotkb/",
            main$uniprot_id,
            "/entry' target='_blank'>",
            main$uniprot_id,
            "</a>"
          ),
          case_when(
            main$strand == -1 ~ "Reverse",
            main$strand == 1 ~ "Forward",
            TRUE ~ as.character(main$strand)
          ),
          paste0(main$hg38_chr, ":", main$gene_start_bp, "-", main$gene_end_bp)
        )
      )
    ) %>% mutate(across(everything(), as.character))
  }

  output$gene_summary <- renderReactable({
    req(has_active_gene())
    gene_summary_df <- build_gene_summary_table(active_gene())
    req(!is.null(gene_summary_df))

    reactable(
      gene_summary_df,
      borderless = TRUE,
      defaultPageSize = 25,
      compact = TRUE,
      columns = list(
        field = colDef(
          name = "",
          align = "right",
          html = TRUE,
          style = function(value) list(fontWeight = "bold")
        ),
        value = colDef(
          name = "",
          html = TRUE,
          align = "left"
        )
      )
    )
  })

  output$gene_clinvar_summary <- renderReactable({
    req(has_active_gene())
    df <- gene_clinvar_counts[[active_gene()]]
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df %>%
      mutate(
        VarClass = factor(VarClass, levels = c("P/LP", "VUS", "B/LB"))
      ) %>%
      arrange(VarClass) %>%
      select(VarClass, n) %>%
      pivot_wider(names_from = VarClass, values_from = n, values_fill = 0) %>%
      mutate(across(everything(), as.character)) %>% 
      bind_rows(setNames(as.list(colnames(.)), colnames(.)), .) %>% 
      reactable(
        borderless = TRUE,
        compact = TRUE,
        pagination = FALSE,
        defaultColDef = reactable::colDef(
          align = "center",
          html = TRUE,
          header = "",
              cell = function(value, index, name) {
                if (index == 1) {
              add_pill(value)
                } else {
              value
                }
              }
        )
      )

  })

  outputOptions(output, "gene_view_ui", suspendWhenHidden = FALSE)
  outputOptions(output, "gene_title", suspendWhenHidden = FALSE)
  outputOptions(output, "gene_summary", suspendWhenHidden = FALSE)
  outputOptions(output, "gene_clinvar_summary", suspendWhenHidden = FALSE)
  outputOptions(output, "gene_plot", suspendWhenHidden = FALSE)
  outputOptions(output, "variant_table", suspendWhenHidden = FALSE)

  
}


shinyApp(ui, server)

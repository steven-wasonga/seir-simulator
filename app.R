# ==============================================================================
# SEIR Epidemic Simulator
# Interactive outbreak modelling for public health & epidemiology
# ==============================================================================
# Install deps (run once):
# install.packages(c("shiny","bslib","plotly","deSolve","dplyr","tidyr","scales","DT"))
# ==============================================================================

library(shiny)
library(bslib)
library(plotly)
library(deSolve)
library(dplyr)
library(tidyr)
library(scales)
library(DT)

# ==============================================================================
# DISEASE PRESETS
# ==============================================================================

PRESETS <- list(
  "COVID-19 (Ancestral)"  = list(R0 = 2.8,  incubation = 5,  infectious = 5,  ifr = 1.2,  hosp_rate = 8),
  "COVID-19 (Omicron)"    = list(R0 = 10.0, incubation = 3,  infectious = 4,  ifr = 0.3,  hosp_rate = 3),
  "Influenza (Seasonal)"  = list(R0 = 1.3,  incubation = 2,  infectious = 4,  ifr = 0.1,  hosp_rate = 2),
  "Influenza (1918 H1N1)" = list(R0 = 2.2,  incubation = 2,  infectious = 4,  ifr = 2.5,  hosp_rate = 12),
  "Measles"               = list(R0 = 15.0, incubation = 12, infectious = 8,  ifr = 0.2,  hosp_rate = 5),
  "Ebola (West Africa)"   = list(R0 = 1.8,  incubation = 9,  infectious = 7,  ifr = 50.0, hosp_rate = 80),
  "SARS-CoV-1"            = list(R0 = 3.0,  incubation = 5,  infectious = 7,  ifr = 9.6,  hosp_rate = 25),
  "Mpox"                  = list(R0 = 2.1,  incubation = 10, infectious = 14, ifr = 1.0,  hosp_rate = 10),
  "Custom"                = list(R0 = 2.5,  incubation = 5,  infectious = 7,  ifr = 1.0,  hosp_rate = 5)
)

# ==============================================================================
# ODE MODEL
# ==============================================================================

seir_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    N <- S + E + I + R + V
    
    beta_t     <- if (t >= t_npi) beta * (1 - npi_eff) else beta
    vax_rate_t <- if (t >= t_vax && t < (t_vax + vax_duration)) vax_daily_rate else 0
    
    foi <- beta_t * I / N
    
    dS <- -foi * S - vax_rate_t * S
    dE <-  foi * S - sigma * E
    dI <-  sigma * E - gamma * I
    dR <-  gamma * I * (1 - ifr_prob)
    dD <-  gamma * I * ifr_prob
    dV <-  vax_rate_t * S
    
    list(c(dS, dE, dI, dR, dD, dV))
  })
}

run_model <- function(N, R0, incubation_days, infectious_days,
                      ifr_pct, hosp_rate_pct,
                      seed_cases, t_max,
                      npi_day, npi_reduction,
                      vax_day, vax_coverage, vax_duration_days) {
  
  beta  <- R0 / infectious_days
  sigma <- 1 / incubation_days
  gamma <- 1 / infectious_days
  
  vax_daily_rate <- if (vax_duration_days > 0) vax_coverage / vax_duration_days else 0
  
  params <- c(
    beta           = beta,
    sigma          = sigma,
    gamma          = gamma,
    ifr_prob       = ifr_pct / 100,
    t_npi          = npi_day,
    npi_eff        = npi_reduction / 100,
    t_vax          = vax_day,
    vax_daily_rate = vax_daily_rate,
    vax_duration   = vax_duration_days
  )
  
  state <- c(S = N - seed_cases, E = 0, I = seed_cases, R = 0, D = 0, V = 0)
  times <- seq(0, t_max, by = 0.5)
  
  out <- deSolve::ode(y = state, times = times, func = seir_ode,
                      parms = params, method = "lsoda")
  
  as.data.frame(out) %>%
    mutate(
      N         = S + E + I + R + D + V,
      incidence = pmax(0, -c(0, diff(S + V))),
      daily_inc = incidence * 2,
      hosp      = I * (hosp_rate_pct / 100),
      Reff      = R0 * (S / N) * if_else(time >= npi_day, 1 - npi_reduction / 100, 1)
    )
}

# ==============================================================================
# SUMMARY METRICS
# ==============================================================================

compute_metrics <- function(df, N, R0, vax_coverage, npi_reduction) {
  peak_I      <- max(df$I, na.rm = TRUE)
  peak_day    <- df$time[which.max(df$I)]
  total_inf   <- max(df$R + df$D, na.rm = TRUE)
  total_dead  <- max(df$D, na.rm = TRUE)
  attack_rate <- total_inf / N * 100
  peak_hosp   <- max(df$hosp, na.rm = TRUE)
  herd_thresh <- (1 - 1/R0) * 100
  epi_end     <- df$time[max(which(df$I > 1))]
  
  list(
    peak_I      = round(peak_I),
    peak_day    = round(peak_day),
    total_inf   = round(total_inf),
    total_dead  = round(total_dead),
    attack_rate = round(attack_rate, 1),
    peak_hosp   = round(peak_hosp),
    herd_thresh = round(herd_thresh, 1),
    epi_end     = round(epi_end)
  )
}

# ==============================================================================
# COLOUR PALETTE
# ==============================================================================

COLS <- list(
  S = "#4e9af1", E = "#e8c33a", I = "#e8523a", R = "#4bc98a",
  D = "#c94b4b", V = "#b06de8", hosp = "#ff8c42", Reff = "#ffffff",
  bg = "#0b0f14", card = "#111820", border = "#1e2d3d",
  text = "#cdd9e5", muted = "#637080", accent = "#4e9af1"
)

# ==============================================================================
# UI
# ==============================================================================

ui <- page_navbar(
  title = tags$span(
    style = "font-family: 'DM Mono', monospace; letter-spacing: 3px; font-size: 15px;
             color: #4e9af1; font-weight: 500;",
    "\u2b21 SEIR EPIDEMIC SIMULATOR"
  ),
  window_title = "SEIR Epidemic Simulator",
  
  theme = bs_theme(
    version = 5, bg = "#0b0f14", fg = "#cdd9e5",
    primary = "#4e9af1", secondary = "#1e2d3d",
    success = "#4bc98a", danger = "#e8523a", warning = "#e8c33a", info = "#b06de8",
    base_font    = font_google("Barlow"),
    heading_font = font_google("Barlow Condensed"),
    code_font    = font_google("DM Mono"),
    "navbar-bg"          = "#080c10",
    "card-bg"            = "#111820",
    "card-border-color"  = "#1e2d3d",
    "input-bg"           = "#0d1520",
    "input-border-color" = "#1e2d3d",
    "input-color"        = "#cdd9e5",
    "form-label-color"   = "#637080"
  ),
  
  header = tags$head(tags$style(HTML('
    @import url("https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=Barlow:wght@300;400;500;600&family=Barlow+Condensed:wght@300;400;500;600;700&display=swap");
    body { background: #0b0f14 !important; }
    .nav-tabs .nav-link { color: #637080 !important; border-color: transparent !important; font-family: "Barlow Condensed"; letter-spacing: 1px; font-size: 13px; }
    .nav-tabs .nav-link.active { color: #4e9af1 !important; border-bottom: 2px solid #4e9af1 !important; background: transparent !important; }
    .kpi-box { background: #111820; border: 1px solid #1e2d3d; border-radius: 4px; padding: 14px 18px; position: relative; overflow: hidden; }
    .kpi-box::before { content: ""; position: absolute; top: 0; left: 0; width: 3px; height: 100%; background: var(--kpi-accent, #4e9af1); }
    .kpi-label { font-family: "DM Mono", monospace; font-size: 10px; letter-spacing: 2px; text-transform: uppercase; color: #637080; margin-bottom: 4px; }
    .kpi-value { font-family: "Barlow Condensed", sans-serif; font-size: 26px; font-weight: 600; color: #cdd9e5; line-height: 1; }
    .kpi-sub { font-size: 11px; color: #637080; margin-top: 3px; font-family: "DM Mono", monospace; }
    .sidebar-section-label { font-family: "DM Mono", monospace; font-size: 10px; letter-spacing: 2px; text-transform: uppercase; color: #4e9af1; margin: 16px 0 8px 0; border-bottom: 1px solid #1e2d3d; padding-bottom: 5px; }
    .form-range::-webkit-slider-thumb { background: #4e9af1; }
    .form-range::-webkit-slider-runnable-track { background: #1e2d3d; }
    .selectize-control .selectize-input { background: #0d1520 !important; border-color: #1e2d3d !important; color: #cdd9e5 !important; font-family: "Barlow", sans-serif; font-size: 13px; }
    .selectize-control .selectize-dropdown { background: #111820 !important; border-color: #1e2d3d !important; color: #cdd9e5 !important; }
    .selectize-control .selectize-dropdown-content .option:hover, .selectize-control .selectize-dropdown-content .option.active { background: #1e2d3d !important; }
    .card { border-radius: 4px !important; }
    .card-header { font-family: "Barlow Condensed"; letter-spacing: 1px; font-size: 13px; text-transform: uppercase; }
    .intervention-tag { display: inline-block; font-family: "DM Mono", monospace; font-size: 10px; padding: 2px 7px; border-radius: 3px; margin: 2px; letter-spacing: 1px; }
    .scenario-badge { font-family: "DM Mono"; font-size: 10px; letter-spacing: 1px; background: #1e2d3d; border-radius: 3px; padding: 2px 8px; color: #4e9af1; text-transform: uppercase; }
    ::-webkit-scrollbar { width: 5px; height: 5px; }
    ::-webkit-scrollbar-track { background: #0b0f14; }
    ::-webkit-scrollbar-thumb { background: #1e2d3d; border-radius: 3px; }
    .shiny-input-container label { font-size: 12px; }
    .irs--shiny .irs-bar { background: #4e9af1; border-color: #4e9af1; }
    .irs--shiny .irs-handle { background: #4e9af1; border-color: #4e9af1; }
    .irs--shiny .irs-single { background: #1e2d3d; font-family: "DM Mono"; font-size: 11px; }
    .irs--shiny .irs-from, .irs--shiny .irs-to { background: #1e2d3d; }
    .irs--shiny .irs-line { background: #1e2d3d; }
    .irs--shiny .irs-grid-text { color: #637080; font-size: 9px; }
    .download-btn { background: transparent !important; border: 1px solid #1e2d3d !important; color: #637080 !important; font-family: "DM Mono", monospace !important; font-size: 11px !important; letter-spacing: 1px !important; padding: 5px 12px !important; border-radius: 3px !important; transition: all 0.2s; }
    .download-btn:hover { border-color: #4e9af1 !important; color: #4e9af1 !important; }
  '))),
  
  # ============================================================
  # TAB 1: Simulator
  # ============================================================
  nav_panel("SIMULATOR",
            layout_sidebar(
              sidebar = sidebar(
                width = 290, bg = "#080c10",
                style = "border-right: 1px solid #1e2d3d; overflow-y: auto; height: 100%;",
                
                tags$div(class = "sidebar-section-label", "DISEASE PRESET"),
                selectInput("preset", NULL, choices = names(PRESETS), selected = "COVID-19 (Ancestral)"),
                
                tags$div(class = "sidebar-section-label", "POPULATION"),
                sliderInput("pop_N",      "Population size (N)",         10000, 10000000, 1000000, 10000, pre = "", sep = ","),
                sliderInput("seed_cases", "Initial infectious cases",    1, 500, 10, 1),
                
                tags$div(class = "sidebar-section-label", "PATHOGEN"),
                sliderInput("R0",         "Basic reproduction number (R\u2080)", 0.5, 20,   2.8, 0.1),
                sliderInput("incubation", "Incubation period (days)",            1,   21,   5,   0.5),
                sliderInput("infectious", "Infectious period (days)",            1,   21,   5,   0.5),
                sliderInput("ifr",        "Infection fatality rate (%)",         0.01, 60,  1.2, 0.01),
                sliderInput("hosp_rate",  "Hospitalisation rate (%)",            0.1,  90,  8,   0.1),
                
                sliderInput("t_max", "Simulation length (days)", 60, 730, 365, 10),
                
                tags$div(class = "sidebar-section-label", "NON-PHARMACEUTICAL INTERVENTION"),
                sliderInput("npi_day", "Implement NPI on day",           1, 365, 60, 1),
                sliderInput("npi_eff", "NPI effectiveness (% \u03b2 reduction)", 0, 90, 0, 1),
                
                tags$div(class = "sidebar-section-label", "VACCINATION"),
                sliderInput("vax_day",      "Start vaccination on day",         1, 365, 90,  1),
                sliderInput("vax_coverage", "Final coverage (% of susceptibles)", 0, 100, 0, 1),
                sliderInput("vax_duration", "Rollout duration (days)",          7, 365, 90,  7),
                
                tags$hr(style = "border-color: #1e2d3d; margin: 16px 0;"),
                downloadButton("dl_data", "\u2b07 EXPORT CSV", class = "download-btn", style = "width:100%;")
              ),
              
              layout_columns(
                col_widths = c(12),
                uiOutput("kpi_row"),
                
                card(
                  full_screen = TRUE,
                  style = "background:#111820; border:1px solid #1e2d3d;",
                  card_header(
                    layout_columns(
                      col_widths = c(6, 6),
                      tags$span("EPIDEMIC CURVE", style = "color:#4e9af1;"),
                      tags$div(style = "text-align:right;", uiOutput("intervention_tags"))
                    ),
                    style = "background:#0d1520; border-bottom:1px solid #1e2d3d;"
                  ),
                  plotlyOutput("epi_curve", height = "340px")
                ),
                
                layout_columns(
                  col_widths = c(6, 6),
                  card(style = "background:#111820; border:1px solid #1e2d3d;",
                       card_header("COMPARTMENTS OVER TIME",
                                   style = "background:#0d1520; border-bottom:1px solid #1e2d3d; color:#4e9af1;"),
                       plotlyOutput("compartment_plot", height = "280px")),
                  card(style = "background:#111820; border:1px solid #1e2d3d;",
                       card_header("EFFECTIVE REPRODUCTION NUMBER (R\u1eefff)",
                                   style = "background:#0d1520; border-bottom:1px solid #1e2d3d; color:#4e9af1;"),
                       plotlyOutput("reff_plot", height = "280px"))
                ),
                
                layout_columns(
                  col_widths = c(6, 6),
                  card(style = "background:#111820; border:1px solid #1e2d3d;",
                       card_header("DAILY HOSPITALISATION BURDEN",
                                   style = "background:#0d1520; border-bottom:1px solid #1e2d3d; color:#4e9af1;"),
                       plotlyOutput("hosp_plot", height = "240px")),
                  card(style = "background:#111820; border:1px solid #1e2d3d;",
                       card_header("FINAL STATE BREAKDOWN",
                                   style = "background:#0d1520; border-bottom:1px solid #1e2d3d; color:#4e9af1;"),
                       plotlyOutput("final_pie", height = "240px"))
                )
              )
            )
  ),
  
  # ============================================================
  # TAB 2: Scenario Comparison
  # ============================================================
  nav_panel("SCENARIO COMPARE",
            layout_columns(
              col_widths = c(12),
              
              tags$div(
                style = "padding: 8px 0 0 0; color:#637080; font-family:'DM Mono'; font-size:11px; letter-spacing:1px;",
                "COMPARING: current sidebar settings vs. counterfactual scenarios below"
              ),
              
              layout_columns(
                col_widths = c(4, 4, 4),
                
                card(style = "background:#111820; border:1px solid #1e2d3d;",
                     card_header(tags$span(class = "scenario-badge", "SCENARIO B"),
                                 style = "background:#0d1520; border-bottom:1px solid #1e2d3d;"),
                     sliderInput("b_npi_eff", "NPI effectiveness (%)", 0, 90,  40, 1),
                     sliderInput("b_vax_cov", "Vaccine coverage (%)",  0, 100, 50, 1),
                     sliderInput("b_npi_day", "NPI start day",         1, 200, 30, 1)),
                
                card(style = "background:#111820; border:1px solid #1e2d3d;",
                     card_header(tags$span(class = "scenario-badge", style = "color:#4bc98a;", "SCENARIO C"),
                                 style = "background:#0d1520; border-bottom:1px solid #1e2d3d;"),
                     sliderInput("c_npi_eff", "NPI effectiveness (%)", 0, 90,  70, 1),
                     sliderInput("c_vax_cov", "Vaccine coverage (%)",  0, 100, 80, 1),
                     sliderInput("c_npi_day", "NPI start day",         1, 200, 15, 1)),
                
                card(style = "background:#111820; border:1px solid #1e2d3d;",
                     card_header(tags$span(class = "scenario-badge", style = "color:#e8523a;",
                                           "SCENARIO D \u2014 NO INTERVENTION"),
                                 style = "background:#0d1520; border-bottom:1px solid #1e2d3d;"),
                     tags$p(style = "color:#637080; font-size:12px; padding:8px; font-family:'DM Mono';",
                            "Unmitigated epidemic \u2014 all parameters from sidebar, NPI = 0%, Vax = 0%."))
              ),
              
              card(full_screen = TRUE, style = "background:#111820; border:1px solid #1e2d3d;",
                   card_header("PEAK INFECTIOUS COMPARISON",
                               style = "background:#0d1520; border-bottom:1px solid #1e2d3d; color:#4e9af1;"),
                   plotlyOutput("compare_plot", height = "360px")),
              
              card(style = "background:#111820; border:1px solid #1e2d3d;",
                   card_header("SCENARIO METRICS TABLE",
                               style = "background:#0d1520; border-bottom:1px solid #1e2d3d; color:#4e9af1;"),
                   DTOutput("compare_table"))
            )
  ),
  
  # ============================================================
  # TAB 3: Methods
  # ============================================================
  nav_panel("METHODS",
            tags$div(
              style = "max-width:760px; margin:32px auto; font-family:'Barlow'; line-height:1.7; color:#cdd9e5;",
              tags$h2("Model Description",
                      style = "font-family:'Barlow Condensed'; color:#4e9af1; font-weight:600; letter-spacing:2px;"),
              tags$p("This app implements a deterministic SEIVRD compartmental model solved as a system of ODEs using the ",
                     tags$code("deSolve"), " package."),
              
              tags$h4("Compartments", style = "font-family:'Barlow Condensed'; color:#cdd9e5; margin-top:24px;"),
              tags$ul(
                tags$li(tags$strong("S"), " \u2014 Susceptible"),
                tags$li(tags$strong("E"), " \u2014 Exposed (infected, pre-infectious)"),
                tags$li(tags$strong("I"), " \u2014 Infectious"),
                tags$li(tags$strong("R"), " \u2014 Recovered (immune)"),
                tags$li(tags$strong("D"), " \u2014 Dead (disease-specific mortality)"),
                tags$li(tags$strong("V"), " \u2014 Vaccinated (moved from S, fully protected)")
              ),
              
              tags$h4("Equations", style = "font-family:'Barlow Condensed'; color:#cdd9e5; margin-top:24px;"),
              tags$div(
                style = "background:#0d1520; border:1px solid #1e2d3d; padding:16px; border-radius:4px;
                 font-family:'DM Mono'; font-size:13px; color:#4bc98a;",
                tags$p("dS/dt = \u2212\u03b2(t) \u00b7 I/N \u00b7 S \u2212 \u03bd(t) \u00b7 S"),
                tags$p("dE/dt =  \u03b2(t) \u00b7 I/N \u00b7 S \u2212 \u03c3 \u00b7 E"),
                tags$p("dI/dt =  \u03c3 \u00b7 E \u2212 \u03b3 \u00b7 I"),
                tags$p("dR/dt =  \u03b3 \u00b7 I \u00b7 (1 \u2212 IFR)"),
                tags$p("dD/dt =  \u03b3 \u00b7 I \u00b7 IFR"),
                tags$p("dV/dt =  \u03bd(t) \u00b7 S")
              ),
              
              tags$h4("Key Parameters", style = "font-family:'Barlow Condensed'; color:#cdd9e5; margin-top:24px;"),
              tags$table(class = "table", style = "font-size:13px;",
                         tags$thead(tags$tr(tags$th("Symbol"), tags$th("Meaning"), tags$th("Derived from"))),
                         tags$tbody(
                           tags$tr(tags$td(tags$code("\u03b2")),    tags$td("Transmission rate"),           tags$td("R\u2080 / infectious period")),
                           tags$tr(tags$td(tags$code("\u03c3")),    tags$td("Rate of becoming infectious"), tags$td("1 / incubation period")),
                           tags$tr(tags$td(tags$code("\u03b3")),    tags$td("Recovery rate"),               tags$td("1 / infectious period")),
                           tags$tr(tags$td(tags$code("\u03b2(t)")), tags$td("Time-varying transmission"),   tags$td("\u03b2 \u00b7 (1 \u2212 NPI%) after NPI day")),
                           tags$tr(tags$td(tags$code("\u03bd(t)")), tags$td("Vaccination rate"),            tags$td("Coverage% \u00f7 rollout days, during campaign"))
                         )
              ),
              
              tags$h4("Assumptions & Limitations",
                      style = "font-family:'Barlow Condensed'; color:#cdd9e5; margin-top:24px;"),
              tags$ul(
                tags$li("Homogeneous mixing \u2014 all individuals have equal contact rates."),
                tags$li("No age structure, risk groups, or spatial heterogeneity."),
                tags$li("Vaccination provides immediate, complete protection."),
                tags$li("No waning immunity or reinfection."),
                tags$li("NPI effect is instantaneous and constant after implementation."),
                tags$li("IFR and hospitalisation rate are constant across groups and time.")
              ),
              
              tags$h4("References", style = "font-family:'Barlow Condensed'; color:#cdd9e5; margin-top:24px;"),
              tags$ul(
                tags$li("Kermack & McKendrick (1927). A contribution to the mathematical theory of epidemics. Proc. R. Soc. A."),
                tags$li("Soetaert, Petzoldt & Setzer (2010). Solving differential equations in R. R Journal."),
                tags$li("Keeling & Rohani (2008). Modeling Infectious Diseases in Humans and Animals. Princeton UP.")
              )
            )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  observeEvent(input$preset, {
    p <- PRESETS[[input$preset]]
    if (!is.null(p)) {
      updateSliderInput(session, "R0",         value = p$R0)
      updateSliderInput(session, "incubation", value = p$incubation)
      updateSliderInput(session, "infectious", value = p$infectious)
      updateSliderInput(session, "ifr",        value = p$ifr)
      updateSliderInput(session, "hosp_rate",  value = p$hosp_rate)
    }
  })
  
  model_data <- reactive({
    run_model(
      N                 = input$pop_N,
      R0                = input$R0,
      incubation_days   = input$incubation,
      infectious_days   = input$infectious,
      ifr_pct           = input$ifr,
      hosp_rate_pct     = input$hosp_rate,
      seed_cases        = input$seed_cases,
      t_max             = input$t_max,
      npi_day           = input$npi_day,
      npi_reduction     = input$npi_eff,
      vax_day           = input$vax_day,
      vax_coverage      = input$vax_coverage / 100,
      vax_duration_days = input$vax_duration
    )
  })
  
  metrics <- reactive({
    compute_metrics(model_data(), input$pop_N, input$R0,
                    input$vax_coverage, input$npi_eff)
  })
  
  # --------------------------------------------------------------------------
  # KPI Row
  # --------------------------------------------------------------------------
  
  kpi <- function(label, value, sub, accent) {
    tags$div(
      class = "kpi-box",
      style = paste0("--kpi-accent:", accent, ";"),
      tags$div(class = "kpi-label", label),
      tags$div(class = "kpi-value", value),
      tags$div(class = "kpi-sub",   sub)
    )
  }
  
  output$kpi_row <- renderUI({
    m <- metrics()
    layout_columns(
      col_widths = c(2, 2, 2, 2, 2, 2),
      kpi("PEAK INFECTIOUS",         format(m$peak_I,     big.mark = ","), paste0("day ", m$peak_day), "#e8523a"),
      kpi("TOTAL INFECTIONS",        format(m$total_inf,  big.mark = ","), paste0(m$attack_rate, "% attack rate"), "#e8c33a"),
      kpi("TOTAL DEATHS",            format(m$total_dead, big.mark = ","), paste0("IFR ", input$ifr, "%"), "#c94b4b"),
      kpi("PEAK HOSPITALISATIONS",   format(m$peak_hosp,  big.mark = ","), "concurrent", "#ff8c42"),
      kpi("HERD IMMUNITY THRESHOLD", paste0(m$herd_thresh, "%"),           paste0("R\u2080 = ", input$R0), "#4e9af1"),
      kpi("EPIDEMIC END",            paste0("Day ", m$epi_end),            paste0(m$epi_end, " days duration"), "#4bc98a")
    )
  })
  
  # --------------------------------------------------------------------------
  # Intervention tags
  # --------------------------------------------------------------------------
  
  output$intervention_tags <- renderUI({
    tags_list <- list()
    if (input$npi_eff > 0)
      tags_list[[1]] <- tags$span(class = "intervention-tag",
                                  style = "background:#1a2a1a; color:#4bc98a; border:1px solid #2a4a2a;",
                                  paste0("NPI day ", input$npi_day, " (", input$npi_eff, "% \u2193\u03b2)"))
    if (input$vax_coverage > 0)
      tags_list[[2]] <- tags$span(class = "intervention-tag",
                                  style = "background:#1a1a2a; color:#b06de8; border:1px solid #2a2a4a;",
                                  paste0("VAX day ", input$vax_day, " (", input$vax_coverage, "% cov.)"))
    if (length(tags_list) == 0)
      tags_list[[1]] <- tags$span(class = "intervention-tag",
                                  style = "background:#1e2d3d; color:#637080;", "NO INTERVENTIONS")
    tags$div(style = "text-align:right;", tagList(tags_list))
  })
  
  # --------------------------------------------------------------------------
  # Epidemic curve
  # --------------------------------------------------------------------------
  
  output$epi_curve <- renderPlotly({
    df_plot <- model_data() %>% filter(time %% 1 == 0)
    
    p <- plot_ly(df_plot, x = ~time, type = "scatter", mode = "none") %>%
      add_trace(y = ~I,    name = "Infectious",  fill = "tozeroy",
                fillcolor = "rgba(232,82,58,0.35)",  line = list(color = COLS$I,    width = 2)) %>%
      add_trace(y = ~E,    name = "Exposed",     fill = "tozeroy",
                fillcolor = "rgba(232,195,58,0.25)", line = list(color = COLS$E,    width = 1.5)) %>%
      add_trace(y = ~hosp, name = "Hospitalised", fill = "tozeroy",
                fillcolor = "rgba(255,140,66,0.3)",  line = list(color = COLS$hosp, width = 1.5, dash = "dot"))
    
    shapes <- list(); annotations <- list()
    if (input$npi_eff > 0) {
      shapes[[1]]      <- list(type = "line", x0 = input$npi_day, x1 = input$npi_day,
                               y0 = 0, y1 = 1, yref = "paper",
                               line = list(color = "#4bc98a", width = 1.5, dash = "dash"))
      annotations[[1]] <- list(x = input$npi_day, y = 0.97, yref = "paper", text = "NPI",
                               showarrow = FALSE, font = list(color = "#4bc98a", size = 10, family = "DM Mono"))
    }
    if (input$vax_coverage > 0) {
      shapes[[length(shapes)+1]]      <- list(type = "line", x0 = input$vax_day, x1 = input$vax_day,
                                              y0 = 0, y1 = 1, yref = "paper",
                                              line = list(color = "#b06de8", width = 1.5, dash = "dash"))
      annotations[[length(annotations)+1]] <- list(x = input$vax_day, y = 0.90, yref = "paper", text = "VAX",
                                                   showarrow = FALSE, font = list(color = "#b06de8", size = 10, family = "DM Mono"))
    }
    
    p %>% layout(
      xaxis = list(title = "Day", color = "#637080", gridcolor = "#1e2d3d", zerolinecolor = "#1e2d3d"),
      yaxis = list(title = "Individuals", color = "#637080", gridcolor = "#1e2d3d", zerolinecolor = "#1e2d3d"),
      paper_bgcolor = "#111820", plot_bgcolor = "#111820",
      font = list(color = "#637080", family = "Barlow"),
      legend = list(font = list(color = "#cdd9e5"), bgcolor = "#0d1520", bordercolor = "#1e2d3d", borderwidth = 1),
      shapes = shapes, annotations = annotations,
      hovermode = "x unified", margin = list(l = 50, r = 20, t = 10, b = 40)
    )
  })
  
  # --------------------------------------------------------------------------
  # Compartments
  # --------------------------------------------------------------------------
  
  output$compartment_plot <- renderPlotly({
    df <- model_data() %>% filter(time %% 1 == 0) %>%
      select(time, S, E, I, R, D, V) %>% pivot_longer(-time)
    col_map <- c(S = COLS$S, E = COLS$E, I = COLS$I, R = COLS$R, D = COLS$D, V = COLS$V)
    
    plot_ly(df, x = ~time, y = ~value, color = ~name,
            colors = unname(col_map[unique(df$name)]),
            type = "scatter", mode = "lines", line = list(width = 1.8)) %>%
      layout(
        xaxis = list(title = "Day", color = "#637080", gridcolor = "#1e2d3d"),
        yaxis = list(title = "N",   color = "#637080", gridcolor = "#1e2d3d"),
        paper_bgcolor = "#111820", plot_bgcolor = "#111820",
        font = list(color = "#637080", family = "Barlow"),
        legend = list(font = list(color = "#cdd9e5"), bgcolor = "#0d1520", bordercolor = "#1e2d3d", borderwidth = 1),
        hovermode = "x unified", margin = list(l = 50, r = 20, t = 10, b = 40)
      )
  })
  
  # --------------------------------------------------------------------------
  # Reff plot
  # --------------------------------------------------------------------------
  
  output$reff_plot <- renderPlotly({
    df <- model_data() %>% filter(time %% 1 == 0)
    
    plot_ly(df, x = ~time, y = ~Reff, type = "scatter", mode = "lines",
            line = list(color = "#ffffff", width = 2),
            fill = "tozeroy", fillcolor = "rgba(255,255,255,0.05)", name = "Reff") %>%
      add_trace(x = range(df$time), y = c(1, 1), type = "scatter", mode = "lines",
                line = list(color = "#e8523a", width = 1.5, dash = "dot"),
                name = "Reff = 1", showlegend = TRUE) %>%
      layout(
        xaxis = list(title = "Day",  color = "#637080", gridcolor = "#1e2d3d"),
        yaxis = list(title = "Reff", color = "#637080", gridcolor = "#1e2d3d"),
        paper_bgcolor = "#111820", plot_bgcolor = "#111820",
        font = list(color = "#637080", family = "Barlow"),
        legend = list(font = list(color = "#cdd9e5"), bgcolor = "#0d1520"),
        annotations = list(list(x = max(df$time) * 0.8, y = 1.1,
                                text = "epidemic control threshold", showarrow = FALSE,
                                font = list(color = "#e8523a", size = 10, family = "DM Mono"))),
        margin = list(l = 50, r = 20, t = 10, b = 40)
      )
  })
  
  # --------------------------------------------------------------------------
  # Hospitalisation
  # --------------------------------------------------------------------------
  
  output$hosp_plot <- renderPlotly({
    df            <- model_data() %>% filter(time %% 1 == 0)
    capacity_line <- input$pop_N * 0.002
    
    plot_ly(df, x = ~time, y = ~hosp, type = "scatter", mode = "none",
            fill = "tozeroy", fillcolor = "rgba(255,140,66,0.5)",
            line = list(color = COLS$hosp), name = "Hospitalised") %>%
      add_trace(x = range(df$time), y = c(capacity_line, capacity_line),
                type = "scatter", mode = "lines",
                line = list(color = "#e8523a", dash = "dot", width = 1.5),
                name = "System capacity (0.2%)", fill = "none") %>%
      layout(
        xaxis  = list(title = "Day",         color = "#637080", gridcolor = "#1e2d3d"),
        yaxis  = list(title = "Individuals", color = "#637080", gridcolor = "#1e2d3d"),
        paper_bgcolor = "#111820", plot_bgcolor = "#111820",
        font   = list(color = "#637080", family = "Barlow"),
        legend = list(font = list(color = "#cdd9e5"), bgcolor = "#0d1520", bordercolor = "#1e2d3d"),
        margin = list(l = 50, r = 20, t = 10, b = 40)
      )
  })
  
  # --------------------------------------------------------------------------
  # Final state pie
  # --------------------------------------------------------------------------
  
  output$final_pie <- renderPlotly({
    df_last <- tail(model_data(), 1)
    labels  <- c("Susceptible", "Recovered", "Dead", "Vaccinated")
    values  <- c(df_last$S, df_last$R, df_last$D, df_last$V)
    colors  <- c(COLS$S, COLS$R, COLS$D, COLS$V)
    
    plot_ly(labels = labels, values = values, type = "pie", hole = 0.5,
            marker = list(colors = colors, line = list(color = "#111820", width = 2)),
            textinfo = "label+percent",
            insidetextfont = list(color = "#cdd9e5", size = 11, family = "Barlow"),
            hovertemplate = "<b>%{label}</b><br>N: %{value:,.0f}<br>%{percent}<extra></extra>") %>%
      layout(showlegend = FALSE, paper_bgcolor = "#111820",
             font = list(color = "#637080", family = "Barlow"),
             margin = list(t = 10, b = 10, l = 10, r = 10))
  })
  
  # --------------------------------------------------------------------------
  # Scenario comparison
  # --------------------------------------------------------------------------
  
  run_scenario <- function(npi_eff, vax_cov, npi_day, vax_day = NULL) {
    vd <- if (is.null(vax_day)) input$vax_day else vax_day
    run_model(
      N = input$pop_N, R0 = input$R0,
      incubation_days = input$incubation, infectious_days = input$infectious,
      ifr_pct = input$ifr, hosp_rate_pct = input$hosp_rate,
      seed_cases = input$seed_cases, t_max = input$t_max,
      npi_day = npi_day, npi_reduction = npi_eff,
      vax_day = vd, vax_coverage = vax_cov / 100,
      vax_duration_days = input$vax_duration
    )
  }
  
  scenarios <- reactive({
    list(
      A = model_data(),
      B = run_scenario(input$b_npi_eff, input$b_vax_cov, input$b_npi_day),
      C = run_scenario(input$c_npi_eff, input$c_vax_cov, input$c_npi_day),
      D = run_scenario(0, 0, 999)
    )
  })
  
  output$compare_plot <- renderPlotly({
    sc       <- scenarios()
    cols     <- c(A = "#4e9af1", B = "#e8c33a", C = "#4bc98a", D = "#e8523a")
    names_sc <- c(
      A = paste0("A: Current (NPI ", input$npi_eff, "%, Vax ", input$vax_coverage, "%)"),
      B = paste0("B: NPI ", input$b_npi_eff, "%, Vax ", input$b_vax_cov, "%"),
      C = paste0("C: NPI ", input$c_npi_eff, "%, Vax ", input$c_vax_cov, "%"),
      D = "D: No intervention"
    )
    p <- plot_ly()
    for (s in c("A","B","C","D")) {
      df_s <- sc[[s]] %>% filter(time %% 1 == 0)
      p <- add_trace(p, data = df_s, x = ~time, y = ~I,
                     type = "scatter", mode = "lines",
                     name = names_sc[s], line = list(color = cols[s], width = 2.5))
    }
    p %>% layout(
      xaxis = list(title = "Day", color = "#637080", gridcolor = "#1e2d3d"),
      yaxis = list(title = "Infectious", color = "#637080", gridcolor = "#1e2d3d"),
      paper_bgcolor = "#111820", plot_bgcolor = "#111820",
      font = list(color = "#637080", family = "Barlow"),
      legend = list(font = list(color = "#cdd9e5"), bgcolor = "#0d1520", bordercolor = "#1e2d3d", borderwidth = 1),
      hovermode = "x unified", margin = list(l = 50, r = 20, t = 10, b = 40)
    )
  })
  
  # FIX: renderDT is now properly closed with }) before output$dl_data.
  # Previously the downloadHandler was accidentally nested inside renderDT,
  # causing the column-name error (the handler object's names were used instead
  # of df_table's names).
  output$compare_table <- renderDT({
    sc <- scenarios()
    N  <- input$pop_N
    
    rows <- purrr::imap(sc, function(df, name) {
      m <- compute_metrics(df, N, input$R0, 0, 0)
      tibble(
        Scenario         = switch(name,
                                  A = "A \u2014 Current",
                                  B = paste0("B \u2014 NPI ", input$b_npi_eff, "%, Vax ", input$b_vax_cov, "%"),
                                  C = paste0("C \u2014 NPI ", input$c_npi_eff, "%, Vax ", input$c_vax_cov, "%"),
                                  D = "D \u2014 No intervention"),
        `Peak I`          = format(m$peak_I,     big.mark = ","),
        `Peak day`        = m$peak_day,
        `Total infected`  = format(m$total_inf,  big.mark = ","),
        `Attack rate %`   = m$attack_rate,
        `Deaths`          = format(m$total_dead, big.mark = ","),
        `Peak hosp.`      = format(m$peak_hosp,  big.mark = ","),
        `Duration (days)` = m$epi_end
      )
    })
    
    df_table <- bind_rows(rows)
    
    datatable(
      df_table,
      rownames = FALSE,
      options  = list(
        dom        = "t",
        pageLength = 4,
        initComplete = JS("function(s,j){$(this.api().table().header()).css({'background':'#0d1520','color':'#4e9af1','font-family':'DM Mono','font-size':'11px','letter-spacing':'1px'});}")
      ),
      class = "cell-border"
    ) %>%
      formatStyle(
        columns         = names(df_table),  # now correctly refers to df_table's column names
        backgroundColor = "#111820",
        color           = "#cdd9e5",
        fontFamily      = "Barlow",
        fontSize        = "13px"
      ) %>%
      formatStyle(
        "Attack rate %",
        color      = styleInterval(c(30, 60), c("#4bc98a", "#e8c33a", "#e8523a")),
        fontWeight = "bold"
      )
  })  # <-- closes renderDT
  
  # --------------------------------------------------------------------------
  # CSV export  (correctly outside renderDT, inside server)
  # --------------------------------------------------------------------------
  
  output$dl_data <- downloadHandler(
    filename = function() paste0("seir_simulation_", Sys.Date(), ".csv"),
    content  = function(file) {
      model_data() %>%
        filter(time %% 1 == 0) %>%
        mutate(across(where(is.numeric), ~round(.x, 2))) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
}  # <-- closes server function

# ==============================================================================
# RUN
# ==============================================================================

shinyApp(ui = ui, server = server)


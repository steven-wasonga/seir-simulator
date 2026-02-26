# ⬡ SEIR Epidemic Simulator

An interactive, single-file R Shiny app for modelling infectious disease outbreaks. Built for public health exploration and epidemiology education

![R](https://img.shields.io/badge/R-%3E%3D4.0-276DC3?style=flat&logo=r&logoColor=white)
![Shiny](https://img.shields.io/badge/Shiny-app-4e9af1?style=flat)
![License](https://img.shields.io/badge/license-MIT-green?style=flat)

---

## Features

- **SEIVRD compartmental model** — Susceptible, Exposed, Infectious, Vaccinated, Recovered, Dead
- **8 disease presets** — COVID-19 (Ancestral & Omicron), Seasonal Flu, 1918 H1N1, Measles, Ebola, SARS-CoV-1, Mpox, plus a fully customisable Custom mode
- **NPI modelling** — simulate non-pharmaceutical interventions (mask mandates, lockdowns, etc.) with a configurable start day and effectiveness
- **Vaccination rollout** — set coverage target, start day, and rollout duration
- **Scenario comparison** — run up to 4 scenarios side-by-side (including an unmitigated baseline) to compare peak infections, deaths, and attack rates
- **Live KPI dashboard** — peak infectious, total infections, total deaths, peak hospitalisations, herd immunity threshold, epidemic end day
- **Interactive charts** — epidemic curve, compartment trajectories, Rₑff over time, hospitalisation burden vs. capacity, final state pie chart
- **CSV export** — download the full simulation output

---

## Screenshots

<img width="1920" height="1080" alt="Screenshot (34)" src="https://github.com/user-attachments/assets/365ceb8d-dde6-4a6e-9a72-4f6774faed03" />
<img width="1920" height="1080" alt="Screenshot (33)" src="https://github.com/user-attachments/assets/9aac6956-57e5-4b1d-90ed-e59bbc4c3cd1" />
<img width="1920" height="1080" alt="Screenshot (36)" src="https://github.com/user-attachments/assets/00e88cc8-4989-4619-a934-151f68896de0" />
<img width="1920" height="1080" alt="Screenshot (35)" src="https://github.com/user-attachments/assets/958e8b2a-6276-411f-b9bb-200564cc0a26" />


## Installation

**1. Clone the repository**
```bash
git clone https://github.com/yourusername/seir-simulator.git
cd seir-simulator
```

**2. Install R dependencies** (run once inside R)
```r
install.packages(c("shiny", "bslib", "plotly", "deSolve", "dplyr", "tidyr", "scales", "DT", "purrr"))
```

⚠️ deSolve is not available on CRAN and must be installed manually from the archived version.
Download the latest archived .tar.gz from the CRAN archive:
https://cran.r-project.org/src/contrib/Archive/deSolve/
Then install it in R pointing to the downloaded file:
```r
rinstall.packages("/path/to/deSolve_x.x.x.tar.gz", repos = NULL, type = "source")
```
**3. Run the app**
```r
shiny::runApp("seir_app_fixed.R")
```

Or open `seir_app_fixed.R` in RStudio and click **Run App**.

---

## Model Description

The app solves a deterministic **SEIVRD** ODE system using the `deSolve` package.

### Compartments

| Symbol | Meaning |
|--------|---------|
| **S** | Susceptible |
| **E** | Exposed (infected, not yet infectious) |
| **I** | Infectious |
| **V** | Vaccinated (fully protected, moved from S) |
| **R** | Recovered (immune) |
| **D** | Dead (disease-specific mortality) |

### Equations

```
dS/dt = −β(t) · I/N · S − ν(t) · S
dE/dt =  β(t) · I/N · S − σ · E
dI/dt =  σ · E − γ · I
dR/dt =  γ · I · (1 − IFR)
dD/dt =  γ · I · IFR
dV/dt =  ν(t) · S
```

### Key Parameters

| Symbol | Meaning | Derived from |
|--------|---------|--------------|
| β | Transmission rate | R₀ / infectious period |
| σ | Rate of becoming infectious | 1 / incubation period |
| γ | Recovery rate | 1 / infectious period |
| β(t) | Time-varying transmission | β · (1 − NPI%) after NPI day |
| ν(t) | Vaccination rate | Coverage% ÷ rollout days, during campaign |

### Assumptions & Limitations

- Homogeneous mixing — all individuals have equal contact rates
- No age structure, risk groups, or spatial heterogeneity
- Vaccination provides immediate, complete protection
- No waning immunity or reinfection
- NPI effect is instantaneous and constant after implementation
- IFR and hospitalisation rate are constant across groups and time

---

## Disease Presets

| Disease | R₀ | Incubation (days) | Infectious (days) | IFR (%) |
|---|---|---|---|---|
| COVID-19 (Ancestral) | 2.8 | 5 | 5 | 1.2 |
| COVID-19 (Omicron) | 10.0 | 3 | 4 | 0.3 |
| Influenza (Seasonal) | 1.3 | 2 | 4 | 0.1 |
| Influenza (1918 H1N1) | 2.2 | 2 | 4 | 2.5 |
| Measles | 15.0 | 12 | 8 | 0.2 |
| Ebola (West Africa) | 1.8 | 9 | 7 | 50.0 |
| SARS-CoV-1 | 3.0 | 5 | 7 | 9.6 |
| Mpox | 2.1 | 10 | 14 | 1.0 |

---

## Dependencies

| Package | Purpose |
|---------|---------|
| `shiny` | App framework |
| `bslib` | Bootstrap 5 theming |
| `plotly` | Interactive charts |
| `deSolve` | ODE solver (lsoda) |
| `dplyr` / `tidyr` | Data wrangling |
| `scales` | Number formatting |
| `DT` | Interactive data tables |
| `purrr` | Scenario iteration |

---

## References

- Kermack & McKendrick (1927). *A contribution to the mathematical theory of epidemics.* Proc. R. Soc. A.
- Soetaert, Petzoldt & Setzer (2010). *Solving differential equations in R.* The R Journal.
- Keeling & Rohani (2008). *Modeling Infectious Diseases in Humans and Animals.* Princeton University Press.

---

## License

MIT License. See `LICENSE` for details.

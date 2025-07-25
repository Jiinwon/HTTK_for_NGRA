#!/usr/bin/env Rscript
# =====================================================================
# HTTK tutorial – BPA (CAS 80-05-7) full pipeline with non‐uniform TK distribution
# 1) load packages
# 2) identify chemical
# 3) parameterize & adjust
# 4) run 1-compartment & PBTK models
# 5) Monte Carlo Css: mg/day scenario
# 6) Monte Carlo Css: manual TK parameter sampling
# 7) IVIVE (95% AED)
# 8) summary & histogram + density
# =====================================================================

library(httk)
library(ggplot2)

# 1) identify chemical
chem.cas <- get_chem_id(chem.cas = "80-05-7")$chem.cas

# 2) parameterize & adjust
use_pbtk <- TRUE
params <- tryCatch(
  parameterize_pbtk(chem.cas = chem.cas, species = "Human", suppress.messages = TRUE),
  error = function(e) {
    message("parameterize_pbtk failed: ", e$message)
    use_pbtk <<- FALSE
    parameterize_1comp(chem.cas = chem.cas, species = "Human", suppress.messages = TRUE)
  }
)
if (is.na(params$Funbound.plasma) || params$Funbound.plasma == 0) params$Funbound.plasma <- 0.03
if (is.na(params$Clint)           || params$Clint == 0)           params$Clint           <- 1.2

# 3-A) 1-compartment model
sim1 <- solve_1comp(dose = 10, interval = 24, n.doses = 5, params = params, chem.cas = chem.cas)
plot(sim1, main = "1-Compartment Plasma Concentration")

# 3-B) PBTK model
if (use_pbtk) {
  sim_pbtk <- solve_pbtk(
    parameters    = params,
    days          = 7,
    daily.dose    = 10,
    doses.per.day = 1,
    tsteps        = 2,
    output.units  = "mg/L"
  )
  plot(sim_pbtk, which = "Cplasma", main = "PBTK Plasma Concentration")
} else {
  message("Skipping PBTK simulation due to insufficient parameters.")
}

# 5) Monte Carlo Css: mg/day scenario
pop       <- httkpop_generate(method = "vi", nsamp = 1000)
weights   <- pop$BW
css_mgday <- vapply(weights, function(w) {
  tmp    <- params; tmp$BW <- w
  sim    <- solve_1comp(dose = 10, interval = 24, n.doses = 1, params = tmp, chem.cas = chem.cas)
  max(sim[ , "Ccompartment"])
}, numeric(1))
MW        <- 228.29                          # g/mol
css_mgday <- css_mgday * (MW / 1000)         # µM → mg/L
med_mgday <- median(css_mgday)

# 6) Monte Carlo Css: manual TK parameter sampling
set.seed(123)
n         <- 1000
# sample fu and Clint with 30% CV
fu_samp   <- rnorm(n, mean = params$Funbound.plasma, sd = params$Funbound.plasma * 0.3)
cl_samp   <- rnorm(n, mean = params$Clint,            sd = params$Clint * 0.3)
css_tk    <- mapply(function(fu, cl) {
  tmp <- params
  tmp$Funbound.plasma <- max(fu, 1e-6)
  tmp$Clint           <- max(cl, 1e-6)
  sim <- solve_1comp(dose = 10, interval = 24, n.doses = 1, params = tmp, chem.cas = chem.cas)
  max(sim[ , "Ccompartment"])
}, fu_samp, cl_samp)
css_tk    <- css_tk * (MW / 1000)
med_tk    <- median(css_tk)

# 7) IVIVE – 95% AED
aed95 <- tryCatch(
  calc_mc_oral_equiv(
    chem.cas       = chem.cas,
    conc           = 10,
    samples        = 1000,
    which.quantile = 0.95
  ),
  error = function(e) {
    message("calc_mc_oral_equiv failed: ", e$message)
    NA
  }
)

# 8) summary & histogram + density
aed_text <- if (!is.na(aed95)) round(aed95, 3) else "NA"
cat("\n===== Results Summary =====\n",
    "* Css mg/day median:       ", round(med_mgday, 3), "mg/L\n",
    "* Css TK sampling median:  ", round(med_tk,     3), "mg/L\n",
    "* 95% AED:                 ", aed_text,        "mg/kg/day\n")

df <- data.frame(
  Css      = c(css_mgday, css_tk),
  Scenario = rep(c("mg/day", "TK sampling"), each = n)
)

ggplot(df, aes(x = Css, fill = Scenario)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 50) +
  geom_density(size = 1) +
  geom_vline(xintercept = c(med_mgday, med_tk), linetype = "dashed", color = c("blue", "red")) +
  labs(title = "Css Distribution Comparison – BPA",
       x     = "Css (mg/L)",
       y     = "Density") +
  theme_minimal()
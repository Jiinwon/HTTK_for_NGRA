#!/usr/bin/env Rscript
# =====================================================================
# HTTK 실습 – Bisphenol A(BPA, CAS 80-05-7)
#   1) 식별자  → 2) PBTK 파라미터  → 3) 1C·PBTK 시뮬(run_tk 별칭)
#   4) 몬테카로 Css → 5) IVIVE 95 퍼센타일 AED  → 6) 결과 요약
# =====================================================================

## 0. 패키지 준비 ------------------------------------------------------------
if (!requireNamespace("httk", quietly = TRUE))
  install.packages("httk")
library(httk)
library(ggplot2)

## ── solve_model() 새 이름으로 래핑 ───────────────────────────────────────────
run_tk <- function(...) httk::solve_model(...)

## 1. 화학물질 식별 -----------------------------------------------------------
info      <- get_chem_id(chem.cas = "80-05-7")   # BPA
chem.cas  <- info$chem.cas                       # "80-05-7"

## 2. PBTK 파라미터 -----------------------------------------------------------
params <- parameterize_pbtk(
  chem.cas = chem.cas,
  species  = "Human",
  suppress.messages = TRUE
)

# 결측 파라미터 간단 보정(예시)
if (is.na(params$Funbound.plasma) || params$Funbound.plasma == 0)
  params$Funbound.plasma <- 0.03
if (is.na(params$Clint) || params$Clint == 0)
  params$Clint <- 1.2

## 3-A. 1-컴파트 모델 ---------------------------------------------------------
sim_1c <- solve_1comp(
  dose     = 10,      # mg/kg
  interval = 24,      # h
  n.doses  = 5,
  params   = params,
  chem.cas = chem.cas # solve_1comp 내부 호출용
)
plot(sim_1c, main = "1-Compartment Plasma conc.")

## 3-B. PBTK 모델(run_tk 별칭) -------------------------------------------------
sim_pbtk <- run_tk(
  model        = "pbtk",
  parameters   = params,
  dose         = 10,
  dosing.type  = "oral",
  days         = 7,
  delta.t      = 0.5
)
plot(sim_pbtk, which = "Cplasma", main = "PBTK Plasma conc. (run_tk)")

## 4. 몬테카로 Css(1 000명) ----------------------------------------------------
## Monte-Carlo Css (1 000명, 인구 + TK 파라미터 변동)
css <- calc_mc_css(
  chem.cas    = chem.cas,
  samples     = 1000,    # 몬테카로 인구 수
  httkpop     = TRUE,    # 체중·생리학 변동 반영
  vary.params = TRUE,    # TK 파라미터(fu, Clint 등) 분포 사용
  daily.dose  = 10       # ← mg kg⁻¹ day⁻¹ 단위
)
css_median <- median(css)

## 5. IVIVE – 95 퍼센타일 AED --------------------------------------------------
invitro_AC50 <- 10   # µM
aed95 <- calc_mc_oral_equiv(
  chem.cas       = chem.cas,
  conc           = invitro_AC50,
  samples        = 1000,
  which.quantile = 0.95
)

## 6. 결과 요약 ---------------------------------------------------------------
cat("\n===== 결과 요약 =====\n",
    "* Funbound.plasma :", params$Funbound.plasma,
    "\n* Clint           :", params$Clint,
    "\n* Css 중앙값      :", round(css_median, 3), "mg/L",
    "\n* 95% AED         :", round(aed95, 3), "mg/kg/day\n")

## 7. Css 분포 히스토그램(선택) ----------------------------------------------
ggplot(data.frame(css), aes(css)) +
  geom_histogram(bins = 40, fill = "steelblue4") +
  geom_vline(xintercept = css_median, linetype = "dashed") +
  labs(title = "Monte-Carlo Css distribution – BPA",
       x = "Css (mg/L)", y = "Count")
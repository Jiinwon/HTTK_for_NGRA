# 1) HTTK 로드 및 버전 체크
library(httk)
stopifnot(packageVersion("httk") >= "2.7.0")  # HTTK 2.7 이상 필수

# 2) 재현성 확보
set.seed(20250725)

# 3) 계산 대상 CAS 리스트
chem.cas <- c("101-86-0", "106-22-9", "127-51-5", "91-64-5")

# 4) 95th Percentile Css 계산 함수
calc_css95_uM <- function(cas) {
  # A) PBTK 파라미터 초기화
  params <- parameterize_pbtk(
    chem.cas             = cas,
    species              = "Human",
    default.to.human     = TRUE,        # 결손 ADME → 인간 평균 대체
    force.human.clint.fub= TRUE,        # 인간 Clint·fup 강제 사용
    class.exclude        = FALSE,
    physchem.exclude     = FALSE,
    suppress.messages    = TRUE
  )

  # B) Monte Carlo 샘플링
  mc.params <- create_mc_samples(
    parameters = params,
    samples    = 1000,                 # 1000개 샘플 생성
    model      = "pbtk",
    suppress.messages = TRUE
  )

  # C) Css 분포 계산 (1 mg/kg/day, mg/L)
  css.dist_mgL <- calc_mc_css(
    parameters   = mc.params,
    model        = "pbtk",
    daily.dose   = 1,
    output.units = "mg/L"
  )

  # D) 95th 백분위수 & µM 환산
  css95_mgL <- quantile(css.dist_mgL, probs = 0.95)
  mw        <- get_physchem_param(param = "MW", chem.cas = cas, default.to.human = TRUE)  # g/mol
  css95_mgL * 1000 / mw                                          # mg/L → µM
}

# 5) 배치 실행 및 결과 출력
results <- sapply(chem.cas, calc_css95_uM)
print(data.frame(CAS = chem.cas, Css95_uM = round(results, 3)))

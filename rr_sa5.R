## 생존분석 (퇴원일 기준 time origin, diff = index_discharge → 첫 재입원)
##  - Cox: LOS는 원본(일수, VSHSP_DD_CNT) 그대로 사용
##  - KM : LOS 이분화(>=5일 vs <5일) - median값이기때문
##  - 기간설정: 퇴원일로부터 최대 3년까지만 추적(3년 초과 사건은 검열)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(survival)
  library(broom)
  library(readr)
  library(ggplot2)
  library(survminer)
  library(gridtext)
  library(grid)
})

index_rows <- final_cohort1 %>%
  inner_join(
    index_info %>% select(RN_INDI, first_admit, index_discharge, index_code),
    by = "RN_INDI"
  ) %>%
  filter(
    FORM_CD == "02",
    MDCARE_START_DATE == first_admit,
    SICK_SYM1 == index_code
  ) %>%
  distinct(RN_INDI, .keep_all = TRUE)

## 1) 동반질환(lookback 180일) 
hx_vars <- c("hx_hf","hx_ihd","hx_af","hx_htn","hx_dm","hx_dlp","hx_ckd","hx_copd","hx_stroke","hx_cancer")
lookback_comorb <- 180

history_prior_comorb <- final_cohort1 %>%
  inner_join(index_info %>% select(RN_INDI, first_admit), by = "RN_INDI") %>%
  filter(
    MDCARE_START_DATE < first_admit,
    MDCARE_START_DATE >= first_admit - days(lookback_comorb)
  )

comorb_180 <- history_prior_comorb %>%
  group_by(RN_INDI) %>%
  summarise(
    hx_hf     = any(str_detect(SICK_SYM1, "^I50")     | str_detect(SICK_SYM2, "^I50"),     na.rm=TRUE),
    hx_ihd    = any(str_detect(SICK_SYM1, "^I2[0-5]") | str_detect(SICK_SYM2, "^I2[0-5]"), na.rm=TRUE),
    hx_af     = any(str_detect(SICK_SYM1, "^I48")     | str_detect(SICK_SYM2, "^I48"),     na.rm=TRUE),
    hx_htn    = any(str_detect(SICK_SYM1, "^I1[0-5]") | str_detect(SICK_SYM2, "^I1[0-5]"), na.rm=TRUE),
    hx_dm     = any(str_detect(SICK_SYM1, "^E1[0-4]") | str_detect(SICK_SYM2, "^E1[0-4]"), na.rm=TRUE),
    hx_dlp    = any(str_detect(SICK_SYM1, "^E78")     | str_detect(SICK_SYM2, "^E78"),     na.rm=TRUE),
    hx_ckd    = any(str_detect(SICK_SYM1, "^N18")     | str_detect(SICK_SYM2, "^N18"),     na.rm=TRUE),
    hx_copd   = any(str_detect(SICK_SYM1, "^J4[0-7]") | str_detect(SICK_SYM2, "^J4[0-7]"), na.rm=TRUE),
    hx_stroke = any(str_detect(SICK_SYM1, "^I6[0-9]") | str_detect(SICK_SYM2, "^I6[0-9]"), na.rm=TRUE),
    hx_cancer = any(str_detect(SICK_SYM1, "^C\\d{2}") | str_detect(SICK_SYM2, "^C\\d{2}"), na.rm=TRUE),
    .groups = "drop"
  )

wd_tbl1 <- wd_tbl1 %>%
  left_join(comorb_180, by = "RN_INDI") %>%
  mutate(
    across(all_of(hx_vars),
           ~ factor(as.integer(replace_na(., FALSE)), levels = c(0,1), labels = c("No","Yes"))),
    age10  = if (!"age10" %in% names(.)) age_index/10 else age10,
    outpt5 = if (!"outpt5" %in% names(.)) (prior_03_n_6m %>% replace_na(0))/5 else outpt5,
    los_q  = if (!"los_q" %in% names(.)) cut(VSHSP_DD_CNT,
                                             breaks = quantile(VSHSP_DD_CNT, probs = c(0,.25,.5,.75,1), na.rm=TRUE),
                                             include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4")) else los_q
  )

## 2) 이벤트(첫 재입원): “퇴원일 이후” 동일상병 첫 입원
same_dx_after <- final_cohort1 %>%
  inner_join(index_rows %>% select(RN_INDI, index_discharge, index_code), by = "RN_INDI") %>%
  filter(
    FORM_CD == "02",
    MDCARE_START_DATE > index_discharge,
    SICK_SYM1 == index_code | SICK_SYM2 == index_code
  )

first_readmit_rows <- same_dx_after %>%
  group_by(RN_INDI) %>%
  slice_min(order_by = MDCARE_START_DATE, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(RN_INDI, first_readmit_date = MDCARE_START_DATE)

## 3) 생존분석 (퇴원일 → 이벤트/검열) + 퇴원일 이후 3년이내
study_end <- as.Date("2019-12-31")

## 사건 후보
surv_event <- index_rows %>%
  inner_join(first_readmit_rows, by = "RN_INDI") %>%
  transmute(
    RN_INDI,
    start = index_discharge,
    end   = first_readmit_date,
    event = 1L
  )

## 비사건(기본 검열: 연구 종료일)
surv_censor <- index_rows %>%
  anti_join(first_readmit_rows, by = "RN_INDI") %>%
  transmute(
    RN_INDI,
    start = index_discharge,
    end   = study_end,
    event = 0L
  )

surv_pre <- bind_rows(surv_event, surv_censor) %>%
  mutate(diff_days_raw = as.integer(end - start))

## 3년검열 적용
surv_cut <- surv_pre %>%
  mutate(
    cutoff   = pmin(start %m+% years(3), study_end),
    end_adj  = pmin(end, cutoff),
    event_adj= if_else(event == 1L & end <= cutoff, 1L, 0L),
    time     = pmax(0L, as.integer(end_adj - start)),
    event    = event_adj
  ) %>%
  select(RN_INDI, start, end = end_adj, time, event, cutoff, diff_days_raw)

## (C) 공변량 결합 + 파생변수 
surv_data <- surv_cut %>%
  left_join(
    wd_tbl1 %>% select(RN_INDI, age_index, age10, sex_f, VSHSP_DD_CNT, outpt5, index_dx_group_f, any_of(hx_vars)),
    by = "RN_INDI"
  ) %>%
  mutate(
    los   = as.numeric(VSHSP_DD_CNT),             
    LOS5  = factor(if_else(VSHSP_DD_CNT >= 5, ">=5d", "<5d")),
    los_q4= factor(if_else(
      VSHSP_DD_CNT >= quantile(VSHSP_DD_CNT, 0.75, na.rm = TRUE),
      "Q4 (High)", "Q1-3 (Low)"
    ))
  )

surv_data %>%
  summarise(
    n = n(),
    events = sum(event),
    median_fu_days = median(time, na.rm = TRUE),
    IQR_L = quantile(time, .25, na.rm = TRUE),
    IQR_U = quantile(time, .75, na.rm = TRUE)
  )

## 4) 단변량 Cox (동반질환)
uni_list <- lapply(hx_vars, function(v) {
  f <- as.formula(paste0("Surv(time, event) ~ ", v))
  fit <- survival::coxph(f, data = surv_data)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(var = v)
})
uni_hx <- dplyr::bind_rows(uni_list) %>%
  dplyr::select(var, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(across(c(estimate, conf.low, conf.high), ~round(., 2)),
                p.value = signif(p.value, 3)) %>%
  dplyr::arrange(p.value)
readr::write_csv(uni_hx, "cox_univariate_comorb.csv")

## 4) 단변량 Cox (동반질환 + 주요 변수들)
base_vars <- c("age10", "sex_f", "los")   
all_vars  <- c(base_vars, hx_vars)       

uni_list <- lapply(all_vars, function(v) {
  f <- as.formula(paste0("Surv(time, event) ~ ", v))
  fit <- survival::coxph(f, data = surv_data)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(var = v)
})

uni_all <- dplyr::bind_rows(uni_list) %>%
  dplyr::select(var, estimate, conf.low, conf.high, p.value) %>%
  dplyr::mutate(across(c(estimate, conf.low, conf.high), ~round(., 2)),
                p.value = signif(p.value, 3)) %>%
  dplyr::arrange(p.value)

readr::write_csv(uni_all, "cox_univariate_all.csv")
## 5) 다변량 Cox 
fit_cox_main <- survival::coxph(
  Surv(time, event) ~ age10 + sex_f + los + hx_htn,
  data = surv_data
)
summary(fit_cox_main)
cox.zph(fit_cox_main)

hr_main <- broom::tidy(fit_cox_main, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::mutate(across(c(estimate, conf.low, conf.high), ~round(., 2)),
                p.value = signif(p.value, 3))
readr::write_csv(hr_main, "cox_main_LOS_results.csv")

##  Cox 메인 모형 성능 요약 
s <- summary(fit_cox_main)
k   <- length(coef(fit_cox_main))
aic <- -2 * s$loglik[2] + 2 * k
bic <- -2 * s$loglik[2] + log(s$n) * k

LR_stat   <- unname(s$logtest[1]);  LR_df <- unname(s$logtest[2]);  LR_p <- unname(s$logtest[3])
Wald_stat <- unname(s$waldtest[1]); Wald_df <- unname(s$waldtest[2]); Wald_p <- unname(s$waldtest[3])
Score_stat<- unname(s$sctest[1]);  Score_df<- unname(s$sctest[2]);  Score_p<- unname(s$sctest[3])

perf_tbl <- tibble::tibble(
  model            = "cox_main_LOS_raw_3y_censor",
  n                = s$n,
  events           = s$nevent,
  concordance      = round(s$concordance[1], 3),
  concordance_se   = round(s$concordance[2], 3),
  loglik_start     = round(s$loglik[1], 3),
  loglik_final     = round(s$loglik[2], 3),
  LR_chisq         = round(LR_stat, 3),
  LR_df            = LR_df,
  LR_p             = signif(LR_p, 3),
  Wald_chisq       = round(Wald_stat, 3),
  Wald_df          = Wald_df,
  Wald_p           = signif(Wald_p, 3),
  Score_chisq      = round(Score_stat, 3),
  Score_df         = Score_df,
  Score_p          = signif(Score_p, 3),
  AIC              = round(aic, 3),
  BIC              = round(bic, 3)
)
readr::write_csv(perf_tbl, "cox_main_performance.csv")

cz <- cox.zph(fit_cox_main)
ph_tbl <- as.data.frame(cz$table)
ph_tbl$term <- rownames(ph_tbl)
ph_tbl <- dplyr::select(ph_tbl, term, chisq, p)
readr::write_csv(ph_tbl, "cox_main_PH_test.csv")

## Kaplan–Meier & 누적위험 (30일, 90일, 1년, 2년, 3년)
fit_all <- survival::survfit(Surv(time, event) ~ 1, data = surv_data)

landmarks <- c(30, 90, 365, 730, 1095)
km_pts  <- summary(fit_all, times = landmarks)
km_tab  <- data.frame(
  time = km_pts$time,
  surv = round(km_pts$surv, 3),
  risk = round(1 - km_pts$surv, 3)
)
readr::write_csv(km_tab, "km_overall_landmarks.csv")


df <- km_tab %>%
  transmute(time = factor(time,
                          levels = landmarks,
                          labels = c("30일","90일","1년","2년","3년")),
            risk_pct = risk * 100)


p_bar <- ggplot(df, aes(time, risk_pct)) +
  geom_col(width = 0.55, fill = "grey40") +
  geom_text(aes(label = sprintf("%.1f%%", risk_pct)),
            vjust = -0.6, size = 4) +
  scale_y_continuous(
    name = "누적 재입원 위험(%)",
    expand = expansion(mult = c(0, 0.08)),
    limits = c(0, max(df$risk_pct) * 1.18)
  ) +
  labs(x = NULL, title = "누적 재입원 위험률 (30일·90일·1년·2년·3년)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("plots/overall_risk_landmarks_bar.png",
       p_bar, width = 7, height = 3.8, dpi = 300)



## LOS 이분화: 5일 컷(los median이 5일)
fit_km_LOS5 <- survival::survfit(Surv(time, event) ~ LOS5, data = surv_data)
diff_LOS5   <- survival::survdiff(Surv(time, event) ~ LOS5, data = surv_data)
p_LOS5      <- survminer::ggsurvplot(fit_km_LOS5, data = surv_data, pval = TRUE, risk.table = TRUE)
ggplot2::ggsave("km_by_LOS5.png", p_LOS5$plot, width = 6, height = 4, dpi = 300)
if (inherits(p_LOS5$table, "ggplot")) {
  ggplot2::ggsave("km_by_LOS5_risktable.png", p_LOS5$table, width = 6, height = 2.5, dpi = 300)
}


png("plots/km_by_LOS5_full.png", width = 7, height = 5.5, units = "in", res = 300)
print(p_LOS5)
dev.off()


sink("sessionInfo.txt"); print(sessionInfo()); sink()


dir.create("plots", showWarnings = FALSE)
safe_ggsave <- function(file, plt, width = 6, height = 4, dpi = 300) {
  tryCatch(
    ggplot2::ggsave(file, plt, width = width, height = height, dpi = dpi),
    error = function(e) message("[ggsave skipped] ", file, " :: ", e$message)
  )
}



##시각화
#  Overall KM (S(t)) 
fit_all <- survival::survfit(Surv(time, event) ~ 1, data = surv_data)
p_overall <- survminer::ggsurvplot(
  fit_all, data = surv_data,
  conf.int = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(), title = "Overall KM (S(t))"
)
safe_ggsave("plots/km_overall.png", p_overall$plot)
if (inherits(p_overall$table, "ggplot"))
  safe_ggsave("plots/km_overall_risktable.png", p_overall$table, width = 6, height = 2.5)

#  Overall 누적발생(=1-S(t)) 
p_overall_event <- survminer::ggsurvplot(
  fit_all, data = surv_data,
  fun = "event", conf.int = FALSE, risk.table = TRUE,
  ggtheme = theme_minimal(), title = "Overall Cumulative Incidence (1 - S(t))"
)
safe_ggsave("plots/km_overall_cumevent.png", p_overall_event$plot)
if (inherits(p_overall_event$table, "ggplot"))
  safe_ggsave("plots/km_overall_cumevent_risktable.png", p_overall_event$table, width = 6, height = 2.5)

#  Overall 누적위험(H(t)) 
p_overall_cumhaz <- survminer::ggsurvplot(
  fit_all, data = surv_data,
  fun = "cumhaz", conf.int = FALSE, risk.table = FALSE,
  ggtheme = theme_minimal(), title = "Overall Cumulative Hazard (H(t))"
)
safe_ggsave("plots/km_overall_cumhaz.png", p_overall_cumhaz$plot)



#  KM: LOS 5일 컷 
fit_km_LOS5 <- survival::survfit(Surv(time, event) ~ LOS5, data = surv_data)
p_LOS5 <- survminer::ggsurvplot(
  fit_km_LOS5, data = surv_data,
  pval = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(), title = "KM by LOS (>=5d vs <5d)",
  legend.title = "LOS5", legend.labs = levels(droplevels(surv_data$LOS5))
)
safe_ggsave("plots/km_by_LOS5.png", p_LOS5$plot)
if (inherits(p_LOS5$table, "ggplot"))
  safe_ggsave("plots/km_by_LOS5_risktable.png", p_LOS5$table, width = 6, height = 2.5)





## 다변량 Cox 결과 그림 저장

dir.create("plots", showWarnings = FALSE)

safe_ggsave <- function(file, plt, width = 7, height = 5, dpi = 300) {
  tryCatch(ggplot2::ggsave(file, plt, width = width, height = height, dpi = dpi),
           error = function(e) message("[ggsave skipped] ", file, " :: ", e$message))
}
safe_png_print <- function(file, plot_obj, width = 7, height = 5, res = 300, units = "in") {
  tryCatch({
    grDevices::png(filename = file, width = width, height = height, units = units, res = res)
    print(plot_obj); grDevices::dev.off()
  }, error = function(e) message("[png print skipped] ", file, " :: ", e$message))
}

df_for_forest <- stats::model.frame(fit_cox_main)

p_forest_auto <- survminer::ggforest(
  fit_cox_main,
  data = df_for_forest,                  
  main = "Cox: age10 + sex + LOS + HTN",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.8
)

ggplot2::ggsave("plots/cox_forest_auto.png", p_forest_auto, width = 7, height = 5, dpi = 300)

#추가 km시각화
landmarks <- c(30, 90, 365, 730, 1094)
fit_all <- survfit(Surv(time, event) ~ 1, data = surv_data)

p_landmark <- ggsurvplot(
  fit_all, data = surv_data,
  conf.int = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(),
  break.time.by = 365,
  title = "Overall Kaplan–Meier (재입원-free survival)"
)

km_pts <- summary(fit_all, times = landmarks, extend = TRUE)
landmark_df <- data.frame(
  time = km_pts$time,
  surv = km_pts$surv,
  label_time = factor(km_pts$time,
                      levels = landmarks,
                      labels = c("30일", "90일", "1년", "2년", "3년"))
)



p_landmark$plot <- p_landmark$plot +
  geom_point(data = landmark_df,
             aes(x = time, y = surv),
             color = "red", size = 2) +
  geom_text(data = landmark_df,
            aes(x = time, y = surv,
                label = paste0(label_time, "\n", sprintf("%.1f%%", (1-surv)*100))),
            vjust = -0.6, size = 3)


ggsave("plots/km_overall_landmarks.png",
       p_landmark$plot, width = 8, height =6, dpi = 300)



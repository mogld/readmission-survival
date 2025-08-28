library(tidyverse)
library(tableone)
# --- 첫 재입원일 추출 ---
first_readmit_date <- same_dx_events_all %>%
  group_by(RN_INDI) %>%
  summarise(
    first_readmit_dt = min(MDCARE_START_DATE, na.rm = TRUE),
    .groups = "drop"
  )

wd <- wd %>%
  left_join(first_readmit_date %>%
              rename(first_readmit_dt_raw = first_readmit_dt),
            by = "RN_INDI") %>%
  mutate(
    first_readmit_dt = if_else(readmit_flag == 1L,
                               first_readmit_dt_raw,
                               as.Date(NA_character_))
  ) %>%
  select(-first_readmit_dt_raw)
# --- 분석용 파생변수: wd -> wd_tbl1 ---
wd_tbl1 <- wd %>%
  mutate(
    group  = factor(if_else(readmit_flag == 1L, "readmit", "nonreadmit"),
                    levels = c("nonreadmit","readmit")),
    sex_f  = factor(SEX, levels = c(1,2), labels = c("M","F")),
    
    index_dx_group_f = case_when(
      str_detect(index_code, "^I(2[0-9]|3[0-9]|4[0-9]|5[0-1])") ~ "심장질환",  # I20–I51
      str_detect(index_code, "^I6[0-9]")                         ~ "뇌혈관질환",# I60–I69
      TRUE ~ "기타"
    ),
    index_dx_group_f = factor(index_dx_group_f, 
                              levels = c("심장질환","뇌혈관질환","기타")),
    
    path_type_f     = factor(HSPTZ_PATH_TYPE),
    er_flag         = factor(as.integer(HSPTZ_PATH_TYPE %in% c(11,21,31)),
                             levels = c(0,1), labels = c("No","Yes")),
    

    mcare_result_f  = factor(MCARE_RSLT_TYPE),
    
    index_year_f    = factor(MDCARE_START_YEAR),
    index_cat       = factor(substr(index_code, 1, 3)),
    death_flag      = factor(!is.na(DTH_YYYYMM), levels = c(FALSE, TRUE),
                             labels = c("No","Yes")),
    
    diff_days = case_when(
      group == "readmit" ~ as.numeric(first_readmit_dt - DISCHARGE_DT),
      TRUE ~ NA_real_
    )
  )

# Table 1 변수 
# 연속형
vars_cont <- c(
  "age_index", "VSHSP_DD_CNT", "diff_days",
  "prior_03_n_6m", "prior_emer_n_6m", "prior_02_n_6m"
)

# 범주형
vars_cat <- c(
  "sex_f","er_flag","path_type_f","mcare_result_f",
  "index_year_f","index_cat","death_flag","index_dx_group_f"
)

vars_all       <- c(vars_cont, vars_cat)
factorVars     <- vars_cat

nonnormalVars  <- c("VSHSP_DD_CNT","diff_days",
                    "prior_03_n_6m","prior_emer_n_6m","prior_02_n_6m")

# --- Table 1 생성  
tab1 <- CreateTableOne(
  vars       = vars_all,
  strata     = "group",     
  data       = wd_tbl1,
  factorVars = factorVars
)

print(tab1,
      showAllLevels = TRUE,   
      addOverall    = TRUE,   
      smd           = TRUE,  
      test          = TRUE,   
      nonnormal     = nonnormalVars,
      exact         = c("sex_f","er_flag","death_flag","index_dx_group_f"),
      quote         = FALSE,
      noSpaces      = TRUE,
      digits        = 3,      
      digits.nonnum = 1        
)


# --- CSV 저장
tab1_mat <- print(tab1,
                  showAllLevels = TRUE,
                  addOverall    = TRUE,   
                  smd           = TRUE,
                  test          = TRUE,
                  nonnormal     = nonnormalVars,
                  exact         = c("sex_f","er_flag","death_flag","index_dx_group_f"),
                  printToggle   = FALSE,
                  noSpaces      = TRUE,
                  digits        = 3,
                  digits.nonnum = 1,
                  includeNA     = TRUE)

write.csv(as.data.frame(tab1_mat), "Table_wd_final.csv", row.names = TRUE)

tab1_overall <- CreateTableOne(
  vars       = vars_all,
  data       = wd_tbl1,         # strata 없이 전체
  factorVars = factorVars
)

overall_tbl <- print(tab1_overall,
                     showAllLevels = TRUE,
                     addOverall    = FALSE, # strata 없으므로 자체가 Overall
                     smd           = FALSE,
                     test          = FALSE,
                     nonnormal     = nonnormalVars,
                     printToggle   = FALSE,
                     digits        = 3,
                     digits.nonnum = 1)

overall_df <- as.data.frame(overall_tbl)
write.csv(overall_df, "Table_overall.csv", row.names = TRUE)

# --- 표본수/결측률 확인
cat("N(wd_tbl1):", nrow(wd_tbl1), 
    " | IDs:", dplyr::n_distinct(wd_tbl1$RN_INDI), "\n")

wd_tbl1 %>% count(group)

wd_tbl1 %>%
  summarise(across(c(age_index, VSHSP_DD_CNT, diff_days),
                   ~ mean(is.na(.))*100, .names="{.col}_NA%"))

# --- summary 전체 기술통계
summary(wd_tbl1)


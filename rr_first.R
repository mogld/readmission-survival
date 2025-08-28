# 진료+건강검진+출생사망+자격데이터 합치기
library(tidyverse)
library(readxl)

m20_data <- read_csv("NSC2_M20.csv")
m20_1619_data <- read_csv("NSC2_M20_1619.csv")

m20 <- bind_rows(m20_data, m20_1619_data)
glimpse(m20)

selected_m20 <- m20 %>%
  select(
    RN_INDI,
    MDCARE_STRT_DT,
    FORM_CD,
    SICK_SYM1,
    SICK_SYM2,
    HSPTZ_PATH_TYPE,
    VSHSP_DD_CNT,
    MCARE_RSLT_TYPE,
    FST_HSPTZ_DT
  )

BD <- read_csv("NSC2_BND.csv")
selected_BD <- BD %>%
  select(
    RN_INDI,
    BTH_YYYY,
    DTH_YYYYMM
  ) %>%
  distinct(RN_INDI, .keep_all = TRUE)

gender1 <- read_csv("NSC2_BNC_V2_1.csv")
gender2 <- read_csv("NSC2_BNC_1619.csv")
gender<- bind_rows(gender1, gender2)
selected_gender <- gender %>%
  select(
    RN_INDI,
    SEX
  ) %>%
  distinct(RN_INDI, .keep_all = TRUE)

all_data <- selected_m20 %>%
  inner_join(selected_BD, by = "RN_INDI") %>%
  inner_join(selected_gender, by = "RN_INDI") 
glimpse(all_data)

# 대상자필터링
# 1. 주상병-심장질환이나 뇌혈관질환때문에 입원한적이 있는 환자
patient_list <- all_data %>%
  filter(
    str_detect(SICK_SYM1, "^I(2[0-9]|3[0-9]|4[0-9]|5[0-1])") |
      str_detect(SICK_SYM1, "^I6[0-9]") 
  ) %>%
  filter(FORM_CD %in% c("02")) %>%
  select(RN_INDI) %>%
  distinct()

# 2. 1에 속하는 환자의 모든 기록  -21126, 62명
final_cohort <- all_data %>%
  semi_join(patient_list, by = "RN_INDI") %>%
  mutate(
    MDCARE_START_DATE = ymd(MDCARE_STRT_DT),
    MDCARE_START_YEAR = year(MDCARE_START_DATE)
  ) %>%
  filter(MDCARE_START_YEAR >= 2006 & MDCARE_START_YEAR <= 2019)

write_csv(final_cohort,"final_cohort.csv")

obs_start <- as.Date("2006-01-01")
icd_pat   <- "^I(2[0-9]|3[0-9]|4[0-9]|5[0-1]|6[0-9])"

# 3. index date 설정: 심장/뇌혈관질환으로 '입원'한 최초 날짜
index_date <- final_cohort %>% 
  filter(FORM_CD == "02", 
         MDCARE_START_DATE >= obs_start, 
         str_detect(SICK_SYM1, "^I(2[0-9]|3[0-9]|4[0-9]|5[0-1])") | str_detect(SICK_SYM1, "^I6[0-9]")) %>% 
        group_by(RN_INDI) %>% slice_min(MDCARE_START_DATE, with_ties = FALSE) %>% 
        ungroup() %>% 
        transmute( RN_INDI, first_admit = MDCARE_START_DATE, index_code = SICK_SYM1, BTH_YYYY )

# 4. 워시아웃 대상자: 2006년 이전에 주/부상병 중 하나라도 index_date 이력이 있던 환자
washout_ids <- final_cohort %>%
  inner_join(index_date, by = "RN_INDI") %>%
  filter(MDCARE_START_DATE < obs_start) %>%
  filter(SICK_SYM1 == index_code | SICK_SYM2 == index_code) %>%
  pull(RN_INDI) %>%
  unique()

final_cohort <- final_cohort %>%
  filter(MDCARE_START_DATE >= obs_start,
         !RN_INDI %in% washout_ids)

index_date <- index_date %>%
  filter(!RN_INDI %in% washout_ids) %>%
  select(RN_INDI, first_admit, index_SICK_SYM1 = index_code, BTH_YYYY)

# 퇴원일 계산
final_cohort <- final_cohort %>%
  mutate(
    DISCHARGE_DT = if_else(
      FORM_CD %in% c("02"),
      MDCARE_START_DATE + days(VSHSP_DD_CNT),  
      as.Date(NA)  
    )
  )

# 4. 입원 중 사망 여부 확인 ('월'기준)
death <- final_cohort %>%
  filter(FORM_CD %in% c("02","07","12")) %>%
  filter(!is.na(DTH_YYYYMM)) %>%
  mutate(death_month = ym(DTH_YYYYMM)) %>%
  mutate(
    admit_month = floor_date(MDCARE_START_DATE, "month"),
    discharge_month = floor_date(DISCHARGE_DT, "month")
  ) %>%
  filter(death_month >= admit_month & death_month <= discharge_month) %>%
  pull(RN_INDI) %>%
  unique()

# 5. 퇴원 기록 불분명한 환자 제외
unclear_discharge <- final_cohort %>%
  filter(FORM_CD == "02") %>%  
  filter(is.na(DISCHARGE_DT) | is.na(VSHSP_DD_CNT) | DISCHARGE_DT < MDCARE_START_DATE) %>%
  pull(RN_INDI) %>%
  unique()

excluded_ids <- union(death, unclear_discharge)
final_cohort <- final_cohort %>%
  filter(!RN_INDI %in% excluded_ids)

# 필요한 변수생성: 첫입원일 기준 나이계산
lookback_days <- 180  

# 수정된 age_at_index 계산 코드
age_at_index <- index_date %>%
  mutate(
    first_admit = as.Date(first_admit),
    age_index = year(first_admit) - as.integer(BTH_YYYY)
  ) %>%
  select(RN_INDI, age_index)


# age na이거나 0세 이하이면 연구대상자에서 제외 -2명제외, 20386,60명
invalid_age_ids <- age_at_index %>%
  filter(is.na(age_index) | age_index < 0) %>%
  pull(RN_INDI) %>%
  unique()

final_cohort <- final_cohort %>%
  filter(!RN_INDI %in% invalid_age_ids)
index_date <- index_date %>%
  filter(!RN_INDI %in% invalid_age_ids)






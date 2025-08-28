
# 0) index 입원
index_info <- final_cohort %>%
  inner_join(index_date, by = "RN_INDI") %>%
  filter(FORM_CD == "02", MDCARE_START_DATE == first_admit) %>%
  arrange(RN_INDI, DISCHARGE_DT) %>%
  group_by(RN_INDI) %>% slice(1) %>% ungroup() %>%
  transmute(
    RN_INDI,
    first_admit,
    index_discharge = DISCHARGE_DT,
    index_code      = index_SICK_SYM1
  )


# 3) index 이후 '동일 상병코드' 입원 횟수
same_dx_inclusive <- final_cohort %>%
  inner_join(index_info, by = "RN_INDI") %>%
  filter(
    FORM_CD == "02",
    MDCARE_START_DATE >= first_admit,
    SICK_SYM1 == index_code | SICK_SYM2 == index_code
  ) %>%
  group_by(RN_INDI) %>%
  summarise(
    same_dx_index = n_distinct(as.Date(MDCARE_START_DATE)),
    .groups = "drop"
  )

# 1) 환자 단위 그룹 라벨링 (1=재입원, 0=비재입원)
readmit01 <- same_dx_inclusive %>%
  transmute(
    RN_INDI,
    same_dx_count = as.integer(same_dx_index >= 2L)  
  )
index_info <- index_info %>% left_join(readmit01, by = "RN_INDI")

readmit01 <- index_info %>%
  select(RN_INDI) %>%
  left_join(same_dx_inclusive, by = "RN_INDI") %>%
  mutate(
    n_same_dx    = coalesce(same_dx_index, 1L),  
    readmit_flag = as.integer(n_same_dx >= 2L)               
  ) %>%
  select(RN_INDI, readmit_flag, n_same_dx)

# 1) 그룹 분리 
patients_readmit    <- readmit01 %>% filter(readmit_flag == 1L) %>% select(RN_INDI)
patients_nonreadmit <- readmit01 %>% filter(readmit_flag == 0L) %>% select(RN_INDI)


##최최종 연구대상자 
final_cohort1 <- final_cohort %>%
  semi_join(index_date %>% select(RN_INDI) %>% distinct(), by = "RN_INDI") %>%
  arrange(RN_INDI, MDCARE_START_DATE)
final_cohort1 <- final_cohort1 %>%
  left_join(age_at_index %>% select(RN_INDI, age_index), by = "RN_INDI")

# 2) 사건 단위 테이블 (id당 최대 2행: index, 첫 재입원) 
# 2-1) index 입원 행 (모든 대상자)
index_rows <- final_cohort1 %>%
  inner_join(index_info %>% select(RN_INDI, first_admit, index_discharge, index_code),
             by = "RN_INDI") %>%
  filter(
    FORM_CD == "02",
    MDCARE_START_DATE == first_admit,
    SICK_SYM1 == index_code            
  ) %>%
  distinct(RN_INDI, .keep_all = TRUE)  

#최종 그룹별 전체테이블 
final_readmit <- final_cohort1 %>%
  semi_join(patients_readmit, by = "RN_INDI") %>%
  arrange(RN_INDI, MDCARE_START_DATE)

final_nonreadmit <- final_cohort1 %>%
  semi_join(patients_nonreadmit, by = "RN_INDI") %>%
  arrange(RN_INDI, MDCARE_START_DATE)


# 재입원군: index + 첫 재입원(최대 2행)
same_dx_events_all <- final_cohort1 %>%
  inner_join(index_info %>% select(RN_INDI, first_admit, index_code), by = "RN_INDI") %>%
  filter(
    FORM_CD == "02",
    MDCARE_START_DATE > first_admit
  ) %>%
  group_by(RN_INDI, MDCARE_START_DATE, DISCHARGE_DT) %>%
  summarise(
    same_dx = any(SICK_SYM1 == index_code | SICK_SYM2 == index_code, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  ungroup() %>%
  filter(same_dx) %>%
  arrange(RN_INDI, MDCARE_START_DATE, DISCHARGE_DT)

last_readmit_events <- same_dx_events_all %>%
  group_by(RN_INDI) %>%
  slice_max(order_by = MDCARE_START_DATE, with_ties = FALSE) %>%
  ungroup()
last_readmit_rows <- final_cohort1 %>%
  inner_join(last_readmit_events,
             by = c("RN_INDI", "MDCARE_START_DATE", "DISCHARGE_DT")) %>%
  filter(FORM_CD == "02") %>%
  group_by(RN_INDI) %>% slice(1) %>% ungroup()


# 재입원군 pair(인덱스 + 최종 재입원)
readmit_pairs <- bind_rows(
  index_rows %>% semi_join(patients_readmit, by = "RN_INDI"),
  last_readmit_rows %>% semi_join(patients_readmit, by = "RN_INDI")
) %>%
  arrange(RN_INDI, MDCARE_START_DATE)

#비재입원군
nonreadmit_rows <- index_rows %>%
  semi_join(patients_nonreadmit, by = "RN_INDI")
#재입원군 
readmit_index_only <- index_rows %>%
  semi_join(patients_readmit, by = "RN_INDI") %>%
  arrange(RN_INDI, MDCARE_START_DATE) %>%
  distinct(RN_INDI, .keep_all = TRUE)

# 비재입원군: index만
nonreadmit_index_only <- index_rows %>%
  semi_join(patients_nonreadmit, by = "RN_INDI") %>%
  arrange(RN_INDI, MDCARE_START_DATE) %>%
  distinct(RN_INDI, .keep_all = TRUE)

# 두 그룹 최종 데이터테이블 
wd <- bind_rows(readmit_index_only, nonreadmit_index_only) %>%
  distinct(RN_INDI, .keep_all = TRUE) %>%
  left_join(readmit01 %>% select(RN_INDI, readmit_flag), by = "RN_INDI")

# 3) 재입원일-입원일 간격 계산 
intervals_readmit <- readmit_pairs %>%
  arrange(RN_INDI, MDCARE_START_DATE) %>%
  group_by(RN_INDI) %>%
  mutate(
    diff_days = lead(MDCARE_START_DATE) - MDCARE_START_DATE
  ) %>%
  ungroup() 

end_date <- ymd("2019-12-31")
intervals_nonreadmit <- nonreadmit_rows %>%
  mutate(
  diff_days = end_date - MDCARE_START_DATE
  )

diff_readmit <- intervals_readmit %>%
  filter(!is.na(diff_days)) %>%            # 인덱스행만(재입원행은 NA)
  select(RN_INDI, diff_days) %>%
  distinct(RN_INDI, .keep_all = TRUE)

diff_nonreadmit <- intervals_nonreadmit %>%
  select(RN_INDI, diff_days) %>%
  distinct(RN_INDI, .keep_all = TRUE)

# 1) time_df 생성 
time_df <- bind_rows(
  diff_readmit,
  diff_nonreadmit   
) %>%
  distinct(RN_INDI, .keep_all = TRUE) %>%
  mutate(diff_days = as.integer(diff_days))

wd <- wd %>%
  left_join(time_df, by = "RN_INDI")

# 1) 과거 180일 (lookback) 진료 이력만 추출
history_prior_6m <- final_cohort1 %>%
  inner_join(index_info, by = "RN_INDI") %>%
  filter(
    MDCARE_START_DATE < first_admit,
    MDCARE_START_DATE >= first_admit - days(lookback_days)
  )

# 2) 과거 6개월: 외래/응급/입원 횟수
prior_6m <- history_prior_6m %>%
  left_join(index_date %>% select(RN_INDI, index_SICK_SYM1), by = "RN_INDI") %>%
  group_by(RN_INDI) %>%
  summarise(
    prior_03_n_6m  = sum(FORM_CD != "02", na.rm = TRUE),
    prior_emer_n_6m = sum(HSPTZ_PATH_TYPE %in% c(11, 21, 31), na.rm = TRUE),
    prior_02_n_6m = sum(
      FORM_CD == "02" &
        !is.na(SICK_SYM1) & !is.na(index_SICK_SYM1) &
        SICK_SYM1 != index_SICK_SYM1,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

wd <- wd %>%
  left_join(prior_6m, by = "RN_INDI")








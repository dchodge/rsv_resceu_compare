

calc_impact <- function(base, programme, filename) {
    age_groups <- c("0-3mo", "4-6mo", "7–11mo", "12–23mo", "24-59mo", "0-11mo", "0-59mo")

    # Cases averted
    asymp_cases_avert <- (base$cases$asymp$under5[,-1] %>% rowSums) - (programme$cases$asymp$under5[,-1] %>% rowSums)
    symp_nhc_cases_avert <- (base$cases$symp_nhc$under5[,-1] %>% rowSums) - (programme$cases$symp_nhc$under5[,-1] %>% rowSums)
    pri_cases_avert <- (base$cases$pri$under5[,-1] %>% rowSums) - (programme$cases$pri$under5[,-1] %>% rowSums)
    sec_cases_avert <- (base$cases$sec$under5[,-1] %>% rowSums) - (programme$cases$sec$under5[,-1] %>% rowSums)
    hosp_cases_avert <- (base$cases$hosp$under5[,-1] %>% rowSums) - (programme$cases$hosp$under5[,-1] %>% rowSums)
    icu_cases_avert <- (base$cases$icu$under5[,-1] %>% rowSums) - (programme$cases$icu$under5[,-1] %>% rowSums)
    death_cases_avert <- (base$cases$death$under5[,-1] %>% rowSums) - (programme$cases$death$under5[,-1] %>% rowSums)

    # Incremental QALY gained (undiscounted)
    qaly_undist <- base$QALY$under5$undiscounted_total - programme$QALY$under5$undiscounted_total

    # Incremental direct costs (undiscsounted, without intervention)
    direct_cost_undist <- -(programme$cost$under5$undiscounted_direct - base$cost$under5$undiscounted_direct)

    # Incremental indirect costs (undiscsounted, without intervention)
    indirect_cost_undist <- -(programme$cost$under5$undiscounted_indirect  - base$cost$under5$undiscounted_indirect)

    # Incremental QALY gained (discounted)
    qaly_dist <- base$QALY$under5$discounted_total - programme$QALY$under5$discounted_total

    # Incremental direct costs (discsounted, without intervention)
    direct_cost_dist <- -(programme$cost$under5$discounted_direct - base$cost$under5$discounted_direct)

    # Intervention costs
    inter_cost <- programme$cost$under5$discounted_inter - base$cost$under5$discounted_inter

    # Total costs
    total_cost <- inter_cost - direct_cost_dist

    # acer (payers perspective)
    acer_nhs <- c(rep(NA, 5), (total_cost / qaly_dist)[6:7])

    # indirect costs
    indirect_cost <- - (programme$cost$under5$discounted_indirect  - base$cost$under5$discounted_indirect)

    total_cost_societal <- total_cost - indirect_cost
    # acer (societal perspective)
    acer_societal <- c(rep(NA, 5), ((total_cost_societal) / qaly_dist)[6:7])

    output_df <- bind_cols(age_groups,
        asymp_cases_avert, symp_nhc_cases_avert,
        pri_cases_avert, sec_cases_avert, hosp_cases_avert, icu_cases_avert, death_cases_avert,
        qaly_undist, direct_cost_undist, indirect_cost_undist,
        qaly_dist, direct_cost_dist, inter_cost, total_cost, 
        acer_nhs, indirect_cost, total_cost_societal, acer_societal) %>% as.data.frame
    colnames(output_df)  <- c(
        "age_group", "asymp_averted", "symp_nhc_averted", "pri_averted", "sec_averted", "hosp_averted", "icu_averted", "death_averted", 
        "incr_QALY_undis", "incr_direct_cost_undis", "indirect_cost_undist", 
        "incr_QALY_dis", "incr_direct_cost_dis", "inter_cost", "total_costd_dist", 
        "acer_payers", "indirect_cost", "total_cost_societal", "acer_societal"
        )

    write.csv(output_df, file = here::here("outputs", "impact", paste0(filename, ".csv") ))

}

save_excel <- function(none_final, filename) {
    asymp_cases <- none_final$cases$asymp$under5
    symp_nhc_cases <- none_final$cases$symp_nhc$under5

  pri_cases <- none_final$cases$pri$under5

  pri_cases <- none_final$cases$pri$under5
  sec_cases <- none_final$cases$sec$under5
  hosp_cases <- none_final$cases$hosp$under5
  icu_cases <- none_final$cases$icu$under5
  death_cases <- none_final$cases$death$under5
  costs <- none_final$cost$under5
  qaly <- none_final$QALY$under5

  write.csv(asymp_cases, here::here("outputs", filename, "asymp_cases.csv"))
  write.csv(symp_nhc_cases, here::here("outputs", filename, "symp_nhc_cases.csv"))
  write.csv(pri_cases, here::here("outputs", filename, "pri_cases.csv"))
  write.csv(sec_cases, here::here("outputs", filename, "sec_cases.csv"))
  write.csv(hosp_cases, here::here("outputs", filename, "hosp_cases.csv"))
  write.csv(icu_cases, here::here("outputs", filename, "icu_cases.csv"))
  write.csv(death_cases, here::here("outputs", filename, "death_cases.csv"))
  write.csv(costs, here::here("outputs", filename, "costs.csv"))
  write.csv(qaly, here::here("outputs", filename, "qaly.csv"))

}

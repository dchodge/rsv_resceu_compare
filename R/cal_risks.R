cal_risks <- function(outputBase, resceudata) {
    all_cases <- outputBase$outcomes_month_age %>% filter(outcome == "all_cases")
    all_cases_mean <- all_cases %>% group_by(time, age_group) %>% summarise(inci_all = mean(incidence, na.rm = TRUE))
    all_cases_mean <- all_cases_mean %>% mutate(time  = time - (min(all_cases_mean$time) - 1) )

    ######################
    ## Hospital Inpatients ##
    ######################
    # see R/gen_data_resceu.R for cleaning of this data

    # Add to dataframe
    hosp_cases <- resceudata$observationalData %>%
        pivot_longer(everything(), values_to = "inci_h", names_to = "cols") %>%
        mutate(age_group = rep(1:25, 12), time = unlist(1:12 %>% map(~rep(.x, 25)))) %>%
        dplyr::select(inci_h, age_group, time)

    all_probs_df <- all_cases_mean %>% left_join(hosp_cases) %>% mutate(prop_h = inci_h / inci_all)

    ######################
    ## ICU visits ##
    ######################
    # % of hosp admission resulting in ICU
    prob_ICU_hosp <- c(1.85, 0.85, 0.81, 0.76, 0.39, 0.87, 0.27, rep(0.2, 5), rep(.7, 4), rep(0, 9)) / 100
    
    # Add to dataframe
    all_probs_df_trim <- all_probs_df %>%
        mutate(prob_ICU = prob_ICU_hosp) %>%
        ungroup %>%
        as.data.frame %>%
        mutate(prob_ICU_cases = prop_h * prob_ICU) %>%
        dplyr::select(time, age_group, prop_h, prob_ICU_cases)

    ######################
    ## Hospital inpatients ##
    ######################
    ageGroupBoundary <- resceudata$ageGroupBoundary
    population <- resceudata$population

    # Import from excel data 
    raw_hospout_age <- read.csv(file = here("inst", "extdata", "hosp_in_raw.csv"))[, 1]
    raw_hospout <- read.csv(file = here("inst", "extdata", "hosp_in_raw.csv"))[, 2:13]

    empty_older <- as.data.frame(matrix(0, 14, 12))
    names(empty_older) <- names(raw_hospout)

    rates_raw_hospout <- bind_rows(raw_hospout[1:12, ],
        raw_hospout[13, ],
        raw_hospout[14, ],
        raw_hospout[14, ],
        raw_hospout[14, ],
        empty_older)

    rates_real_hospout <- bind_cols(rates_raw_hospout[, 7:12], rates_raw_hospout[, 1:6])

    rates_real_hospout_list <- list()
    for (i in 1:25) {
        if (i < 25) {
            rates_real_hospout_list[[i]] <- rates_real_hospout[i, ] / 1000 * 100000 * (ageGroupBoundary[i + 1] - ageGroupBoundary[i]) / 12
        }
        else {
            rates_real_hospout_list[[i]] <- rates_real_hospout[i, ] / 1000 * (population - 100000*ageGroupBoundary[i-1]) / 12
        }
    }
    incidence_hospout <- rates_real_hospout_list %>% bind_rows %>% t %>% as.data.frame


    hosp_out_cases <- incidence_hospout %>%
        pivot_longer(everything(), values_to = "inci_h_out", names_to = "cols") %>%
        mutate(age_group = rep(1:25, 12), time = unlist(1:12 %>% map(~rep(.x, 25)))) %>%
        dplyr::select(age_group, time, inci_h_out)

    all_probs_df <- all_cases_mean %>% left_join(hosp_out_cases) %>% mutate(prop_h_out = inci_h_out / inci_all)

    all_probs_df_trim <- all_probs_df %>% left_join(all_probs_df_trim) %>% dplyr::select(time, age_group, prop_h, prop_h_out, prob_ICU_cases)
    
    ################
    ## AMBULATORY ##
    ################
    # Get the hospital distributions over calendar months
    calmonthdist <- map(1:12, ~(resceudata$observationalData[ , .x])/sum(resceudata$observationalData[ , .x]))
    # Get total ambulatory cases under 6 months of age assumptions (given in excel documents)
    tot_amb_u6mo <- sum(resceudata$observationalData[ , 1:6]) * 5
    # Age distribution of abulatory events
    dist_amb_age <- tot_amb_u6mo * c(0.0689, 0.1248, 0.1763, 0.1601, 0.1956, 0.2743)

    dist_amb_age_month_list <- list()
    for (a in 1:6) { 
        dist_amb_age_month_list[[a]] <- data.frame(
            amb_inci = dist_amb_age[a] * calmonthdist[[a]],
            age_group = a,
            time = 1:12
        )
    }
    for (a in 7:16) { 
        dist_amb_age_month_list[[a]] <- data.frame(
            amb_inci  = resceudata$observationalData[ ,a] * 12,
            age_group = a,
            time = 1:12
        )
    }
    for (a in 17:25) { 
        dist_amb_age_month_list[[a]] <- data.frame(
            amb_inci  = resceudata$observationalData[ ,a] * 0,
            age_group = a,
            time = 1:12
        )
    }
    dist_amb_age_month <- dist_amb_age_month_list %>% bind_rows

    all_probs_df <- all_cases_mean %>% left_join(dist_amb_age_month) %>% mutate(prop_am = amb_inci / inci_all)
    all_probs_df <- all_probs_df %>% mutate(prop_am = case_when(age_group >= 24~0.31, age_group < 24~prop_am))
    all_probs_df_trim <- all_probs_df %>% left_join(all_probs_df_trim) %>% dplyr::select(time, age_group, prop_h, prop_h_out, prob_ICU_cases, prop_am)




    # Deaths
    raw_death_age <- read.csv(file = here("inst", "extdata", "death_raw.csv"))[, 1]
    raw_death <- read.csv(file = here("inst", "extdata", "death_raw.csv"))[, 2:13]

    rates_raw_death <- bind_rows(
        raw_death[rep(1, 6), ],
        raw_death[rep(2, 6), ],
        raw_death[rep(3, 4), ],
        raw_death[rep(4, 7), ],
        (raw_death[rep(5, 1), ] + raw_death[rep(6, 1), ]) / 2,
        (raw_death[rep(7, 1), ] + raw_death[rep(8, 1), ]) / 2,
    )
    rates_real_death <- bind_cols(rates_raw_death[, 7:12], rates_raw_death[, 1:6])

    rates_real_death_list <- list()
    for (i in 1:25) {
        if (i < 25) {
            rates_real_death_list[[i]] <- rates_real_death[i, ] / 1000 * 100000 * (ageGroupBoundary[i + 1] - ageGroupBoundary[i]) / 12
        }
        else {
            rates_real_death_list[[i]] <- rates_real_death[i, ] / 1000 * (population - 100000 * ageGroupBoundary[i-1]) / 12
        }
    }
    incidence_death <- rates_real_death_list %>% bind_rows %>% t %>% as.data.frame


    death_cases <- incidence_death %>%
        pivot_longer(everything(), values_to = "death_prop", names_to = "cols") %>%
        mutate(age_group = rep(1:25, 12), time = unlist(1:12 %>% map(~rep(.x, 25)))) %>%
        dplyr::select(age_group, time, death_prop)

    all_probs_df <- all_cases_mean %>% left_join(death_cases) %>% mutate(prop_d = death_prop / inci_all)

    all_probs_df_trim <- all_probs_df %>% left_join(all_probs_df_trim) %>% dplyr::select(time, age_group, prop_h, prop_h_out, prob_ICU_cases, prop_am, prop_d)
    all_probs_outcome <- all_probs_df_trim
    save(all_probs_df_trim, file = here("data", "all_probs_outcome.RData"))

    outcome_probs <- all_probs_df_trim %>% group_by(age_group) %>%
        summarise(prop_h = mean(prop_h, na.rm = TRUE), prop_h_out = mean(prop_h_out, na.rm = TRUE), prob_ICU_cases = mean(prob_ICU_cases, na.rm = TRUE), 
            prop_am = mean(prop_am, na.rm = TRUE), prop_d = mean(prop_d, na.rm = TRUE)) %>%
        as.data.frame
    save(outcome_probs, file = here("data", "outcome_probs.RData"))
}

calculate_discount_death <- function(life_expectancy = 83) {
    require(stats)

    f <- function(x) {
        exp(-0.03 * x)
    }
    age_groups <- c(0:11/12, 1:5, 10, 15, 25, 35, 45, 55, 65, 75) 
    discount_life_exp <- age_groups %>% map(~integrate(f, 0, life_expectancy - .x)[[1]]) %>% unlist
    # [1] 30.569668 30.562750 30.555815 30.548862 30.541892 30.534905 30.527900 30.520878 30.513838 30.506780
#[11] 30.499705 30.492612 30.485502 30.398772 30.309402 30.217309 30.122412 29.602775 28.999043 27.482653
#[21] 25.435741 22.672699 18.942983 13.908392  7.112405
}
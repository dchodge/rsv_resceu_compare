#' Generate a daily vaccination calendar for Palivizumab administration
#'
#' @param pal_info A list of parameters assocaited with Palivizumab administration
#' @return A list, first element is the when the doses are given, the second is when individuals seroconvert.
get_pal_calendar <- function(pal_info) { 

    calendar <- matrix(0, 365, 25)
    sero <- matrix(0, 365, 25)

    for (i in (pal_info$t_start * 7):(pal_info$t_end * 7)) {
        for (j in 1:25) {
           calendar[(i - 1) %% 365 + 1, j] <- pal_info$age_id[j] * pal_info$cov / 30.0
           sero[(i - 1) %% 365 + 1, j] <- calendar[(i - 1) %% 365 + 1, j] * pal_info$eff;
        }
    }
    list(calendar = calendar, sero = sero)
}

#' Generate a daily vaccination calendar for long-acting monnoclonal ab administration for VHR infants
#'
#' @param mab_vhr_info A list of parameters associated with long-acting monoclonal administration
#' @return A list, first element is the when the doses are given, the second is when individuals seroconvert.
get_vhr_mAB_calendar <- function(mab_vhr_info) { 

    calendar <- matrix(0, 365, 25)
    sero <- matrix(0, 365, 25)

    for (i in (mab_vhr_info$t_start * 7):(mab_vhr_info$t_end * 7)) {
        for (j in 1:25) {
           calendar[(i - 1) %% 365 + 1, j] <- mab_vhr_info$age_id[j] * mab_vhr_info$cov / 30.0
           sero[(i - 1) %% 365 + 1, j] <- calendar[(i - 1) %% 365 + 1, j] * mab_vhr_info$eff;
        }
    }
    list(calendar = calendar, sero = sero)
}

#' Generate a daily vaccination calendar for long-acting monnoclonal ab administration for HR or LR infants
#'
#' @param mab_info A list of parameters associated with long-acting monoclonal administration
#' @return A list, first element is the when the doses are given, the second is when individuals seroconvert.
get_mAB_calendar <- function(mab_info) { 

    calendar <- matrix(0.0, 365, 25)
    sero <- matrix(0.0, 365, 25)

    for (i in (mab_info$t_start * 7):(mab_info$t_end * 7)) {
        for (j in 1:25) {
           calendar[(i - 1) %% 365 + 1, j] <- mab_info$age_id[j] * mab_info$cov / 30.0
        }
    }
    if (mab_info$catchup) {
        for (i in (mab_info$t_start * 7):(mab_info$t_start * 7 + 28)) {
            for (j in 1:25) {
                calendar[(i - 1) %% 365 + 1, j] <- mab_info$age_id_catchup[j] * mab_info$cov / 30.0
            }
        }
    }
    for (i in (mab_info$t_start * 7):(mab_info$t_end * 7)) {
        for (j in 1:25) {
           sero[(i - 1) %% 365 + 1, j] <- calendar[(i - 1) %% 365 + 1, j] * mab_info$eff;
        }
    }

    list(calendar = calendar, sero = sero)
}

#' Generate a daily vaccination calendar for normal vaccination administration for HR or LR infants
#'
#' @param lav_info A list of parameters associated with lav administration
#' @return A list, first element is the when the doses are given, the second is when individuals seroconvert.
get_LAV_calendar <- function(lav_info) { 
    
    calendar <- matrix(0, 365, 25)
    sero <- matrix(0, 365, 25)

    if (lav_info$uptake_type == 1)
    {
        for (i in 1:365)
            for (j in 1:25)
                calendar[(i - 1) %% 365 + 1, j] <- lav_info$uptake_daily[i] * lav_info$age_id[j] * lav_info$cov;
    }
    else if (lav_info$uptake_type == 2)
    {
        for (i in (lav_info$t_start * 7):(lav_info$t_end * 7))
            for (j in 1:25)
                calendar[(i - 1) %% 365 + 1, j] <- lav_info$age_id[j] * lav_info$cov / 30.0;
    } else {
        stop("Value for uptake_type not recognised.")
    }
    
    for (i in 1:365)
        for (j in 1:25)
            for(k in 1:30)
                sero[(i + k - 1) %% 365 + 1, j] <- sero[(i + k - 1) %% 365 + 1, j] +
                    dweibull(k, 2.42, 12.87) * calendar[(i - 1) %% 365 + 1, j] * lav_info$eff;
    
    list(calendar = calendar, sero = sero)

}

#' Generate a daily vaccination calendar for maternal vaccination administration.
#'
#' @param mat_info A list of parameters associated with maternal administration administration
#' @return A list, first element is the when the doses are given, the second is when individuals seroconvert.
get_mat_calendar <- function(mat_info) {

    calendar <- matrix(0, 365, 25)
    sero <- matrix(0, 365, 25)

    for (i in (mat_info$t_start * 7):(mat_info$t_end * 7))
        for (j in 1:25)
            calendar[(i - 1) %% 365 + 1, j] <- mat_info$age_id[j] / 365.0
    

    for(i in 1:365) {
        for(j in 1:25) {
            for(k in 1:30) {
                sero[(i + k - 1) %% 365 + 1, j] <- sero[(i + k - 1) %% 365 + 1, j] + dweibull(k, 2.42, 12.87) * calendar[(i - 1) %% 365 + 1, j] * mat_info$eff_mat;
            }
        }
        for(k in (8*7):(12*7)) {
            sero[(i + k - 1) %% 365 + 1, 1] <- sero[(i + k - 1) %% 365 + 1, 1] + dunif(k, 8 * 7, 12 * 7) * calendar[(i - 1) %% 365 + 1, 20] * mat_info$eff_inf * 365;
        }
        
    }
    list(calendar = calendar, sero = sero)
}

#' Function to convert a custom weekly uptake vector into a daily uptake calendar
#'
#' @param up_week_raw A vector of weekly uptake proportions
#' @param start_time Start week when administration occurs
#' @return A vector of the daily proportion of a target group which are vaccinated.
get_daily_uptake <- function(up_week_raw, start_time)
{
    up_week <- vector("numeric", length = 52)
    up_day <- vector("numeric", length = 365)
    for (k in (start_time):(start_time + 20)) { 
        pos_week <- (k - 1) %% 52 + 1
        pos_ist <- (k - start_time) + 1
        up_week[pos_week] <- up_week_raw[pos_ist]
    }
            
    for (i in 1:365) {
        pos_day <- (i - 1 + start_time*7) %% 365 + 1;
        up_day[pos_day] <- (up_week[((i / 7 + start_time) %% 52) + 1] - up_week[((i/7 - 1 + start_time) %% 52) + 1]) / 7.0;
        if (up_day[pos_day] < 0)
            up_day[pos_day] <- 0
    }
    up_day
}

#' Function to check the custom entries of the vaccination calendars.
#'
#' @param info A list of paramter values
check_all_entries <- function(info) {
    if (is.null(info$age_id))
        stop("PARAMETER MISSING: Binary vector indicating target age groups not provided. (age_id)")
    if (is.null(info$t_start))
        stop("PARAMETER MISSING: Start time of vaccination programme not provided. (t_start)")
    if (is.null(info$t_end))
        stop("PARAMETER MISSING: End time of vaccination programme not provided. (t_end)")
    if (is.null(info$eff))
        stop("PARAMETER MISSING: Efficacy not provided. (eff)")
    if (is.null(info$cov))
        stop("PARAMETER MISSING: Coverage not provided. (cov)")
}

#' Function to check the custom entries, specific to maternal vaccination calendars.
#'
#' @param info A list of paramter values
check_mat_entries <- function(info) {
    if (is.null(info$age_id))
        stop("PARAMETER MISSING: Binary vector indicating target age groups not provided. (age_id)")
    if (is.null(info$t_start))
        stop("PARAMETER MISSING: Start time of vaccination programme not provided. (t_start)")
    if (is.null(info$t_end))
        stop("PARAMETER MISSING: End time of vaccination programme not provided. (t_end)")
    if (is.null(info$eff_inf))
        stop("PARAMETER MISSING: Efficacy for infants not provided. (eff_inf)")
    if (is.null(info$eff_mat))
        stop("PARAMETER MISSING: Efficacy for pregnant women not provided. (eff_mat)")
    if (is.null(info$cov))
        stop("PARAMETER MISSING: Coverage not provided. (cov)")
}

#' Function to check the custom entries, specific to long-acting monoclonal antibody calendars.
#'
#' @param info A list of paramter values
check_mab_entries <- function(info) {
    if (is.null(info$catchup))
        stop("PARAMETER MISSING: Must indicate boolean for including a catchup programme at start of season. (catchup)")
    else {
        if(info$catchup) {
            if (is.null(info$age_id_catchup))
                stop("PARAMETER MISSING: Must provide a binary vector indicating target age groups for the catchup programme. (age_id_catchup)")
        }
    }
}

#' Function to check the custom entries, specific to lav calendars.
#'
#' @param info A list of paramter value
check_lav_entries <- function(info) {
    if (is.null(info$age_id))
        stop("PARAMETER MISSING: Binary vector indicating target age groups not provided. (age_id)")
    if (is.null(info$uptake_daily))
        stop("PARAMETER MISSING: Daily uptake proportions missing. (uptake_daily)")
    if (is.null(info$eff))
        stop("PARAMETER MISSING: Efficacy not provided. (eff)")
    if (is.null(info$cov))
        stop("PARAMETER MISSING: Coverage not provided. (cov)")
}

#' Function to convert input parameters into a daily uptake calendar for each age group in the transmission model
#'
#' @param vac_program_info A list of details outlining the parmeters for each of the 7 calendars
#' @param metric The metric (default week)
#' @param cov_c The proportion of infants who are in the targeted group (equal to the coverage of maternal vaccination.)
#' @return Daily vaccination calendars for use in the dynamic transmission
create_calendar <- function(vac_program_info, metric = "week", cov_c) {
    # need to include week and month options
    all_cal <- list(
        pal = matrix(0, 365, 25),
        mAB_VHR = matrix(0, 365, 25),
        mAB_HR = matrix(0, 365, 25),
        mAB_LR = matrix(0, 365, 25),
        LAV_HR = matrix(0, 365, 25),
        LAV_LR = matrix(0, 365, 25),
        mat_LR = matrix(0, 365, 25)
    )   
    if(!is.null(vac_program_info$pal$id) && vac_program_info$pal$id) {
        check_all_entries(vac_program_info$pal)
        all_cal$pal <- get_pal_calendar(vac_program_info$pal)
    }
    if(!is.null(vac_program_info$mAB_VHR$id) && vac_program_info$mAB_VHR$id ) {
        check_all_entries(vac_program_info$mAB_VHR)
        all_cal$mAB_VHR <- get_vhr_mAB_calendar(vac_program_info$mAB_VHR)
    }
    if(!is.null(vac_program_info$mAB_HR$id) && vac_program_info$mAB_HR$id ) {
        check_all_entries(vac_program_info$mAB_HR)
        check_mab_entries(vac_program_info$mAB_HR)
        all_cal$mAB_HR <- get_mAB_calendar(vac_program_info$mAB_HR)
    }
    if(!is.null(vac_program_info$mAB_LR$id) && vac_program_info$mAB_LR$id ) {
        check_all_entries(vac_program_info$mAB_LR)
        check_mab_entries(vac_program_info$mAB_LR)
        all_cal$mAB_LR <- get_mAB_calendar(vac_program_info$mAB_LR)
    }
    if(!is.null(vac_program_info$LAV_HR$id) && vac_program_info$LAV_HR$id ) {
        if (vac_program_info$LAV_HR$uptake_type == 2) {
            check_all_entries(vac_program_info$LAV_HR)
        } else {
            check_lav_entries(vac_program_info$LAV_HR)
        }
        all_cal$LAV_HR <- get_LAV_calendar(vac_program_info$LAV_HR)
    }
    if(!is.null(vac_program_info$LAV_LR$id) && vac_program_info$LAV_LR$id) {
        if (vac_program_info$LAV_LR$uptake_type == 2) {
            check_all_entries(vac_program_info$LAV_LR)
        } else {
            check_lav_entries(vac_program_info$LAV_LR)
        }
        all_cal$LAV_LR <- get_LAV_calendar(vac_program_info$LAV_LR)
    }
    if(!is.null(vac_program_info$mat$id) && vac_program_info$mat$id ) { 
        check_mat_entries(vac_program_info$mat)
        all_cal$mat_LR <- get_mat_calendar(vac_program_info$mat)
    }


    list(all_cal = all_cal, cov_c = cov_c)
}

#' Function to get the serorevertion times given an output from create_calendar
#' 
#' @param output A list outputs from create_calendar
#' @return A list of matrices showing the seroconversion proportions for each of the seven vaccination calendars.
get_cal <- function(output) {
    names <- c("pal", "mAB_VHR", "mAB_HR", "mAB_LR", "LAV_HR", "LAV_LR", "mat_LR")
    out <- vector(mode = "list", length = length(output[[1]])) %>% setNames(names)
    j <- 1
    for (i in output[[1]]) {
        if (is.matrix(i))
            out[[j]] <- i
        else
            out[[j]] <- i[[2]]
        j <- j + 1
    }
    out
}

#' Function to get the dosing times given an output from create_calendar
#' 
#' @param output A list outputs from create_calendar
#' @return A list of matrices showing the dosing proportions for each of the seven vaccination calendars.
get_dose <- function(output) {
    names <- c("pal", "mAB_VHR", "mAB_HR", "mAB_LR", "LAV_HR", "LAV_LR", "mat_LR")
    out <- vector(mode = "list", length = length(output[[1]])) %>% setNames(names)
    j <- 1
    for (i in output[[1]]) {
        if (is.matrix(i))
            out[[j]] <- i
        else
            out[[j]] <- i[[1]]
        j <- j + 1
    }
    out
}

#' Function to run a custom vaccination calendar
#' 
#' @param seeds A vector of seed value (integers)
#' @param func_vac A function which generate the vaccination calendar given a seed value
#' @param vac_par_info A list of values which help define the intervention programmes
#' @param cov_c The proportion of infants who are in the targeted group (equal to the coverage of maternal vaccination.)
#' @return A list of two dataframes, the first shows the annual incidence for health outcomes, the second shows the economic metrics.
run_sample_custom <- function(seeds, func_vac, vac_par_info, cov_c, post, cost_imp, discount_rate = 0.03, disease_mult = 1) {
    QALY_list <- list(); cost_list <- list();
    outcomes_month_age_list <- list(); econmetric_list <- list(); outcomes_annual_list <- list(); i <- 1;
    for (seed in seeds) {
    
      vac_program_info_custom <- func_vac(seed)

      vac_cal <- create_calendar(vac_program_info_custom, "week", cov_c)
      cal <- get_cal(vac_cal)
      dose <- get_dose(vac_cal)
      out <- classRunInterventions$SampleMonthly(cal, dose, vac_cal[["cov_c"]],
        vac_par_info, as.matrix(post)[seed, ])
      outcomes_cea <- get_outcomes(out, as.matrix(post)[seed, ], cost_imp, discount_rate, seed, disease_mult)
      
      QALY_list[[i]] <- outcomes_cea$QALY
      cost_list[[i]] <- outcomes_cea$cost
      outcomes_month_age_list[[i]] <- outcomes_cea$outcomes_age_month

      cat("sample_no: ", i, "\n"); i <- i + 1;
    }
    list(
        QALY = bind_rows(QALY_list),
        cost = bind_rows(cost_list),
        outcomes_month_age = bind_rows(outcomes_month_age_list),
        vac_cal = out$doses,
        cal_pre = cal,
        dose_pre = dose
    )
}



get_output_format <- function(outcomes_age_month, outcomes_name) {
    w <- c(rep(1 / 12, 12), 1, 1, 1, 1, 5, 5, 10, 10, 10, 10, 10, 10, 10)

    outputA <- outcomes_age_month %>%
        filter(outcome == outcomes_name) %>%
        group_by(time, age_group) %>%
        summarise(incidence = mean(incidence, na.rm = TRUE)) %>%
        pivot_wider(everything(), names_from = "time", values_from = "incidence") %>%
        setNames(c("age_group", "jul", "aug", "sep", "oct", "nov", "dec", "jan", "feb", "mar", "apr", "may", "jun")) %>%
        dplyr::select(age_group, jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec)

    # Under 5s
    outputu5 <- bind_rows(
        c(outputA[1:3, 2:13] %>% apply(2, sum)), # 0–2 months   
        c(outputA[4:6, 2:13] %>% apply(2, sum)),  # 3–5 months
        c(outputA[7:12, 2:13] %>% apply(2, sum)),   # 6–11 months
        c(outputA[13, 2:13]),                       # 1 year
        c(outputA[14:16, 2:13] %>% apply(2, sum)),  # 2–4 years
        c(outputA[1:12, 2:13] %>% apply(2, sum)),  # <1 years 
        c(outputA[1:16, 2:13] %>% apply(2, sum))  # 0-4 years
    )

    outputu5_final <- bind_cols(c("0-3mo", "4-6mo", "7-11mo", "12-23mo", "24-59mo", "0-12mo", "0-59mo"), outputu5)
    
    # over 5s
    outputo5 <- bind_rows(
        c(bind_rows(outputA[17:22, 2:13], outputA[23, 2:13] / 2) %>% apply(2, sum)), # 5-59 yrs   
        c(outputA[23, 2:13] / 2),  # 60–64 yrs
        c(outputA[24, 2:13] / 2),   # 65–69 yrs
        c(outputA[24, 2:13] / 2),   # 70–74 yrs
        c(outputA[25, 2:13] / 3),   # 75-79 yrs
        c(outputA[25, 2:13] / 3),   # 80-84 yrs 
        c(outputA[25, 2:13] / 3),   # 85+ yrs 
        c(bind_rows(outputA[23, 2:13] / 2, outputA[24:25, 2:13]) %>% apply(2, sum)) # 60 + 
    )

    outputo5_final <- bind_cols(c("5-59yrs", "60-64yrs", "65-69yrs", "70-74yrs", "75-79yrs", "80-84yrs", "85+yrs", "60+yrs"), outputo5)

    list(under5 = outputu5_final, over5 = outputo5_final)

}

get_output_format_full <- function(outcomes_age_month, outcomes_name) {
    w <- c(rep(1 / 12, 12), 1, 1, 1, 1, 5, 5, 10, 10, 10, 10, 10, 10, 10)

    outputA <- outcomes_age_month %>%
        filter(outcome == outcomes_name) %>%
        group_by(time, age_group) %>%
        summarise(incidence = mean(incidence, na.rm = TRUE)) %>%
        pivot_wider(everything(), names_from = "time", values_from = "incidence") %>%
        setNames(c("age_group", "jul", "aug", "sep", "oct", "nov", "dec", "jan", "feb", "mar", "apr", "may", "jun")) %>%
        dplyr::select(age_group, jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec)

    outputA

}

qaly_output <- function(QALY) {
  QALY_mean <- QALY %>%
    group_by(age_group, type, metric) %>%
    summarise(value = mean(value, na.rm = TRUE))
    # under 5s
  df_u5_qaly <- bind_rows(
    QALY_mean %>% filter(age_group %in% 1:3) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "0-3m"),
    QALY_mean %>% filter(age_group %in% 4:6) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "4-6m"),
    QALY_mean %>% filter(age_group %in% 7:12) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "7-11m"),
    QALY_mean %>% filter(age_group %in% 13) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "12-23m"),
    QALY_mean %>% filter(age_group %in% 14:16) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "24-59m"),
    QALY_mean %>% filter(age_group %in% 1:12) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "0-11m"),
    QALY_mean %>% filter(age_group %in% 1:16) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "0-59m")
  )

  # over 5s
  df_o5_qaly <- bind_rows(
    bind_rows(
      filter(QALY_mean, age_group %in% 17:22),
      mutate(filter(QALY_mean, age_group == 23), value = value / 2)
    ) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "5-59yrs"),
    QALY_mean %>% filter(age_group %in% 23) %>% group_by(type, metric) %>% summarise(value = value / 2) %>% mutate(age_group = "60-64yrs"),
    QALY_mean %>% filter(age_group %in% 24) %>% group_by(type, metric) %>% summarise(value = value / 2) %>% mutate(age_group = "65-69yrs"),
    QALY_mean %>% filter(age_group %in% 24) %>% group_by(type, metric) %>% summarise(value = value / 2) %>% mutate(age_group = "70-74yrs"),

    QALY_mean %>% filter(age_group %in% 25) %>% group_by(type, metric) %>% summarise(value = value / 3) %>% mutate(age_group = "75-79yrs"),
    QALY_mean %>% filter(age_group %in% 25) %>% group_by(type, metric) %>% summarise(value = value / 3) %>% mutate(age_group = "80-84yrs"),
    QALY_mean %>% filter(age_group %in% 25) %>% group_by(type, metric) %>% summarise(value = value / 3) %>% mutate(age_group = "85+yrs"),
    bind_rows(
      mutate(filter(QALY_mean, age_group == 23), value = value / 2),
      filter(QALY_mean, age_group %in% 24:25)
    ) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "60+yrs")
  )

  list(
    under5 = (df_u5_qaly %>% unite("metric", type, metric) %>% pivot_wider(everything(), names_from = metric, values_from = value)),
    over5 = (df_o5_qaly %>% unite("metric", type, metric) %>% pivot_wider(everything(), names_from = metric, values_from = value))
  )
}

cost_output <- function(cost) {
  cost_mean <- cost %>%
    group_by(age_group, type, metric) %>%
    summarise(value = mean(value, na.rm = TRUE))
    # under 5s
  df_u5_cost <- bind_rows(
    cost_mean %>% filter(age_group %in% 1:3) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "0-3m"),
    cost_mean %>% filter(age_group %in% 4:6) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "4-6m"),
    cost_mean %>% filter(age_group %in% 7:12) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "7-11m"),
    cost_mean %>% filter(age_group %in% 13) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "12-23m"),
    cost_mean %>% filter(age_group %in% 14:16) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "24-59m"),
    cost_mean %>% filter(age_group %in% 1:12) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "0-11m"),
    cost_mean %>% filter(age_group %in% 1:16) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "0-59m")
  )

  # over 5s
  df_o5_cost <- bind_rows(
    bind_rows(
      filter(cost_mean, age_group %in% 17:22),
      mutate(filter(cost_mean, age_group == 23), value = value / 2)
    ) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "5-59yrs"),
    cost_mean %>% filter(age_group %in% 23) %>% group_by(type, metric) %>% summarise(value = value / 2) %>% mutate(age_group = "60-64yrs"),
    cost_mean %>% filter(age_group %in% 24) %>% group_by(type, metric) %>% summarise(value = value / 2) %>% mutate(age_group = "65-69yrs"),
    cost_mean %>% filter(age_group %in% 24) %>% group_by(type, metric) %>% summarise(value = value / 2) %>% mutate(age_group = "70-74yrs"),

    cost_mean %>% filter(age_group %in% 25) %>% group_by(type, metric) %>% summarise(value = value / 3) %>% mutate(age_group = "75-79yrs"),
    cost_mean %>% filter(age_group %in% 25) %>% group_by(type, metric) %>% summarise(value = value / 3) %>% mutate(age_group = "80-84yrs"),
    cost_mean %>% filter(age_group %in% 25) %>% group_by(type, metric) %>% summarise(value = value / 3) %>% mutate(age_group = "85+yrs"),
    bind_rows(
      mutate(filter(cost_mean, age_group == 23), value = value / 2),
      filter(cost_mean, age_group %in% 24:25)
    ) %>% group_by(type, metric) %>% summarise(value = sum(value)) %>% mutate(age_group = "60+yrs")
  )

  list(
    under5 = (df_u5_cost %>% unite("metric", type, metric) %>% pivot_wider(everything(), names_from = metric, values_from = value)),
    over5 = (df_o5_cost %>% unite("metric", type, metric) %>% pivot_wider(everything(), names_from = metric, values_from = value))
  )
}


save_outputs <- function(sample_output, filename, excel_format = TRUE) {

    if (excel_format) {
        cases_outputs <- list(
            asymp = get_output_format(sample_output$outcomes_month_age, "asymptomtic_cases"),
            symp_nhc = get_output_format(sample_output$outcomes_month_age, "nonhealthcare_seeking_symptomatic_cases"),
            pri = get_output_format(sample_output$outcomes_month_age, "pri_cases"),
            sec = get_output_format(sample_output$outcomes_month_age, "sec_cases"),
            hosp = get_output_format(sample_output$outcomes_month_age, "hosp_cases"),
            icu = get_output_format(sample_output$outcomes_month_age, "picu_cases"),
            death = get_output_format(sample_output$outcomes_month_age, "death_cases")
        )

        QALY_outputs <- sample_output$QALY %>% qaly_output
        cost_outputs <- sample_output$cost %>% cost_output

        output_final <- list(cases = cases_outputs,
                QALY = QALY_outputs,
                cost = cost_outputs
        )

        save(output_final, file = here("outputs", "burden", paste0(filename, ".RData")))
        output_final
    } else {
        cases_outputs <- list(
            asymp = get_output_format_full(sample_output$outcomes_month_age, "asymptomtic_cases"),
            symp_nhc = get_output_format_full(sample_output$outcomes_month_age, "nonhealthcare_seeking_symptomatic_cases"),
            pri = get_output_format_full(sample_output$outcomes_month_age, "pri_cases"),
            sec = get_output_format_full(sample_output$outcomes_month_age, "sec_cases"),
            hosp = get_output_format_full(sample_output$outcomes_month_age, "hosp_cases"),
            icu = get_output_format(sample_output$outcomes_month_age, "picu_cases"),
            death = get_output_format_full(sample_output$outcomes_month_age, "death_cases")
        )
        cases_outputs

    }
}

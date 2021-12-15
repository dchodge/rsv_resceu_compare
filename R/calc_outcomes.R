#' Function to convert incidence to the number of symptomatic cases
#' 
#' @param inci A matrix of incidence in each age, social, and risk group.
#' @param seed An integer value
#' @return The number of annual symptomatic cases in each age group vector.
get_Pri <- function(inci, t_w, seed) {
    load(file = here("data", "outcome_probs.RData"))
    set.seed(seed)
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_am
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_am 
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_am 

    data.frame(
        outcome = "pri_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}

get_PICU <- function(inci, t_w, seed) {
    load(file = here("data", "outcome_probs.RData"))
    set.seed(seed)
    
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prob_ICU_cases
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prob_ICU_cases 
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prob_ICU_cases 

    data.frame(
        outcome = "picu_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}

get_S <- function(inci, t_w, posterior, seed) {
    set.seed(seed)
    risk_lr <- vector(mode = "numeric", length = 25)
    for (a in 1:12)
        risk_lr[a] <- (1 - 0.0916)
    for (a in 13:16)
        risk_lr[a] <- (1 - 0.163)
    for (a in 17:18)
        risk_lr[a] <- (1 - 0.516)
    for (a in 19:25)
        risk_lr[a] <- (1 - 0.753)
    
    risk_vhr <- risk_hr <- risk_lr
    
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist) * risk_vhr
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) * risk_hr
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist) * risk_lr

    data.frame(
        outcome = "symptomatic_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}

get_Inc <- function(inci, t_w, seed) {

    
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist)
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) 
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist)

    data.frame(
        outcome = "all_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}

#' Function to convert incidence to GP visits
#' 
#' @param inci A matrix of incidence in each age, social, and risk group.
#' @param seed An integer value
#' @return The number of annual GP cases in each age group vector.
get_Sec <- function(inci, t_w, seed) {
    load(file = here("data", "outcome_probs.RData"))
    set.seed(seed)
    
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_h_out 
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_h_out 
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_h_out 

    data.frame(
        outcome = "sec_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}



#' Function to convert incidence to hospital visits
#' 
#' @param inci A matrix of incidence in each age, social, and risk group.
#' @param seed An integer value
#' @return The number of annual hospital cases in each age group vector.
get_Hosp <- function(inci, t_w, seed) {
    load(file = here("data", "outcome_probs.RData"))
    set.seed(seed)
    
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_h 
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_h 
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_h 

    data.frame(
        outcome = "hosp_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}


#' Function to convert incidence to the number of deaths
#' 
#' @param inci A matrix of incidence in each age, social, and risk group.
#' @param seed An integer value
#' @return The number of annual death cases in each age group vector.
get_Death <- function(inci, t_w, seed) {
    load(file = here("data", "outcome_probs.RData"))
    set.seed(seed)
    
    out_vhr <- (1:25 %>% map(~sum(inci[(c(1, 4, 7) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_d
    out_hr <- (1:25 %>% map(~sum(inci[(c(2, 5, 8) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_d
    out_lr <- (1:25 %>% map(~sum(inci[(c(3, 6, 9) + 9 * (.x - 1))])) %>% unlist) * outcome_probs$prop_d

    data.frame(
        outcome = "death_cases",
        seed = seed,
        time = t_w,
        age_group = 1:25,
        incidence = out_vhr + out_hr + out_lr
    )
}

#' Function to convert incidence of outcomes into the QALY loss
#' 
#' @param symp A vector of number of symptomatic cases per age group
#' @param gp A vector of number of GP cases per age group
#' @param hosp A vector of number of hospital cases per age group
#' @param bd A vector of number of bed days per age group
#' @param death A vector of number of deaths per age group
#' @param seed An integer value
#' @return The total number of QALYs lost
get_QALY <- function(symp, pri, sec, hosp, death, seed) {
    set.seed(seed)
    nhseek_QALY <-  vector(mode = "numeric", length = 25)
    hseek_QALY <-  vector(mode = "numeric", length = 25)
    death_QALY <-  vector(mode = "numeric", length = 25)

    for (a in 1:12) 
        death_QALY[a] <- rnorm(1, 30.569668, 30.569668 * 0.1)

    death_QALY[13] <- rnorm(1, 30.485502, 30.485502 * 0.1)
    death_QALY[14] <- rnorm(1, 30.398772, 30.398772 * 0.1)
    death_QALY[15] <- rnorm(1, 30.309402, 30.309402 * 0.1)
    death_QALY[16] <- rnorm(1, 30.217309, 30.217309 * 0.1)
    death_QALY[17] <- rnorm(1, 30.122412, 30.122412 * 0.1)
    death_QALY[18] <- rnorm(1, 29.602775, 29.602775 * 0.1)
    death_QALY[19] <- rnorm(1, 28.999043, 28.999043 * 0.1)
    death_QALY[20] <- rnorm(1, 27.482653, 27.482653 * 0.1)
    death_QALY[21] <- rnorm(1, 25.435741, 25.435741 * 0.1)
    death_QALY[22] <- rnorm(1, 22.672699, 22.672699 * 0.1)
    death_QALY[23] <- rnorm(1, 18.942983, 18.942983 * 0.1)
    death_QALY[24] <- rnorm(1, 13.908392, 13.908392 * 0.1)
    death_QALY[25] <- rnorm(1, 7.112405,  7.112405 * 0.1)

    tot_Q_c <- rep(0, 25)
    tot_Q_d <- rep(0, 25)

    for (a in 1:16) { # <5 years
        nhseek_QALY[a] <- 0.003024
        hseek_QALY[a] <- 0.003823
       # nhseek_QALY[a] <- rgamma(1, 1.6578, scale = 0.0018241) # non-healthcare seeking QALY loss (3.024 × 10−3)
        #hseek_QALY[a] <- rgamma(1, 1.7927, scale = 0.00213254) # healthcare seeking QALY loss (3.823 × 10−3)
        tot_Q_c[a] <-  tot_Q_c[a] + (symp[a] - pri[a] - sec[a] - hosp[a]) * nhseek_QALY[a]
        tot_Q_c[a] <- tot_Q_c[a] + (pri[a] + sec[a] + hosp[a]) * hseek_QALY[a]
        tot_Q_d[a] <- tot_Q_d[a] + death[a] * death_QALY[a];
    }
    for (a in 17:25) { # >5 years
        nhseek_QALY[a] <- 0.001543
        hseek_QALY[a] <- 0.001950
       # nhseek_QALY[a] <- rgamma(1, 1.36973, scale = 0.0011265)  # non-healthcare seeking QALY loss (1.543 × 10−3)
       # hseek_QALY[a] <- rlnorm(1, -6.23993, 0.933905) # healthcare seeking QALY loss (1.950 × 10−3)
        tot_Q_c[a] <- tot_Q_c[a] + (symp[a] - pri[a] - sec[a] - hosp[a]) * nhseek_QALY[a]
        tot_Q_c[a] <- tot_Q_c[a] + (pri[a] + sec[a] + hosp[a]) * hseek_QALY[a]
        tot_Q_d[a] <- tot_Q_d[a] + death[a] * death_QALY[a];
    }

    list(qaly_cases = tot_Q_c,
        qaly_death = tot_Q_d,
        qaly_total = tot_Q_c + tot_Q_d
    )
}

#' Function to convert incidence of outcomes into the cost of treatment
#' 
#' @param gp A vector of number of GP cases per age group
#' @param bd A vector of number of bed days per age group
#' @param seed An integer value
#' @return The total cost of treatment  
get_costT <- function(pri, sec, hosp, picu, seed) {
    set.seed(seed)
     
    costT <- rep(0, 25)
    pri_c <- 33
    sec_c <- 104
    hosp_days <- 5.8
    hosp_c_u5 <- hosp_days *  627
    picu_days <- 8.1
    PICU_c <- picu_days *  2015 + 613
   # hosp_c_u5 <- rgamma(1, 4.8 / 5.8, scale = (5.8) ^ 2 / 4.8) * 627
    hosp_c_o5 <- 4658
   # PICU_c <- rgamma(1, 8.0 / 8.1, scale = (8.1) ^ 2 / 8.0) * 2015 + 613

    for (a in 1:25) {
        costT[a] <- costT[a] + pri[a] * pri_c;
        costT[a] <- costT[a] + sec[a] * sec_c;
        costT[a] <- costT[a] + picu[a] * PICU_c;
    }
    
    for (a in 1:16)
        costT[a] <- costT[a] + hosp[a] * hosp_c_u5
    
    for (a in 17:25)
        costT[a] <- costT[a] + hosp[a] * hosp_c_o5

    costT
}

#' Function to convert dose schedule into the cost of administration
#' 
#' @param doses_w A matrix of doses
#' @return The total cost of administration
#' # pal, mab, mat, lav
get_costA <- function(doses_w) {
    
    pal_dose <- doses_w[seq(1, 100, 4)]
    mab_dose <- doses_w[seq(2, 100, 4)]
    lav_dose <- doses_w[seq(3, 100, 4)]
    mat_dose <- doses_w[seq(4, 100, 4)]

    pal_c <- 51.82 * 5.08
    mab_c <- 8.32
    lav_c <- 11.36
    mat_c <- 5

    costA <- pal_dose * pal_c + mab_dose * mab_c + lav_dose * lav_c + mat_dose * mat_c
    costA
}

get_costP <- function(doses_w) {
    
    pal_dose <- doses_w[seq(1, 100, 4)]
    mab_dose <- doses_w[seq(2, 100, 4)]
    lav_dose <- doses_w[seq(3, 100, 4)]
    mat_dose <- doses_w[seq(4, 100, 4)]

    pal_c <- 928.6 * 5.08
    mab_c <- 50
    lav_c_y <- 60
    lav_c_o <- 40

    mat_c <- 37.5

    costP <- pal_dose * pal_c + mab_dose * mab_c + c(lav_dose[1:16] * lav_c_y, lav_dose[17:25] * lav_c_o) + mat_dose * mat_c
    costP
}

get_costIndirect <- function(symp, sec, hosp, picu, doses_w, seed) {
    set.seed(seed)

    cost_outpatient_day <- 139
    cost_hosp_day <- 139
    cost_picu_day <- 139
    hosp_days <- 5.8
    picu_days <- 8.1
   # hosp_days <- rgamma(1, 4.8 / 5.8, scale = (5.8) ^ 2 / 4.8)
   # picu_days <- rgamma(1, 8.0 / 8.1, scale = (8.1) ^ 2 / 8.0)
    costIn <- rep(0, 25)

    o5_day_cost <- 36.16
    days_off_inf <- 3.5
    days_off_hosp <- 6.92

    for (a in 1:16) {
        costIn[a] <- costIn[a] + sec[a] * cost_outpatient_day
        costIn[a] <- costIn[a] + hosp[a] * 189 * 0.19
        costIn[a] <- costIn[a] + hosp[a] * hosp_days * cost_hosp_day
        costIn[a] <- costIn[a] + picu[a] * 189 * 0.19
        costIn[a] <- costIn[a] + picu[a] * picu_days * cost_picu_day
    }
    for (a in 17:25) {
        costIn[a] <- costIn[a] + sec[a] * cost_outpatient_day
        costIn[a] <- costIn[a] + hosp[a] * 189 * 0.19
        costIn[a] <- costIn[a] + hosp[a] * hosp_days * cost_hosp_day
        costIn[a] <- costIn[a] + picu[a] * 189 * 0.19
        costIn[a] <- costIn[a] + picu[a] * picu_days * cost_picu_day
        costIn[a] <- costIn[a] + symp[a] * days_off_inf * o5_day_cost
        costIn[a] <- costIn[a] + hosp[a] * days_off_hosp * o5_day_cost
    }
    costIn
}

get_Costs <- function(symp, pri, sec, hosp, picu, death, cost_imp, doses_w, seed) {
    costT <- get_costT(pri, sec, hosp, picu, seed)
    costA <- get_costA(doses_w)
    costP <- get_costP(doses_w)


    cost_direct <- costT
    cost_intervention <- costA + costP
    cost_intervention[1] <- cost_intervention[1] + cost_intervention[19] + cost_intervention[20] + cost_intervention[21] # add maternal costs to first age group
    cost_intervention[19] <- cost_intervention[20] <- cost_intervention[21] <- 0
    cost_indirect <- get_costIndirect(symp, sec, hosp, picu, doses_w, seed)
    cost_total <- cost_direct + cost_intervention + cost_indirect

    list(direct = cost_direct,
        intervention = cost_intervention,
        indirect = cost_indirect,
        total = cost_total
    )
}

#' Function to convert outputs form RunInterventions model into outcomes
#' 
#' @param outputs Output from the RunInterventions model
#' @param posterior Posterior distributions from calibration
#' @param seed integer value
#' @return list of two dataframes, once detailing the health outcomes, one the economic outcomes
get_outcomes <- function(outputs, posterior, cost_imp, discount_rate, seed, disease_mult) {
    inci <- outputs$inci
    doses <- outputs$doses

    undiscount_qaly <- undiscount_cost <- list(rep(0, 25), rep(0, 25), rep(0, 25), rep(0, 25))
    discount_qaly <- discount_cost <- list(rep(0, 25), rep(0, 25), rep(0, 25), rep(0, 25))

    QALY <- 0
    costP <- 0
    costA <- 0
    costT <- 0

    death_tot <- picu_tot <- hosp_tot <- sec_tot <- pri_tot <- symp_nhs_tot <- symp_tot <- asymp_tot <- cases_tot <- 0
    outcomes_age_month <- data.frame()
    r <- 0
    for (t_m in 1:nrow(inci)) {
        inci_tm <- inci[(t_m - 1) %% 12 + 1, ]
        cases <- get_Inc(inci_tm, t_m, seed)
        symp <- get_S(inci_tm, t_m, posterior, seed)
        asymp <- cases %>%
            rename(cases = incidence) %>%
            mutate(symp = symp$incidence, incidence = cases - symp, outcome = "asymptomtic_cases") %>%
            dplyr::select(outcome, seed, time, age_group, incidence)

        pri <- get_Pri(inci_tm, t_m, seed)
        sec <- get_Sec(inci_tm, t_m, seed) %>% mutate(incidence = incidence * disease_mult)
        hosp <- get_Hosp(inci_tm, t_m, seed)  %>% mutate(incidence = incidence * disease_mult)
        picu <- get_PICU(inci_tm, t_m, seed)  %>% mutate(incidence = incidence * disease_mult)
        death <- get_Death(inci_tm, t_m, seed)

        symp_nhs <- symp %>%
            rename(symp = incidence) %>%
            mutate(hcs = pri$incidence + sec$incidence + hosp$incidence, incidence = symp - hcs, outcome = "nonhealthcare_seeking_symptomatic_cases") %>%
            dplyr::select(outcome, seed, time, age_group, incidence)

        if ((t_m >= (1)) & (t_m <= (12))) {
            cases_tot <- cases_tot + sum(cases$incidence)
            symp_tot <- symp_tot + sum(symp$incidence)
            pri_tot <- pri_tot + sum(pri$incidence)
            sec_tot <- sec_tot + sum(sec$incidence)
            hosp_tot <- hosp_tot + sum(hosp$incidence)
            picu_tot <- picu_tot + sum(picu$incidence)
            death_tot <- death_tot + sum(death$incidence)
            outcomes_age <- bind_rows(cases, asymp, symp, symp_nhs, pri, sec, hosp, picu, death)
            outcomes_age_month <- bind_rows(outcomes_age_month, outcomes_age)
        }

        undiscount_qaly_tm <- get_QALY(symp$incidence, pri$incidence, sec$incidence * disease_mult, hosp$incidence * disease_mult, death$incidence, seed)
        discount_qaly_tm <- undiscount_qaly_tm %>% map(~.x * exp(-(t_m - 1) * discount_rate / 12.0))
        undiscount_qaly <- ((1:3 %>% map(~undiscount_qaly[[.x]] + undiscount_qaly_tm[[.x]])) %>% setNames(c("qaly_cases", "qaly_death", "qaly_total")))
        discount_qaly <- ((1:3 %>% map(~discount_qaly[[.x]] + discount_qaly_tm[[.x]])) %>% setNames(c("qaly_cases", "qaly_death", "qaly_total")))
        
        undiscount_cost_tm <- get_Costs(symp$incidence, pri$incidence, sec$incidence * disease_mult, hosp$incidence * disease_mult, picu$incidence * disease_mult, death$incidence, cost_imp, doses[t_m, ], seed)
        discount_cost_tm <- undiscount_cost_tm %>% map(~.x * exp(-(t_m - 1) * discount_rate / 12.0))
        undiscount_cost <- ((1:4 %>% map(~undiscount_cost[[.x]] + undiscount_cost_tm[[.x]])) %>% setNames(c("cost_direct", "cost_inter", "cost_indirect", "cost_total")))
        discount_cost <- ((1:4 %>% map(~discount_cost[[.x]] + discount_cost_tm[[.x]])) %>% setNames(c("cost_direct", "cost_inter", "cost_indirect", "cost_total")))

    }
    undiscount_cost$inter[1] <- undiscount_cost$inter[1] + cost_imp
    discount_cost$inter[1] <- undiscount_cost$inter[1] + cost_imp

    #cost_imp needs to be added
    QALY <- data.frame(
        seed = seed,
        age_group = rep(1:25, 6),
        type = c(rep("undiscounted", 75), rep("discounted", 75)),
        metric = c(rep("cases", 25), rep("deaths", 25), rep("total", 25), rep("cases", 25), rep("deaths", 25), rep("total", 25)),
        value = c(undiscount_qaly$qaly_cases, undiscount_qaly$qaly_death, undiscount_qaly$qaly_total,
            discount_qaly$qaly_cases, discount_qaly$qaly_death, discount_qaly$qaly_total)
    )
    
    cost <- data.frame(
        seed = seed,
        age_group = rep(1:25, 8),
        type = c(rep("undiscounted", 100), rep("discounted", 100)),
        metric = c(rep("direct", 25), rep("inter", 25), rep("indirect", 25), rep("total", 25), rep("direct", 25), rep("inter", 25), rep("indirect", 25), rep("total", 25)),
        value = c(undiscount_cost$cost_direct, undiscount_cost$cost_inter, undiscount_cost$cost_indirect, undiscount_cost$cost_total,
            discount_cost$cost_direct,  discount_cost$cost_inter, discount_cost$cost_indirect, discount_cost$cost_total)
    )
 
    list(
        QALY = QALY,
        cost = cost,
        outcomes_age_month = outcomes_age_month
    )
 
}


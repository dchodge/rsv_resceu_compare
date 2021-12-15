#` The data used in the epidemic intervention model is as follows
#` * rsv_data_uk, demographic data from England and Wales
    # * numberDailyLiveBirths---number of daily live births in England and Wales
    # * population---total population
    # * ageGroupBoundary---age groups considered
#` data_inter_uk, data from England and Wales used in the intervention model
    #  * p_mat---proportion of each age group which are new mothers
    #  * u_p---proportion of each new mother age group
    #  * pVHR---proportion of each age groups which are considered very-high-risk
    #  * pHR---proportion of each age groups which are considered high-risk
    #  * pLR---proportion of each age groups which are considered low-risk
    #  * cnt_matrix_p---contact matrix showing mean number of total physical contacts between each age group
    #  * cnt_matrix_p_h---contact matrix showing mean number of household physical contacts between each age group
    #  * cnt_matrix_c---contact matrix showing mean number of total conversational contacts between each age group
    #  * cnt_matrix_c_h---contact matrix showing mean number of household conversational contacts between each age group
#` run_start, step number to start the model on (0)
#` run_burn, number of days to run the burn in model
#` run_full, number of days, after burn in, to run the model for


#' Initialise the RunInterventions class
#'
#' @param RunInterventions An empty RunInterventions class
#' @param resceudata Demographic data needed to define the intervention ODEs
#' @param data_inter_uk Intervention data needed to further defined the intervention ODEs
#' @param run_burn Burn-in for the model (in days)
#' @param run_full Number of days to run the model for post burn-in
#' @return A classRunInterventions which is fully parameterised
make_RunInterventions <- function(RunInterventions, resceudata, data_inter_uk, run_burn, run_full) {
    classRunInterventions <- new(RunInterventions, resceudata$numberDailyLiveBirths, resceudata$population, resceudata$ageGroupBoundary) # Calls class
    classRunInterventions$p_mat <- data_inter_uk$prop_mat
    classRunInterventions$pVHR <- data_inter_uk$pVHR
    classRunInterventions$pHR <- data_inter_uk$pHR
    classRunInterventions$pLR <- data_inter_uk$pLR

    classRunInterventions$cnt_matrix_p <- data_inter_uk$cnt_matrix_p
    classRunInterventions$cnt_matrix_p_h <- data_inter_uk$cnt_matrix_p_h
    classRunInterventions$pwp_p <- data_inter_uk$pwp_p
    classRunInterventions$pwn_p <- data_inter_uk$pwn_p
    classRunInterventions$nwp_p <- data_inter_uk$nwp_p
    classRunInterventions$nwn_p <- data_inter_uk$nwn_p

    classRunInterventions$cnt_matrix_c <- data_inter_uk$cnt_matrix_c
    classRunInterventions$cnt_matrix_c_h <- data_inter_uk$cnt_matrix_c_h
    classRunInterventions$pwp_c <- data_inter_uk$pwp_c
    classRunInterventions$pwn_c <- data_inter_uk$pwn_c
    classRunInterventions$nwp_c <- data_inter_uk$nwp_c
    classRunInterventions$nwn_c <- data_inter_uk$nwn_c

    classRunInterventions$u_p <- data_inter_uk$prop_mat[19:21] / sum(data_inter_uk$prop_mat[19:21])

    classRunInterventions$run_start <- 0.0
    classRunInterventions$run_burn <- run_burn
    classRunInterventions$run_full <- run_full
    classRunInterventions
 }


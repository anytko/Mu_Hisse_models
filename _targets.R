library(targets)
library(tarchetypes)


tar_option_set(
         packages = c("ggplot2", "corHMM", "geiger", "ape", "dplyr", "nloptr", "dentist", "phytools", "hisse")
 )

source("_functions.R")

list(
    tar_target(data, get_data_tree()),
    tar_target(geo_data, geo_rarity(data$trait)),
    tar_target(phy_data, phy_rarity(data$trait)),
    tar_target(fun_data, fun_rarity(data$trait)),
    tar_target(cr_data, cr_rarity(data$trait)),
    tar_target(geo_fun_data, geo_fun_rarity(data$trait)),
    tar_target(geo_phy_data, geo_phy_rarity(data$trait)),
    tar_target(fun_phy_data, fun_phy_rarity(data$trait)),

    tar_target(trans.rates, hisse::TransMatMakerHiSSE(hidden.traits = 0)),
    tar_target(trans.rates_2, hisse::TransMatMakerHiSSE(hidden.traits = 1)),
    tar_target(trans.rates_3, hisse::TransMatMakerHiSSE(hidden.traits = 2)),
    tar_target(trans.rates_4, hisse::TransMatMakerHiSSE(hidden.traits = 3)),
    tar_target(trans.rates_5, hisse::TransMatMakerHiSSE(hidden.traits = 4)),

    tar_target(trans.rates_mu, hisse::TransMatMakerMuHiSSE(hidden.traits = 0)),
    tar_target(trans.rates_mu_2, hisse::TransMatMakerMuHiSSE(hidden.traits = 1)),
    tar_target(trans.rates_mu_3, hisse::TransMatMakerMuHiSSE(hidden.traits = 2)),
    tar_target(trans.rates_mu_4, hisse::TransMatMakerMuHiSSE(hidden.traits = 3)),
    tar_target(trans.rates_mu_5, hisse::TransMatMakerMuHiSSE(hidden.traits = 4)),


    tar_target(geo_1_hisse, hisse::hisse(phy = data$phy, data = geo_data, trans.rate = trans.rates, turnover = c(1,2), eps = c(1,2), hidden.states = FALSE)), 
    tar_target(geo_2_hisse, hisse::hisse(phy = data$phy, data = geo_data, trans.rate = trans.rates_2, turnover = c(1, 2, 3, 4), eps = c(1,2, 3, 4), hidden.states = TRUE)),
    tar_target(geo_3_hisse, hisse::hisse(phy = data$phy, data = geo_data, trans.rate = trans.rates_3, turnover = c(1,2,3,4,5,6), eps = c(1,2,3,4,5,6), hidden.states = TRUE)), 
    tar_target(geo_4_hisse, hisse::hisse(phy = data$phy, data = geo_data, trans.rate = trans.rates_4, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)),
    tar_target(geo_5_hisse,hisse::hisse(phy = data$phy, data = geo_data, trans.rate = trans.rates_5, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), eps = c(1 ,2, 3, 4, 5, 6, 7, 8, 9, 10), hidden.states = TRUE)),

    # Functional only hisse
    tar_target(fun_1_hisse, hisse::hisse(phy = data$phy, data = fun_data, trans.rate = trans.rates, turnover = c(1,2), eps = c(1,2), hidden.states = FALSE)), 
    tar_target(fun_2_hisse, hisse::hisse(phy = data$phy, data = fun_data, trans.rate = trans.rates_2, turnover = c(1, 2, 3, 4), eps = c(1,2, 3, 4), hidden.states = TRUE)),
    tar_target(fun_3_hisse,hisse::hisse(phy = data$phy, data = fun_data, trans.rate = trans.rates_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)),
    tar_target(fun_4_hisse, hisse::hisse(phy = data$phy, data = fun_data, trans.rate = trans.rates_4, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)),
    tar_target(fun_5_hisse, hisse::hisse(phy = data$phy, data = fun_data, trans.rate = trans.rates_5, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), eps = c(1 ,2, 3, 4, 5, 6, 7, 8, 9, 10), hidden.states = TRUE)),

    # Phylogenetic only hisse
    tar_target(phy_1_hisse,hisse::hisse(phy = data$phy, data = phy_data, trans.rate = trans.rates, turnover = c(1,2), eps = c(1,2), hidden.states = FALSE)), 
    tar_target(phy_2_hisse, hisse::hisse(phy = data$phy, data = phy_data, trans.rate = trans.rates_2, turnover = c(1, 2, 3, 4), eps = c(1,2, 3, 4), hidden.states = TRUE)),
    tar_target(phy_3_hisse, hisse::hisse(phy = data$phy, data = phy_data, trans.rate = trans.rates_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)),
    tar_target(phy_4_hisse, hisse::hisse(phy = data$phy, data = phy_data, trans.rate = trans.rates_4, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)),
    tar_target(phy_5_hisse, hisse::hisse(phy = data$phy, data = phy_data, trans.rate = trans.rates_5, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), eps = c(1 ,2, 3, 4, 5, 6, 7, 8, 9, 10), hidden.states = TRUE)),


    # Classically Rare only hisse
    tar_target(cr_1_hisse, hisse::hisse(phy = data$phy, data = cr_data, trans.rate = trans.rates, turnover = c(1,2), eps = c(1,2), hidden.states = FALSE)), 
    tar_target(cr_2_hisse, hisse::hisse(phy = data$phy, data = cr_data, trans.rate = trans.rates_2, turnover = c(1, 2, 3, 4), eps = c(1,2, 3, 4), hidden.states = TRUE)),
    tar_target(cr_3_hisse, hisse::hisse(phy = data$phy, data = cr_data, trans.rate = trans.rates_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)),
    tar_target(cr_4_hisse, hisse::hisse(phy = data$phy, data = cr_data, trans.rate = trans.rates_4, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)),
    tar_target(cr_5_hisse, hisse::hisse(phy = data$phy, data = cr_data, trans.rate = trans.rates_5, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), eps = c(1 ,2, 3, 4, 5, 6, 7, 8, 9, 10), hidden.states = TRUE)),

    # Geographic and Functional MuHisse
    tar_target(geo_fun_1_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_fun_data, trans.rate = trans.rates_mu, turnover = c(1,2,3,4), eps = c(1,2,3,4), hidden.states = FALSE)),
    tar_target(geo_fun_2_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_fun_data, trans.rate = trans.rates_mu_2, turnover = c(1,2,3,4,5,6,7,8), eps = c(1,2,3,4,5,6,7,8), hidden.states = TRUE)),
    tar_target(geo_fun_3_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_fun_data, trans.rate = trans.rates_mu_3, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12), eps = c(1,2,3,4,5,6,7,8,9,10,11,12), hidden.states = TRUE)),
    tar_target(geo_fun_4_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_fun_data, trans.rate = trans.rates_mu_4, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), eps = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), hidden.states = TRUE)),
    tar_target(geo_fun_5_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_fun_data, trans.rate = trans.rates_mu_5, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), eps = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), hidden.states = TRUE)),

    # Geographic and phylogenetic MuHisse 
    tar_target(geo_phy_1_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_phy_data, trans.rate = trans.rates_mu, turnover = c(1,2,3,4), eps = c(1,2,3,4), hidden.states = FALSE)),
    tar_target(geo_phy_2_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_phy_data, trans.rate = trans.rates_mu_2, turnover = c(1,2,3,4,5,6,7,8), eps = c(1,2,3,4,5,6,7,8), hidden.states = TRUE)),
    tar_target(geo_phy_3_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_phy_data, trans.rate = trans.rates_mu_3, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12), eps = c(1,2,3,4,5,6,7,8,9,10,11,12), hidden.states = TRUE)),
    tar_target(geo_phy_4_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_phy_data, trans.rate = trans.rates_mu_4, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), eps = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), hidden.states = TRUE)),
    tar_target(geo_phy_5_muhisse, hisse::MuHiSSE(phy = data$phy, data = geo_phy_data, trans.rate = trans.rates_mu_5, turnover = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), eps = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), hidden.states = TRUE)),

    # Functional and phylogenetic MuHisse 
    tar_target(fun_phy_1_muhisse, hisse::MuHiSSE(phy = data$phy, data = fun_phy_data, trans.rate = trans.rates_mu, turnover = c(1,2,3,4), eps = c(1,2,3,4), hidden.states = FALSE)),

    tar_target(geo_hisse_summary, summarize_results(list(geo_1_hisse, geo_3_hisse, geo_5_hisse))),
    tar_target(fun_hisse_summary, summarize_results(list(fun_1_hisse, fun_3_hisse, fun_5_hisse))),
    tar_target(phy_hisse_summary, summarize_results(list(phy_1_hisse, phy_3_hisse, phy_5_hisse))),
    tar_target(cr_hisse_summary, summarize_results(list(cr_1_hisse, cr_3_hisse, cr_5_hisse))),
    tar_target(geo_fun_muhisse_summary, summarize_results(list(geo_fun_1_muhisse, geo_fun_3_muhisse, geo_fun_5_muhisse))),
    tar_target(geo_phy_muhisse_summary, summarize_results(list(geo_phy_1_muhisse, geo_phy_3_muhisse, geo_phy_5_muhisse)))
    
    )



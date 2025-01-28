# Necessary libraries
library(ape)
library(hisse)
library(dplyr)

# First - retrieve and clean the phylogeny (the pruned/cleaned phylogeny is included in files)

pruned_phylogeny <- read.tree("data/pruned_phylogeny.newick")


# Second input (rarity tip) the data

full_rarity_data <- read.csv("data/full_rarity_data.csv")

# Preparing the data

# 1 rarity type 

# Geographic rarity (geo rare vs common)
geo_data <- full_rarity_data

geo_data$geographic_rarity <- ifelse(
  geo_data$classifications %in% c("Classically Rare", "Relict", "Endemic", "Environmentally Rare"), 
  1, 
  0
)
geo_data <- geo_data %>%
  select(-classifications)

# Phylogenetic rarity (phylo rare vs common)
phy_data <- full_rarity_data

phy_data$phylogenetic_rarity <- ifelse(
  phy_data$classifications %in% c("Classically Rare", "Relict", "Indicator", "Adaptable Survivor"), 
  1, 
  0
)
phy_data <- phy_data %>%
  select(-classifications)


# Functional rarity (fun rare vs common)
fun_data <- full_rarity_data

fun_data$functional_rarity <- ifelse(
  fun_data$classifications %in% c("Classically Rare", "Endemic", "Indicator", "High Invasive Potential"), 
  1, 
  0
)
fun_data <- fun_data %>%
  select(-classifications)


# All rarity / classical rarity (classically rare vs common)
cr_data <- full_rarity_data

cr_data$geographic_rarity <- ifelse(
  cr_data$classifications %in% c("Classically Rare"), 
  1, 
  0
)
cr_data <- cr_data %>%
  select(-classifications)

# 2 rarity types 

# Geographic (first (10)) and Phylogenetic rarity (second (01))
geo_phy_data <- full_rarity_data

geo_phy_data <- geo_phy_data %>%
  mutate(
    geographic_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Relict", "Environmental Rarity"), 1, 0),
    phylogenetic_rarity = ifelse(classifications %in% c("Classically Rare", "Relict", "Indicator", "Adaptable Survivor"), 1, 0)
  )

geo_phy_data <- geo_phy_data %>%
  select(-classifications)

# Geographic (first (10)) and Functional rairty (second (01))
geo_fun_data <- full_rarity_data

geo_fun_data <- geo_fun_data %>%
  mutate(
    geographic_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Relict", "Environmental Rarity"), 1, 0),
    functional_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Indicator", "High Invasive Potential"), 1, 0)
  )

geo_fun_data <- geo_fun_data %>%
  select(-classifications)

# Functional (first (10)) and Phylogenetic rarity (second (01))
fun_phy_data <- full_rarity_data

fun_phy_data <- fun_phy_data %>%
  mutate(
    functional_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Indicator", "High Invasive Potential"), 1, 0),
    phylogenetic_rarity = ifelse(classifications %in% c("Classically Rare", "Relict", "Indicator", "Adaptable Survivor"), 1, 0)
  )

fun_phy_data <- fun_phy_data %>%
  select(-classifications)

# Run Hisse models (MuHisse next) *** All models are included as .rda objects by the same name below (i.e. load(hisse_model_geo_only.rda) -> hisse_model_geo_only$solution or hisse_model_geo_only$AIC

# No hidden states 
trans.rates <- TransMatMakerHiSSE(hidden.traits = 1)

hisse_model_geo_only <- hisse(phy = pruned_phylogeny, data = geo_data, trans.rate = trans.rates, turnover = c(1, 2), eps = c(1, 2), hidden.states = FALSE)

hisse_model_phy_only <- hisse(phy = pruned_phylogeny, data = phy_data, trans.rate = trans.rates, turnover = c(1, 2), eps = c(1, 2), hidden.states = FALSE)

hisse_model_fun_only <- hisse(phy = pruned_phylogeny, data = fun_data, trans.rate = trans.rates, turnover = c(1, 2), eps = c(1, 2), hidden.states = FALSE)

hisse_model_cr_only <- hisse(phy = pruned_phylogeny, data = cr_data, trans.rate = trans.rates, turnover = c(1, 2), eps = c(1, 2), hidden.states = FALSE)

# 2 hidden states 
trans.rates_hidden_2 <- TransMatMakerHiSSE(hidden.traits = 1)

geo_2_hisse <- hisse(phy = pruned_phylogeny, data = geo_data, trans.rate = trans.rates_hidden_2, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = TRUE)

phy_2_hisse <- hisse(phy = pruned_phylogeny, data = phy_data, trans.rate = trans.rates_hidden_2, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = TRUE)

fun_2_hisse <- hisse(phy = pruned_phylogeny, data = fun_data, trans.rate = trans.rates_hidden_2, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = TRUE)

cr_2_hisse <- hisse(phy = pruned_phylogeny, data = cr_data, trans.rate = trans.rates_hidden_2, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = TRUE)

# 3 hidden states
trans.rates_hidden_3 <- TransMatMakerHiSSE(hidden.traits = 2)

geo_3_hisse <- hisse(phy = pruned_phylogeny, data = geo_data, trans.rate = trans.rates_hidden_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)

phy_3_hisse <- hisse(phy = pruned_phylogeny, data = phy_data, trans.rate = trans.rates_hidden_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)

fun_3_hisse <- hisse(phy = pruned_phylogeny, data = fun_data, trans.rate = trans.rates_hidden_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)

cr_3_hisse <- hisse(phy = pruned_phylogeny, data = cr_data, trans.rate = trans.rates_hidden_3, turnover = c(1, 2, 3, 4, 5, 6), eps = c(1, 2, 3, 4, 5, 6), hidden.states = TRUE)

# Now the MuHisse models 

# No hidden states 
trans.rates_mu <- TransMatMakerMuHiSSE(hidden.traits=0)

muhisse_model_geo_phy <- MuHiSSE(phy = pruned_phylogeny, data = geo_phy_data, trans.rate = trans.rates_mu, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = FALSE)

muhisse_model_geo_fun <- MuHiSSE(phy = pruned_phylogeny, data = geo_fun_data, trans.rate = trans.rates_mu, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = FALSE)

muhisse_model_fun_phy <- MuHiSSE(phy = pruned_phylogeny, data = fun_phy_data, trans.rate = trans.rates_mu, turnover = c(1, 2, 3, 4), eps = c(1, 2, 3, 4), hidden.states = FALSE)

# 2 hidden states
trans.rates_mu_2 <- TransMatMakerMuHiSSE(hidden.traits=1)

geo_phy_2_muhisse <- MuHiSSE(phy = pruned_phylogeny, data = geo_phy_data, trans.rate = trans.rates_mu_2, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)

geo_fun_2_muhisse <- MuHiSSE(phy = pruned_phylogeny, data = geo_fun_data, trans.rate = trans.rates_mu_2, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)

fun_phy_2_muhisse <- MuHiSSE(phy = pruned_phylogeny, data = fun_phy_data, trans.rate = trans.rates_mu_2, turnover = c(1, 2, 3, 4, 5, 6, 7, 8), eps = c(1, 2, 3, 4, 5, 6, 7, 8), hidden.states = TRUE)

# 3 hidden states 
trans.rates_mu_3 <- TransMatMakerMuHiSSE(hidden.traits=2)

geo_phy_3_muhisse <- MuHiSSE(phy = pruned_phylogeny, data = geo_phy_data, trans.rate = trans.rates_mu_3, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), eps = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), hidden.states = TRUE)

geo_fun_3_muhisse <- MuHiSSE(phy = pruned_phylogeny, data = geo_fun_data, trans.rate = trans.rates_mu_3, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), eps = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), hidden.states = TRUE)

fun_phy_3_muhisse <- MuHiSSE(phy = pruned_phylogeny, data = fun_phy_data, trans.rate = trans.rates_mu_3, turnover = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), eps = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), hidden.states = TRUE)


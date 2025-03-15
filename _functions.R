get_data_tree <- function() {
    chrono_phylogeny <- read.tree("data/chrono_tree.newick")
    full_rarity_data <- read.csv("data/full_rarity_data.csv")
    full_rarity_data <- full_rarity_data %>% select(-X)

    return(list(trait=full_rarity_data, phy=chrono_phylogeny))
}

geo_rarity <- function(trait_data) {
  trait_data %>%
    mutate(geographic_rarity = ifelse(
      classifications %in% c("Classically Rare", "Relict", "Endemic", "Environmentally Rare"), 
      1, 
      0
    )) %>%
    select(-classifications)
}


phy_rarity <- function(trait_data) {
  trait_data %>%
    mutate(phylogenetic_rarity = ifelse(
      classifications %in% c("Classically Rare", "Relict", "Indicator", "Adaptable Survivor"), 
      1, 
      0
    )) %>%
    select(-classifications)
}


fun_rarity <- function(trait_data) {
  trait_data %>%
    mutate(functional_rarity = ifelse(
      classifications %in% c("Classically Rare", "Endemic", "Indicator", "High Invasive Potential"), 
      1, 
      0
    )) %>%
    select(-classifications)
}


cr_rarity <- function(trait_data) {
  trait_data %>%
    mutate(classical_rarity = ifelse(
      classifications %in% c("Classically Rare"), 
      1, 
      0
    )) %>%
    select(-classifications)
}



geo_fun_rarity <- function(trait_data) {
  trait_data %>%
    mutate(
      geographic_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Relict", "Environmentally Rare"), 1, 0),
      functional_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Indicator", "Potentially Invasive"), 1, 0)
    ) %>%
    select(-classifications)
}

geo_phy_rarity <- function(trait_data) {
  trait_data %>%
    mutate(
      geographic_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Relict", "Environmentally Rare"), 1, 0),
      phylogenetic_rarity = ifelse(classifications %in% c("Classically Rare", "Relict", "Indicator", "Adaptable Survivor"), 1, 0)
    ) %>%
    select(-classifications)
}

fun_phy_rarity <- function(trait_data) {
  trait_data %>%
    mutate(
      functional_rarity = ifelse(classifications %in% c("Classically Rare", "Endemic", "Indicator", "High Invasive Potential"), 1, 0),
      phylogenetic_rarity = ifelse(classifications %in% c("Classically Rare", "Relict", "Indicator", "Adaptable Survivor"), 1, 0)
    ) %>%
    select(-classifications)
}


summarize_results <- function(models) {
  # Extract AICc values from each model
  aic_values <- sapply(models, function(model) model$AICc)
  
  # Find the minimum AICc value (reference model)
  min_aic <- min(aic_values)
  
  # Calculate ΔAICc (difference from the minimum AICc)
  delta_aic <- aic_values - min_aic
  
  # Create a summary data frame
  aic_comparison <- data.frame(
    Model = names(aic_values),
    AICc = aic_values,
    Delta_AICc = delta_aic
  )
  
  return(aic_comparison)
}



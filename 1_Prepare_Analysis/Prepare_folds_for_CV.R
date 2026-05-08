# In this script we assign outer- and inner- folds to the samples, so we can use the 10-fold CV and tuning based on the same fold structure for the different models

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

# Load Data
Full <- read.csv("Data/Full_Dataset.csv",   header = TRUE, stringsAsFactors = FALSE)
Full

# Assign seed
set.seed(123)
# For 10-fold CV which accounts for the vertical autocorrelation, we need to the splits per site, not fully random, which would be a leak:
n <- max(Full$ID_Location)
folds <- sample(rep(1:10, length.out = n))

# Join assigned outer fold by IDs
Fold_match.df <- data.frame(
  ID_Location = 1:length(folds),
  Outer_Fold = folds
)

Full_Fold <- Full %>%
  select(-any_of("Outer_Fold")) %>%
  left_join(Fold_match.df, by = "ID_Location")


Main_Folds <- Full_Fold[c("ID","ID_Location","Outer_Fold")]
Main_Folds

# Store outer folds with ID
write.csv(Main_Folds, "Data/Folds/Outer/Outer_Folds.csv", row.names = FALSE)



# For each outer fold, we do a split for an inner 10-fold CV.

make_nested_folds_CV <- function(Main_Folds){

  nested_list <- list()

  for(x in 1:10){

    # --- outer split ---
    outer_split <- Main_Folds %>%
      filter(Outer_Fold == x)

    df_outer_train <- Main_Folds %>%
      filter(Outer_Fold != x)

    # --- location-level inner folds ---
    loc_ids <- sort(unique(df_outer_train$ID_Location))
    n_loc   <- length(loc_ids)

    set.seed(123 + x)
    inner_folds <- sample(rep(1:10, length.out = n_loc))

    Inner_match.df <- data.frame(
      ID_Location = loc_ids,
      Inner_Fold  = inner_folds
    )

    # assign inner folds to training data
    df_outer_train <- df_outer_train %>%
      select(-any_of("Inner_Fold")) %>%
      left_join(Inner_match.df, by = "ID_Location")

    # test set gets no inner fold
    outer_split <- outer_split %>%
      mutate(Inner_Fold = NA)

    # combine once (this is the key change)
    nested_df <- bind_rows(df_outer_train, outer_split)

    # save ONE file per outer fold
    write.csv(
      nested_df,
      paste0("Data/Folds/Inner_CV/Nested_Outer_", x, ".csv"),
      row.names = FALSE
    )

    nested_list[[paste0("Outer_", x)]] <- nested_df
  }

  return(nested_list)
}

nested_result <- make_nested_folds_CV(Main_Folds)

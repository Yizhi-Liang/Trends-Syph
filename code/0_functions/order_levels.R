order_levels <- function(df, categories) {
  # Define the order for 'year' and 'race'
  year_order <- 2014:2022
  race_order <- c(
    "AIAN",
    "Asian",
    "Black",
    "Hispanic",
    "Multiracial",
    "NHPI",
    "White"
  )
  
  # Function to apply ordering based on category
  apply_ordering <- function(df, category) {
    order_vector <- eval(rlang::sym(paste0(category, "_order")))
    df[[category]] <- factor(df[[category]], levels = order_vector)
    df
  }
  
  # Apply the ordering to each category using map
  df_ordered <- purrr::reduce(categories,
                              apply_ordering,
                              .init = dplyr::as_tibble(df))
  
  # Dynamically create the arranging instructions
  arrange_instructions <- purrr::map(categories, ~rlang::sym(.x))
  
  # Arrange the df based on the categories
  df_ordered <- df_ordered |>  dplyr::arrange(!!!arrange_instructions)
  
  return(df_ordered)
}

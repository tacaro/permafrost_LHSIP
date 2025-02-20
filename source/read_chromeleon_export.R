
#' read a chromleon export file
#' @param file_path path to the excel file to be read
#' @return a list with two data frames for injection_details and peak analysis
read_chromeleon_export <- function(file_path) {
  
  # read raw data
  raw_data <- suppressMessages(readxl::read_excel(file_path, sheet = "Integration"))
  names(raw_data) <- paste0("x", 1:ncol(raw_data))
  
  
  injection_details <- raw_data %>% 
    # remove rows that have no information
    filter(!is.na(x1)) %>% 
    # focus on rows between Injection Details and Chromatogram sections
    filter(row_number() > which(x1 == "Injection Details")[1] & row_number() < which(x1 == "Chromatogram")[1]) %>% 
    # remove columns that have no information in them at all
    { .[map_lgl(., ~!all(is.na(.x)))] }
  injection_details %>% knitr::kable()
  
  all_results <- raw_data %>% 
    filter(row_number() > which(x1 == "Integration Results")[1]) 
  
  integration_results <- 
    setNames(all_results[-c(1:3),], t(all_results[1,])[,1]) %>%
    mutate_at(vars(-`Peak Name`), function(x) suppressWarnings(as.numeric(x)))
  
  list(injection_details = injection_details, integration_results = integration_results)
}
library(dplyr)

.data <- ches_qol_data 

.data <- .data %>% 
  dplyr::select(v25_sm2, v25_task, v25_wellbeing, RAge, myopia_b, myopia_w, bvi3us_b, bvi3us_w, catsurg) %>% 
  dplyr::mutate(myopia_b = factor(myopia_b, levels = c(0, 1), labels = c("Normal", "Myopic")),
         myopia_w = factor(myopia_w, levels = c(0, 1), labels = c("Normal", "Myopic")),
         bvi3us_b = factor(bvi3us_b, levels = c(0, 1, 2), labels = c("Normal", "VI", "Blind")),
         bvi3us_w = factor(bvi3us_w, levels = c(0, 1, 2), labels = c("Normal", "VI", "Blind")),
         catsurg = factor(catsurg, levels = c(0, 1), labels = c("No Surgery", "Cataract Surgery"))) %>% 
  as.data.frame()

# TODO: Figure out why labels aren't working with dplyr. Should it use a different function?


Hmisc::label(.data$myopia_b) <- "Myopic (Better Eye)"
Hmisc::label(.data$myopia_w) <- "Myopic (Worse Eye)"
Hmisc::label(.data$bvi3us_b) <- "Visual Impairment (Better Eye)"
Hmisc::label(.data$bvi3us_w) <- "Visual Impairment (Worse Eye)"
Hmisc::label(.data$catsurg) <- "Cateract Surgery"
Hmisc::label(.data$v25_sm2) <- "Composite Score"
Hmisc::label(.data$v25_task) <- "Task Subscale"
Hmisc::label(.data$v25_wellbeing) <- "Wellbeing Subscale"

.m <- mesa(.data, -RAge, -v25_sm2, -v25_task, -v25_wellbeing) %>% 
        set_outcome(v25_sm2, v25_task, v25_wellbeing) %>% 
        column_summary() %>% 
        column_outcome() %>% 
        format_table()

.m <- mesa(.data, -myopia_b, -myopia_w) %>% 
  set_outcome(myopia_b, myopia_w) %>% 
  column_summary() %>% 
  column_outcome() %>% 
  format_table()

test <- .m %>% as.data.frame()
test %>% knitr::kable(col.names = names(test))

# missing

.m <- mesa(.data, -RAge, -v25_sm2, -v25_task, -v25_wellbeing) %>% 
  set_outcome(v25_sm2, v25_task, v25_wellbeing) %>% 
  column_summary() %>% 
  column_outcome(.missing_f = c("--", "(omitted)")) %>% 
  format_table()

library(magrittr)

missing_ttest <- function(.m, .x, .y) {
  qx <- rlang::sym(.x)
  .m$data %<>% dplyr::mutate(missing_data = is.na(rlang::UQE(qx)))
  t_test_results <- t.test(.m$data[[.y]] ~ .m$data$missing_data)
  difference <- t_test_results$estimate %>% purrr::reduce(`-`)
  p_val <- t_test_results$p.value
  data.frame(results = round_with_zeros(difference), p_val = clean_pval(p_val))
}

.m <- mesa(.data, -RAge, -v25_sm2, -v25_task, -v25_wellbeing) %>% 
  set_outcome(v25_sm2, v25_task, v25_wellbeing) %>% 
  column_summary() %>% 
  column_outcome(.missing_f = missing_ttest) %>% 
  format_table()

mesa(.data, -RAge, -myopia_b, -myopia_w, -v25_sm2) %>% 
  set_outcome(v25_sm2) %>% 
  column_summary() %>% 
  column_outcome(missing_y = "(All missing)")

#formatting 

mesa(.data, -RAge, -v25_sm2, -v25_task, -v25_wellbeing) %>% 
  set_outcome(v25_sm2, v25_task, v25_wellbeing) %>% 
  column_summary() %>% 
  column_outcome() %>% 
  format_table(indent_after = ":")

mesa(.data, -RAge, -v25_sm2, -v25_task, -v25_wellbeing) %>% 
  set_outcome(v25_sm2, v25_task, v25_wellbeing) %>% 
  column_summary() %>% 
  column_outcome() %>% 
  format_table(indent = "<div id = 'level_row'>", indent_after = "</div>")

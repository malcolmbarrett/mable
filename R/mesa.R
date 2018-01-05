mesa <- function(.data, ..., .set_id = NULL, include_missing = TRUE) {
  table_data <- .data
  if (!is.null(.set_id)) {
    .set_id <- enquo(.set_id)
    table_data <- dplyr::select(table_data, !!.set_id)
  }
  .vars <- quos(...)
  if (!is.null(.vars)) table_data <- dplyr::select(table_data, !!!.vars)
  
  #  TODO: Figure out a better approach for this
  var_labels <- Hmisc::label(table_data)
  var_names <- names(table_data)
  for (i in seq_along(var_labels)) {
    if (var_labels[i] == "") var_labels[i] <- var_names[i]
  }
 
  .m_tbl <- purrr::map2(table_data, var_labels, lay_out_var_names, include_missing = include_missing)
  
  .m <- structure(list(table = .m_tbl, table_data = table_data, data = .data, 
                       include_missing = include_missing), class = "mesa")
  .m
}

lay_out_var_names <- function(.x, .label, include_missing) {
  any_missing <- anyNA(.x)
  if (is.factor(.x)) {
    row_labels <- c(.label, levels(.x))
    if (include_missing & any_missing) row_labels <- c(row_labels, "Missing")
    return(data.frame(Variable = row_labels, stringsAsFactors = FALSE))
  } else if (is.numeric(.x)) {
    row_labels <- .label
    if (include_missing & any_missing) row_labels <- c(row_labels, "Missing")
    return(data.frame(Variable = row_labels, stringsAsFactors = FALSE))
  }
}

set_outcome <- function(.m, ...) {
  .m$outcome <- quos(...)
  .m
}

column_summary <- function(.m, .f = mean_freq, ..., cont_label = "Mean ± SD", cat_label = "Freq. (%)", include_missing = NULL) {
  if (is.null(include_missing)) include_missing <- .m$include_missing
  any_vars_cont <- any(purrr::map_lgl(.m$table_data, is.numeric))
  any_vars_cat <- any(purrr::map_lgl(.m$table_data, is.factor))
  col_label <- ifelse(any_vars_cat & any_vars_cont, paste0(cont_label, "/", cat_label),
                      ifelse(any_vars_cat & !any_vars_cont, cat_label,
                          ifelse(!any_vars_cat & any_vars_cont, cont_label, NA)))
  .m$table[] <- purrr::map(names(.m$table), function(.x) {
     .m_tbl <- .m$table[[.x]]
     
     summary_values <- .f(.m, .x, ...)
     
     .m_tbl[, length(.m_tbl) + 1] <- summary_values
     names(.m_tbl)[length(.m_tbl)] <- col_label
     .m_tbl
  })
  
  .m
}

column_outcome <- function(.m, .f = means_freqs, ..., include_missing = NULL, .missing_f = NULL) {
  if (is.null(include_missing)) include_missing <- .m$include_missing
  outcome_vars <- purrr::map_chr(.m$outcome, rlang::quo_text)
  .m$table[] <- purrr::map(names(.m$table), function(.x) {
    outcome_dfs <- purrr::map_dfc(outcome_vars, function(.y) {
        # TODO: figure out what to return when x is the outcome
        # if (.x == .y) return()
        any_missing <- anyNA(.m$data[[.x]])
        
        if (!is.null(.missing_f) & include_missing & any_missing) {
          outcome_results <- .f(.m, .x, .y, ..., include_missing = FALSE)
          if (is.function(.missing_f)){
            missing_results <- .missing_f(.m, .x, .y, ...)
          } else if (is.character(.missing_f)) {
            missing_results <- matrix(.missing_f, ncol = length(.missing_f))
            missing_results <- as.data.frame(missing_results, stringsAsFactors = FALSE)
          } else {
            stop("`.missing_f` must be either of type function or character")
          }
          
          if (ncol(outcome_results) != ncol(missing_results)) {
            stop("The results of `.missing_f` should be the same number of columns as the results of `.f`")
          }
          
          col_names <- names(outcome_results)
          names(missing_results) <- col_names
          return(rbind(outcome_results, missing_results))
          
        } else {
          
          return(.f(.m, .x, .y, ...))
        }
          
        })
    
    cbind(.m$table[[.x]], outcome_dfs)
  })
  .m
}

means_freqs <- function(.m, x, y, digits = 2, margin = 2, include_missing = NULL, 
                        include_pvalue = TRUE, missing_y = NA_character_) {
  
  if (is.null(include_missing)) include_missing <- .m$include_missing
  # categorical outcome 
  if (is.factor(.m$data[[y]])) {
    
    if (is.factor(.m$data[[x]])) {
      sum_table <- freq_perc(.m$data[c(x, y)], include_missing = include_missing, margin = margin)
      
      chi_sq_results <- chisq.test(.m$data[[x]], .m$data[[y]])
      p_val <- clean_pval(chi_sq_results$p.value)
      sum_table$p_val <- c(p_val, rep("", nrow(sum_table) - 1))
      
    } else if (is.numeric(.m$data[[x]])) {
      
      sum_table <- purrr::map_dfc(levels(.m$data[[y]]), function(.lvl) {
        qy <- rlang::sym(y)
        .filt_data <- dplyr::filter(.m$data, rlang::UQE(qy) == .lvl)
        results <- mean_sd(.filt_data[[x]], na.rm = TRUE, include_missing = FALSE)
        data.frame(results = results, stringsAsFactors = FALSE)
      }) 
      any_missing <- anyNA(.m$data[[x]])
      if (include_missing & any_missing) {
        qx <- rlang::sym(x)
        .filt_data <- dplyr::filter(.m$data, is.na(rlang::UQE(qx)))
        missing_results <- freq_perc(.filt_data[[y]], include_missing = FALSE)[-1]
        missing_results <- matrix(missing_results, ncol = length(missing_results))
        missing_results <- as.data.frame(missing_results, stringsAsFactors = FALSE)
        col_names <- names(sum_table)
        names(missing_results) <- col_names
        sum_table <- rbind(sum_table, missing_results, stringsAsFactors = FALSE)
      }
      if (length(levels(.m$data[[y]])) > 2) {
        .test_f <- oneway.test
      } else {
        .test_f <- t.test
      }
      test_results <- .test_f(.m$data[[x]] ~ .m$data[[y]])
      p_val <- clean_pval(test_results$p.value)
      sum_table$p_val <- c(p_val, rep("", nrow(sum_table) - 1))
    }
    
    names(sum_table) <- c(levels(.m$data[[y]]), "P-Value^a^")
  # continuous outcome    
  } else if (is.numeric(.m$data[[y]])) {
    if (is.factor(.m$data[[x]])) {
        x_levels <- levels(.m$data[[x]])
        results <- purrr::map_chr(x_levels, function(.lvl) {
          qx <- rlang::sym(x)
          .filt_data <- dplyr::filter(.m$data, rlang::UQE(qx) == .lvl)
          mean_sd(.filt_data[[y]], na.rm = TRUE, include_missing = FALSE)
        })
        results <- c("", results)
        # TODO add options for non-parametric tests
        p_val <- oneway.test(.m$data[, y] ~ .m$data[, x])$p.value
        p_val = c(clean_pval(p_val), rep("", length(results) - 1))
        sum_table <- data.frame(results = results, p_val = p_val, stringsAsFactors = FALSE)
        
    } else if (is.numeric(.m$data[[x]])) {
      cor_result <- cor.test(.m$data[[x]], .m$data[[y]])
      correlation <- paste("R =", round_with_zeros(cor_result$estimate, 2))
      # TODO add options for non-parametric correlation
      p_val <- clean_pval(cor_result$p.value)
      sum_table <- data.frame(results = correlation, p_val = p_val, stringsAsFactors = FALSE)
    }
    
    any_missing <- anyNA(.m$data[[x]])
    if (include_missing & any_missing) {
      qx <- rlang::sym(x)
      .filt_data <- dplyr::filter(.m$data, is.na(rlang::UQE(qx)))
      all_y_missing <- all(is.na(.filt_data[[y]]))
      
      if (all_y_missing) {
        missing_mean <- data.frame(results = missing_y, p_val = "", stringsAsFactors = FALSE)
      } else {
        missing_mean <- mean_sd(.filt_data[[y]], na.rm = TRUE, include_missing = FALSE)
        missing_mean <- data.frame(results = missing_mean, p_val = "", stringsAsFactors = FALSE)
      }
      sum_table <- rbind(sum_table, missing_mean)
    }
    
    if (!include_pvalue) sum_table <- sum_table[, -2, drop = FALSE]
    
    outcome_name <- Hmisc::label(.m$data[[y]])
    if (outcome_name == "") outcome_name <- y
    col_names <- c(paste0("Mean ", outcome_name, "/correlation^a^"), "P-Value^b^")
    if (!include_pvalue) col_names <- col_names[-2]
    names(sum_table) <- col_names
    
  }
  
  return(sum_table)
}

adjust_means <- function(.m, x, y, adjust_for, digits = 2, include_missing = NULL, include_pvalue = TRUE) {
  if (is.null(include_missing)) include_missing <- .m$include_missing
  
  # TODO: Figure out what to do about continuous vars... default, feed quartiles to effects? 
  # Then make the column summaries blank for the row since the mean will be repeated
  fmla <- as.formula(paste0(y, " ~ ",  x, "+ ", as.character(adjust_for)[2]))
  fmla_null <- as.formula(paste0(y, " ~ ", as.character(adjust_for)[2]))
  mdl <- lm(fmla, data = .m$data)
  mdl_null <- lm(fmla_null, data = .m$data, subset = !is.na(.m$data[[x]]))
  anova_results <- anova(mdl, mdl_null)
  anova_p_val <- anova_results[2, "Pr(>F)"]
  anova_p_val <- clean_pval(anova_p_val)
  efcts <- effects::effect(x, mdl)

  # TODO: think about changing the requirements for this to a list with results and pvalue, then let column_outcome() handle the blank cell issues 
  adj_table <- data.frame(results = c("", round_with_zeros(efcts$fit, digits = digits)), 
             p_value = c(anova_p_val, rep("", length(efcts$fit))), stringsAsFactors = FALSE)
  if (!include_pvalue) adj_table <- adj_table[, -2, drop = FALSE]
  
  any_missing <- anyNA(.m$data[[x]])
  if (include_missing & any_missing) adj_table[nrow(adj_table) + 1, ] <- ""
  
  # TODO: Find a better way to fix effects for zero categories
  if (nrow(adj_table) < length(.m$table[[x]][, 1])) adj_table[nrow(adj_table) + 1, ] <- ""
    
  # TODO: Change these footnotes so it's dynamic
  outcome_name <- Hmisc::label(.m$data[[y]])
  if (outcome_name == "") outcome_name <- y
  
  col_names <- c(paste0("Mean ", outcome_name, "^a^"), "P-Value^b^")
  if (!include_pvalue) col_names <- col_names[-2]
  names(adj_table) <- col_names
  adj_table
}

mean_freq <- function(.m, .var, include_missing = NULL) {
  if (is.null(include_missing)) include_missing <- .m$include_missing
  
  if (is.numeric(.m$data[[.var]])) {
    return(mean_sd(.m$data[[.var]], na.rm = TRUE, include_missing = include_missing))
  } else if (is.factor(.m$data[[.var]])) {
    return(freq_perc(.m$data[[.var]], include_missing = include_missing))
  }
}

mean_sd <- function(x, .m, digits = 2, ..., include_missing = TRUE){

  mean_x <- round_with_zeros(mean(x, ...), digits = digits)
  sd_x <- round_with_zeros(sd(x, ...), digits = digits)
  mean_sd <- paste(mean_x, "±", sd_x)
  any_missing <- anyNA(x)
  
    if (include_missing & any_missing) {
      freq_x_missing <- sum(is.na(x))
      perc_x_missing <- round_with_zeros((freq_x_missing / length(x)) * 100, digits = digits)
      mean_sd <- c(mean_sd, paste0(freq_x_missing, " (", perc_x_missing, "%)"))
    }
    return(mean_sd)
}

freq_perc <- function(x, include_missing = TRUE, digits = 2, margin = NULL) {
  if (include_missing) useNA <- "ifany" else useNA <- "no"
  freq_x <- table(x, useNA = useNA)
  if (all(freq_x == 0)) {
    perc_x <- paste(round_with_zeros(freq_x, 2))
  } else { 
    perc_x <- round_with_zeros(prop.table(freq_x, margin = margin) * 100, digits = digits)
  }
  if (is.factor(x)) {
    c("", paste0(freq_x, " (", perc_x, "%)"))
  } else if (is.data.frame(x)) {
    # TODO: add option for missing y column
   if (anyNA(colnames(freq_x))) freq_x <- freq_x[, which(!is.na(colnames(freq_x)))]
      purrr::map_dfc(colnames(freq_x), function(.col) {
          data.frame(results = c("", paste0(freq_x[, .col], " (", perc_x[, .col], "%)")), stringsAsFactors = FALSE)     
        })
    }
}

qtable <- function(.data, .vars, .outcomes) {
  # TODO
}

format_table <- function(.m, format = "markdown", indent = rep("&nbsp;", 5), 
                         indent_after = "", variable_indent = "", 
                         variable_indent_after = "", variable_bold = TRUE, 
                         italics_levels = FALSE, italics_missing = TRUE) {
  # TODO: change this to accept a custom function, send out to markdown, html, and latex functions
  .m$table[] <- purrr::map(names(.m$table), function(.x){
    .m_tbl <- .m$table[[.x]]
    
    #  style variable labels and levels
    if (variable_bold) .m_tbl[1, 1] <- paste0("**", .m_tbl[1, 1], "**")
    if (italics_levels & is.factor(.m$data[[.x]])) .m_tbl[2:nrow(.m_tbl), 1] <- paste0("*", .m_tbl[2:nrow(.m_tbl), 1], "*")
    if (italics_missing) .m_tbl[.m_tbl[, 1] == "Missing", 1] <- "*Missing*"
    
    # format before and after variable labels and levels
    if (nrow(.m_tbl) > 1) {
      indent <- paste(indent, collapse = "")
      indented <- paste0(indent, .m_tbl[2:nrow(.m_tbl), 1])
      indented <- paste0(indented, indent_after)
      .m_tbl[2:nrow(.m_tbl), 1] <- indented
      var_indented <- paste0(variable_indent, .m_tbl[1, 1])
      var_indented <- paste0(var_indented, variable_indent_after)
      .m_tbl[1, 1] <- var_indented
    }
    
    .m_tbl
  })
  
  .m
}

column_names <- function(.m, .col_names) {
  .m$column_names <- .col_names
  .m
}

row_names <- function(.m, .row_names) {
  .m$row_names <- .row_names
  .m
}

variable_names <- function(.m, .variable_names) {
  .m$variable_names <- .variable_names
  .m
}

as.data.frame.mesa <- function(.m) {
  #if (exists("variable_names", .m)) .m_tbl[, 1] <- .m$row_names
  .m_tbl <- as.data.frame(dplyr::bind_rows(.m$table))
  names(.m_tbl) <- stringr::str_replace(names(.m_tbl), "(\\^)[1-9]" , "^")
  if (exists("row_names", .m)) .m_tbl[, 1] <- .m$row_names
  if (exists("column_names", .m)) names(.m_tbl) <- .m$column_names
  .m_tbl
}

print.mesa <- function(.m) {
  print(as.data.frame(.m))
}

round_with_zeros <- function(.x, digits = 2) {
  format(round(.x, digits = digits), trim = TRUE, nsmall = digits)
}

clean_pval <- function(x, less_than = .001, digits = 3, add_stars = FALSE, 
                       stars = list(pval = c(.1, .05, .01, .001), 
                                    symbol = c(".", "*", "**", "***"))) {
  # todo: make these more dynamic above so user can control, add_stars option
  ifelse(x < less_than, paste0("<", less_than), round_with_zeros(x, digits))
}

make_formula <- function(y, x) {
  as.formula(paste(y, "~", paste(x, collapse = " + ")))
}

ci95 <- "95% CI "
ci90 <- "90% CI "
ci99 <- "99% CI "
minmax <- paste("Min.-Max.")
iqr <- "IQR "
none <- ""
parenthesis <- c("(", ")")
brackets <- c("[", "]")
dash <- "-"
comma <- ", "

est_range <- function(est, lower, upper, bound = parenthesis, divider = dash, descriptor = none, digits = 2) {
  if (bound == "") bound <- c("", "")
  paste0(round_with_zeros(est, digits = digits), " ", 
         bound[1], 
         descriptor, 
         round_with_zeros(lower, digits = digits), 
         divider, 
         round_with_zeros(upper, digits = digits), 
         bound[2])
}

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(emmeans))

############################################################################
# HELPERS                                                                  #
############################################################################

cstats.parse <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col) || is.factor(col)) {
      col_char <- as.character(col)
      col_cleaned <- sub(",", ".", col_char, fixed = TRUE)
      suppressWarnings({col_num <- as.numeric(col_cleaned)})
      if (!all(is.na(col_num)) && sum(!is.na(col_num)) > 0) {
        has_comma <- any(grepl(",", col_char, fixed = TRUE))
        if (has_comma) {
          return(col_num)
        }
      }
    }
    return(col)
  })
  return(df)
}

cstats.read <- function(path, separator = "\t") {
  df <- read.delim(path, header = T, sep = separator)
  df <- cstats.parse(df)
  attr(df, "F") <- NULL
  attr(df, "V") <- NULL
  attr(df, "R") <- NULL
  return(df)
}

cstats.set_factors <- function(df, factors) {
  missing_cols <- setdiff(factors, names(df))
  if (length(missing_cols) > 0) {stop("Factors do no match column names: ", paste(missing_cols, collapse = ", "))}
  attr(df, "F") <- factors
  return(df)
}

cstats.set_vars <- function(df, vars) {
  missing_cols <- setdiff(vars, names(df))
  if (length(missing_cols) > 0) {stop("Variables do no match column names: ", paste(missing_cols, collapse = ", "))}
  attr(df, "V") <- vars
  return(df)
}

cstats.set_reps <- function(df, replicates) {
  missing_cols <- setdiff(replicates, names(df))
  if (length(missing_cols) > 0) {stop("Variables do no match column names: ", paste(missing_cols, collapse = ", "))}
  attr(df, "R") <- replicates
  return(df)
}

cstats.aggregate <- function(df) {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  reps <- attr(df, "R")

  if (is.null(factors) || length(factors) == 0) {stop("No factors defined. Please use cstats.set_factors")}
  if (is.null(vars) || length(vars) == 0) {stop("No variables defined. Please use cstats.set_vars")}
  if (is.null(reps) || length(reps) == 0) {cat("WARNING : no replicates defined. All samples will be used as technical replicates.\n")}

  group_cols <- c(factors, reps)
  df_agg <- df %>% group_by(across(all_of(group_cols))) %>% summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  attr(df_agg, "F") <- factors
  attr(df_agg, "V") <- vars
  attr(df_agg, "R") <- reps
  
  return(as.data.frame(df_agg))
}

############################################################################
# TESTS                                                                    #
############################################################################

# 2-WAY ANOVA
cstats.2wanova <- function(df) {
  factors <- attr(df, "F")
  vars <- attr(df, "V")

  if (is.null(factors) || length(factors) == 0) {stop("No factors defined. Please use cstats.set_factors.")}
  if (is.null(vars) || length(vars) == 0) {stop("No variables defined. Please use cstats.set_vars.")}
  if (length(factors) != 2) {stop("Only one factor defined. Please use cstats.set_factors.")}

  results_list <- list()
  f1 <- factors[1]
  f2 <- factors[2]
  for (var in vars) {
    formula_str <- paste(var, "~", f1, "*", f2)
    stats_formula <- as.formula(formula_str)
    fit <- lm(stats_formula, data = df)
    anova_res <- anova(fit)
    results_list[[var]] <- anova_res
  }

  return(results_list)
}

# PAIRWISE TEST USING TUKEY TEST BASED ON 2-WAY ANOVA MODEL
cstats.tukey <- function(df, var) {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  
  if (is.null(factors) || length(factors) == 0) {stop("No factors defined. Please use cstats.set_factors.")}
  if (is.null(vars) || length(vars) == 0) {stop("No variables defined. Please use cstats.set_vars.")}
  if (length(factors) != 2) {stop("Only one factor defined. Please use cstats.set_factors.")}
  
  f1 <- factors[1]
  f2 <- factors[2]
  v <- var

  formula_anova <- as.formula(paste(v, "~", f1, "*", f2))
  model_anova <- aov(formula_anova, data = df)
  
  formula_emmeans <- as.formula(paste("~", f2, "|", f1))
  marginal_medians <- emmeans(model_anova, formula_emmeans)
  results_tukey <- pairs(marginal_medians, adjust = "tukey")
  
  return(results_tukey)
}

# PAIRWISE TEST USING MANN-WHITNEY U-TEST
cstats.wilcox <- function(df, var, base = NULL, altrn = "two.sided") {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  
  if (is.null(factors) || length(factors) == 0) {stop("No factors defined. Please use cstats.set_factors.")}
  if (is.null(vars) || length(vars) == 0) {stop("No variables defined. Please use cstats.set_vars.")}
  if (missing(var) || is.null(var)) {stop("var parameter cannot be null.")}
  if (!(var %in% vars)) {stop("Valid var parameter needed.")}
  if (length(factors) > 2) {stop("More than 2 factors defined. Please use cstats.set_factors.")}

  results <- list()

  if (length(factors) == 1) {
    facvals <- unique(df[[factors]])
    for (lvl1 in unique(df[[factors]])) {
      facvals <- facvals[facvals != lvl1]
      for (lvl2 in facvals) {
        name <- paste(lvl1, lvl2, sep = "_")
        Population_A <- df[df[[factors]] == lvl1, var]
        Population_B <- df[df[[factors]] == lvl2, var]
        if(!(length(Population_A) > 1 && length(Population_B) > 1)) {
          cat(paste("Omitting ", lvl1, " vs ", lvl2, ". Not enough observations.\n"))
        } else {
          results[[name]] <- wilcox.test(Population_A, Population_B, paired = F, alternative = altrn)
        }
      }
    }
  }
  else {
    if (missing(base) || is.null(base)) {stop("base parameter cannot be null.")}
    if (!(base %in% df[[factors[1]]])) {stop("2 factors defined. Valid base parameter needed.")}
    conda <- factors[1]
    condb <- factors[2]
    facvals <- unique(df[[condb]])
    for (lvl1 in unique(df[[condb]])) {
      facvals <- facvals[facvals != lvl1]
      for (lvl2 in facvals) {
        name <- paste(base, lvl1, lvl2, sep = "_")
        Population_A <- df[df[[conda]] == base & df[[condb]] == lvl1, var]
        Population_B <- df[df[[conda]] == base & df[[condb]] == lvl2, var]
        if(!(length(Population_A) > 1 && length(Population_B) > 1)) {
          cat(paste("Omitting ", lvl1, " vs ", lvl2, ". Not enough observations.\n"))
        } else {
          results[[name]] <- wilcox.test(Population_A, Population_B, paired = F, alternative = altrn)
        }
      }
    }
  }
  return(results)
}

# PAIRWISE TEST USING WELSCH T-TEST
cstats.t <- function(df, var, base = NULL, altrn = "two.sided") {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  
  if (is.null(factors) || length(factors) == 0) {stop("No factors defined. Please use cstats.set_factors.")}
  if (is.null(vars) || length(vars) == 0) {stop("No variables defined. Please use cstats.set_vars.")}
  if (missing(var) || is.null(var)) {stop("var parameter cannot be null.")}
  if (!(var %in% vars)) {stop("Valid var parameter needed.")}
  if (length(factors) > 2) {stop("More than 2 factors defined. Please use cstats.set_factors.")}

  results <- list()

  if (length(factors) == 1) {
    facvals <- unique(df[[factors]])
    for (lvl1 in unique(df[[factors]])) {
      facvals <- facvals[facvals != lvl1]
      for (lvl2 in facvals) {
        name <- paste(lvl1, lvl2, sep = "_")
        Population_A <- df[df[[factors]] == lvl1, var]
        Population_B <- df[df[[factors]] == lvl2, var]
        if(!(length(Population_A) > 1 && length(Population_B) > 1)) {
          cat(paste("Omitting ", lvl1, " vs ", lvl2, ". Not enough observations.\n"))
        } else {
          results[[name]] <- t.test(Population_A, Population_B, paired = F, alternative = altrn)
        }
      }
    }
  }
  else {
    if (missing(base) || is.null(base)) {stop("base parameter cannot be null.")}
    if (!(base %in% df[[factors[1]]])) {stop("2 factors defined. Valid base parameter needed.")}
    conda <- factors[1]
    condb <- factors[2]
    facvals <- unique(df[[condb]])
    for (lvl1 in unique(df[[condb]])) {
      facvals <- facvals[facvals != lvl1]
      for (lvl2 in facvals) {
        name <- paste(base, lvl1, lvl2, sep = "_")
        Population_A <- df[df[[conda]] == base & df[[condb]] == lvl1, var]
        Population_B <- df[df[[conda]] == base & df[[condb]] == lvl2, var]
        if(!(length(Population_A) > 1 && length(Population_B) > 1)) {
          cat(paste("Omitting ", lvl1, " vs ", lvl2, ". Not enough observations.\n"))
        } else {
          results[[name]] <- t.test(Population_A, Population_B, paired = F, alternative = altrn)
        }
      }
    }
  }

  return(results)
}


############################################################################
# PLOTS                                                                    #
############################################################################

# Plot barplots using Tukey test based on 2-way ANOVA model
cstats.barplot.tukey <- function(df, var) {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  
  if (is.null(factors) || length(factors) != 2) {
    stop("Es requereixen exactament 2 factors definits. Utilitza cstats.set_factors.")
  }
  if (!(var %in% vars)) {
    stop("La variable especificada no existeix dins de les variables de 'V'.")
  }
  
  f1 <- factors[1]  # Primer factor (ex: stage)
  f2 <- factors[2]  # Segon factor (ex: condition)
  v <- var
  
  # 1. Ajustar el model ANOVA de 2 vies amb interacció per a tota la dada
  formula_anova <- as.formula(paste(v, "~", f1, "*", f2))
  model_anova <- aov(formula_anova, data = df)
  
  # 2. Obtenir les comparacions de Tukey de f2 per a cada nivell de f1 usant emmeans
  formula_emmeans <- as.formula(paste("~", f2, "|", f1))
  marginal_medians <- emmeans(model_anova, formula_emmeans)
  results_tukey <- as.data.frame(pairs(marginal_medians, adjust = "tukey"))
  
  # Funció auxiliar per convertir p-valor en asterisques
  p_to_stars <- function(p) {
    if (is.na(p)) return("ns")
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  
  stages <- unique(df[[f1]])
  plot_list <- list()
  
  for (stg in stages) {
    # 1. Filtrar dades per a l'etapa/nivell actual de f1
    df_stage <- df %>% filter(.data[[f1]] == stg)
    
    # 2. Calcular mitjana i SD per a cada nivell de f2 en aquest f1
    summary_df <- df_stage %>%
      group_by(across(all_of(f2))) %>%
      summarise(
        mean_val = mean(.data[[v]], na.rm = TRUE),
        sd_val = sd(.data[[v]], na.rm = TRUE),
        .groups = "drop"
      )
    
    # Obtenir l'ordre de les condicions (f2) a l'eix X
    conditions <- unique(as.character(summary_df[[f2]]))
    
    # 3. Filtrar els resultats de Tukey per al nivell actual de f1
    tukey_stage <- results_tukey %>% filter(.data[[f1]] == stg)
    
    sig_annotations <- data.frame(
      c1 = character(), c2 = character(), label = character(), 
      x1 = numeric(), x2 = numeric(), stringsAsFactors = FALSE
    )
    
    if (nrow(tukey_stage) > 0) {
      for (i in 1:nrow(tukey_stage)) {
        row <- tukey_stage[i, ]
        p_val <- row$p.value
        stars <- p_to_stars(p_val)
        
        if (stars != "ns") {
          # Extreure les dues condicions comparades del contrast d'emmeans (ex: "A - B")
          contrast_str <- as.character(row$contrast)
          parts <- trimws(unlist(strsplit(contrast_str, "-")))
          
          if (length(parts) == 2) {
            lvl1 <- parts[1]
            lvl2 <- parts[2]
            
            x1_pos <- match(lvl1, conditions)
            x2_pos <- match(lvl2, conditions)
            
            if (!is.na(x1_pos) && !is.na(x2_pos)) {
              sig_annotations <- rbind(sig_annotations, data.frame(
                c1 = lvl1, c2 = lvl2, label = stars, 
                x1 = x1_pos, x2 = x2_pos, stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
    
    # 4. Calcular límits de l'eix Y i altures apilades per als claudàtors
    max_y <- max(df_stage[[v]], na.rm = TRUE)
    max_sum <- max(summary_df$mean_val + summary_df$sd_val, na.rm = TRUE)
    base_y <- max(max_y, max_sum, na.rm = TRUE)
    
    step_y <- base_y * 0.12 
    tick_h <- base_y * 0.03
    
    if (nrow(sig_annotations) > 0) {
      sig_annotations$y_pos <- base_y + step_y * (1:nrow(sig_annotations))
      upper_limit <- base_y + step_y * (nrow(sig_annotations) + 1.2)
    } else {
      upper_limit <- base_y * 1.15
    }
    
    # 5. Construir el barplot amb ggplot2
    p <- ggplot(summary_df, aes(x = .data[[f2]], y = mean_val, fill = .data[[f2]])) +
      geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
      geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), 
                    width = 0.2, position = position_dodge(0.6)) +
      geom_point(data = df_stage, aes(x = .data[[f2]], y = .data[[v]]), 
                 position = position_jitter(width = 0.1), size = 1.8, alpha = 0.8) +
      scale_y_continuous(limits = c(0, upper_limit)) +
      theme_minimal() +
      labs(
        title = paste(f1, ":", stg),
        subtitle = paste("Variable:", v),
        x = f2,
        y = v,
        fill = f2
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "none"
      )
    
    # 6. Afegir els claudàtors i estrelles de significança
    if (nrow(sig_annotations) > 0) {
      p <- p + 
        geom_segment(data = sig_annotations, aes(x = x1, xend = x2, y = y_pos, yend = y_pos), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        geom_segment(data = sig_annotations, aes(x = x1, xend = x1, y = y_pos, yend = y_pos - tick_h), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        geom_segment(data = sig_annotations, aes(x = x2, xend = x2, y = y_pos, yend = y_pos - tick_h), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        geom_text(data = sig_annotations, aes(x = (x1 + x2) / 2, y = y_pos, label = label), 
                  inherit.aes = FALSE, vjust = -0.4, size = 4.5, fontface = "bold")
    }
    
    plot_list[[as.character(stg)]] <- p
  }
  
  return(plot_list)
}


# Plot barplots using wilcox test
cstats.barplot.wilcox <- function(df, var) {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  
  if (is.null(factors) || length(factors) != 2) {
    stop("Es requereixen exactament 2 factors definits (ex: stage i condition).")
  }
  if (!(var %in% vars)) {
    stop("La variable especificada no existeix dins de les variables de 'V'.")
  }
  
  f_stage <- factors[1]  # Primer factor (ex: stage)
  f_cond  <- factors[2]  # Segon factor (ex: condition)
  
  stages <- unique(df[[f_stage]])
  plot_list <- list()
  
  # Funció auxiliar per convertir p-valor en asterisques
  p_to_stars <- function(p) {
    if (is.na(p)) return("ns")
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  
  for (stg in stages) {
    # 1. Filtrar dades per a l'etapa actual
    df_stage <- df %>% filter(.data[[f_stage]] == stg)
    
    # 2. Calcular mitjana i SD per a cada condició
    summary_df <- df_stage %>%
      group_by(across(all_of(f_cond))) %>%
      summarise(
        mean_val = mean(.data[[var]], na.rm = TRUE),
        sd_val = sd(.data[[var]], na.rm = TRUE),
        .groups = "drop"
      )
    
    # Obtenir l'ordre de les condicions a l'eix X
    conditions <- unique(as.character(summary_df[[f_cond]]))
    
    # 3. Executar test de Wilcoxon per parelles i preparar anotacions
    sig_annotations <- data.frame(
      c1 = character(), c2 = character(), label = character(), 
      x1 = numeric(), x2 = numeric(), stringsAsFactors = FALSE
    )
    
    if (length(conditions) >= 2) {
      for (i in 1:(length(conditions) - 1)) {
        for (j in (i + 1):length(conditions)) {
          lvl1 <- conditions[i]
          lvl2 <- conditions[j]
          
          popA <- df_stage %>% filter(.data[[f_cond]] == lvl1) %>% pull(.data[[var]])
          popB <- df_stage %>% filter(.data[[f_cond]] == lvl2) %>% pull(.data[[var]])
          
          if (length(popA) > 1 && length(popB) > 1) {
            test_res <- wilcox.test(popA, popB, paired = FALSE)
            stars <- p_to_stars(test_res$p.value)
            
            if (stars != "ns") {
              x1_pos <- match(lvl1, conditions)
              x2_pos <- match(lvl2, conditions)
              
              sig_annotations <- rbind(sig_annotations, data.frame(
                c1 = lvl1, c2 = lvl2, label = stars, 
                x1 = x1_pos, x2 = x2_pos, stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
    
    # 4. Calcular límits de l'eix Y i altures apilades per als claudàtors
    max_y <- max(df_stage[[var]], na.rm = TRUE)
    max_sum <- max(summary_df$mean_val + summary_df$sd_val, na.rm = TRUE)
    base_y <- max(max_y, max_sum, na.rm = TRUE)
    
    step_y <- base_y * 0.12 
    tick_h <- base_y * 0.03
    
    if (nrow(sig_annotations) > 0) {
      # Assignar una alçada (y) incremental per a cada claudàtor de manera vectoritzada
      sig_annotations$y_pos <- base_y + step_y * (1:nrow(sig_annotations))
      upper_limit <- base_y + step_y * (nrow(sig_annotations) + 1.2)
    } else {
      upper_limit <- base_y * 1.15
    }
    
    # 5. Construir el barplot amb ggplot2
    p <- ggplot(summary_df, aes(x = .data[[f_cond]], y = mean_val, fill = .data[[f_cond]])) +
      geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
      geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), 
                    width = 0.2, position = position_dodge(0.6)) +
      geom_point(data = df_stage, aes(x = .data[[f_cond]], y = .data[[var]]), 
                 position = position_jitter(width = 0.1), size = 1.8, alpha = 0.8) +
      scale_y_continuous(limits = c(0, upper_limit)) +
      theme_minimal() +
      labs(
        title = paste("Etapa:", stg),
        subtitle = paste("Variable:", var),
        x = "Condició",
        y = var,
        fill = "Condició"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "none"
      )
    
    # 6. Afegir els claudàtors i estrelles de manera VECTORITZADA (sense bucle for)
    if (nrow(sig_annotations) > 0) {
      p <- p + 
        # Línia horitzontal del claudàtor
        geom_segment(data = sig_annotations, aes(x = x1, xend = x2, y = y_pos, yend = y_pos), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        # Pota esquerra
        geom_segment(data = sig_annotations, aes(x = x1, xend = x1, y = y_pos, yend = y_pos - tick_h), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        # Pota dreta
        geom_segment(data = sig_annotations, aes(x = x2, xend = x2, y = y_pos, yend = y_pos - tick_h), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        # Text de les estrelles al damunt
        geom_text(data = sig_annotations, aes(x = (x1 + x2) / 2, y = y_pos, label = label), 
                  inherit.aes = FALSE, vjust = -0.4, size = 4.5, fontface = "bold")
    }
    
    # Guardar a la llista
    plot_list[[as.character(stg)]] <- p
  }
  
  return(plot_list)
}

# Plot barplots using t.test
cstats.barplot.t <- function(df, var) {
  factors <- attr(df, "F")
  vars <- attr(df, "V")
  
  if (is.null(factors) || length(factors) != 2) {
    stop("Es requereixen exactament 2 factors definits (ex: stage i condition).")
  }
  if (!(var %in% vars)) {
    stop("La variable especificada no existeix dins de les variables de 'V'.")
  }
  
  f_stage <- factors[1]  # Primer factor (ex: stage)
  f_cond  <- factors[2]  # Segon factor (ex: condition)
  
  stages <- unique(df[[f_stage]])
  plot_list <- list()
  
  # Funció auxiliar per convertir p-valor en asterisques
  p_to_stars <- function(p) {
    if (is.na(p)) return("ns")
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else return("ns")
  }
  
  for (stg in stages) {
    # 1. Filtrar dades per a l'etapa actual
    df_stage <- df %>% filter(.data[[f_stage]] == stg)
    
    # 2. Calcular mitjana i SD per a cada condició
    summary_df <- df_stage %>%
      group_by(across(all_of(f_cond))) %>%
      summarise(
        mean_val = mean(.data[[var]], na.rm = TRUE),
        sd_val = sd(.data[[var]], na.rm = TRUE),
        .groups = "drop"
      )
    
    # Obtenir l'ordre de les condicions a l'eix X
    conditions <- unique(as.character(summary_df[[f_cond]]))
    
    # 3. Executar test de Wilcoxon per parelles i preparar anotacions
    sig_annotations <- data.frame(
      c1 = character(), c2 = character(), label = character(), 
      x1 = numeric(), x2 = numeric(), stringsAsFactors = FALSE
    )
    
    if (length(conditions) >= 2) {
      for (i in 1:(length(conditions) - 1)) {
        for (j in (i + 1):length(conditions)) {
          lvl1 <- conditions[i]
          lvl2 <- conditions[j]
          
          popA <- df_stage %>% filter(.data[[f_cond]] == lvl1) %>% pull(.data[[var]])
          popB <- df_stage %>% filter(.data[[f_cond]] == lvl2) %>% pull(.data[[var]])
          
          if (length(popA) > 1 && length(popB) > 1) {
            test_res <- t.test(popA, popB, paired = FALSE)
            stars <- p_to_stars(test_res$p.value)
            
            if (stars != "ns") {
              x1_pos <- match(lvl1, conditions)
              x2_pos <- match(lvl2, conditions)
              
              sig_annotations <- rbind(sig_annotations, data.frame(
                c1 = lvl1, c2 = lvl2, label = stars, 
                x1 = x1_pos, x2 = x2_pos, stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
    
    # 4. Calcular límits de l'eix Y i altures apilades per als claudàtors
    max_y <- max(df_stage[[var]], na.rm = TRUE)
    max_sum <- max(summary_df$mean_val + summary_df$sd_val, na.rm = TRUE)
    base_y <- max(max_y, max_sum, na.rm = TRUE)
    
    step_y <- base_y * 0.12 
    tick_h <- base_y * 0.03
    
    if (nrow(sig_annotations) > 0) {
      # Assignar una alçada (y) incremental per a cada claudàtor de manera vectoritzada
      sig_annotations$y_pos <- base_y + step_y * (1:nrow(sig_annotations))
      upper_limit <- base_y + step_y * (nrow(sig_annotations) + 1.2)
    } else {
      upper_limit <- base_y * 1.15
    }
    
    # 5. Construir el barplot amb ggplot2
    p <- ggplot(summary_df, aes(x = .data[[f_cond]], y = mean_val, fill = .data[[f_cond]])) +
      geom_bar(stat = "identity", alpha = 0.7, width = 0.6) +
      geom_errorbar(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), 
                    width = 0.2, position = position_dodge(0.6)) +
      geom_point(data = df_stage, aes(x = .data[[f_cond]], y = .data[[var]]), 
                 position = position_jitter(width = 0.1), size = 1.8, alpha = 0.8) +
      scale_y_continuous(limits = c(0, upper_limit)) +
      theme_minimal() +
      labs(
        title = paste("Etapa:", stg),
        subtitle = paste("Variable:", var),
        x = "Condició",
        y = var,
        fill = "Condició"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "none"
      )
    
    # 6. Afegir els claudàtors i estrelles de manera VECTORITZADA (sense bucle for)
    if (nrow(sig_annotations) > 0) {
      p <- p + 
        # Línia horitzontal del claudàtor
        geom_segment(data = sig_annotations, aes(x = x1, xend = x2, y = y_pos, yend = y_pos), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        # Pota esquerra
        geom_segment(data = sig_annotations, aes(x = x1, xend = x1, y = y_pos, yend = y_pos - tick_h), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        # Pota dreta
        geom_segment(data = sig_annotations, aes(x = x2, xend = x2, y = y_pos, yend = y_pos - tick_h), 
                     inherit.aes = FALSE, linewidth = 0.4) +
        # Text de les estrelles al damunt
        geom_text(data = sig_annotations, aes(x = (x1 + x2) / 2, y = y_pos, label = label), 
                  inherit.aes = FALSE, vjust = -0.4, size = 4.5, fontface = "bold")
    }
    
    # Guardar a la llista
    plot_list[[as.character(stg)]] <- p
  }
  
  return(plot_list)
}

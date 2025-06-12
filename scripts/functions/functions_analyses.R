#### Title: Functions for data analyses 
#### Author: María Bueno Álvez
#### Description: script collecting functions to perform data analyses (DE & ML)
#### Last edited : 12/08/2024

# Data analyses packages
#library(limma)
#library(tidymodels)
#library(recipes)
#library(themis)
#library(vip)
#library(pROC)    
#library(multiROC)
#library(rsample)
#library(future)
#library(furrr)
#library(MatchIt)
#library(pROC)
#library(effectsize)

# Function to run ANOVA for a given protein
do_anova <- function(df, protein) {
  
  df <- 
    df %>%
    filter(Assay == protein)
  
  # Linear model with all variables
  model <- lm(NPX ~ Disease + Sex + Age + BMI, data = df)
  
  # Conduct ANOVA
  anova_res <- anova(model)
  
  # Calculate Eta-squared using effectsize package
  eta_squared_res <- 
    eta_squared(model, partial = TRUE) |> 
    as.data.frame()
  
  # Get tidy model results
  tidy_anova <- broom::tidy(anova_res)
  
  # Add Eta-squared to the result
  tidy_anova <- 
    tidy_anova |> 
    left_join(eta_squared_res, by = c("term" = "Parameter")) |> 
    mutate(Eta2_perc = Eta2_partial * 100) |> 
    mutate(Protein = protein) |> 
    relocate(Protein)
  
  return(tidy_anova)
}

# Function to match case samples to controls
match_samples <-  function(metadata, 
                           case, 
                           control) {
  
  if(case %in% female_diseases) {
    
    metadata <- 
      metadata |> 
      filter(Sex == "F")
    
  } else if (case %in% male_diseases) {
    
    metadata <- 
      metadata |> 
      filter(Sex == "M")
    
  } else {
    metadata <- metadata
  }
  
  dat <-
    metadata |>
    filter(Disease %in% c(case, control),
           !is.na(Sex),
           !is.na(Age)) |> 
    mutate(Disease = factor(Disease, levels = c(control, case)))
  
  if(case %in% c(male_diseases, female_diseases)) {
    
    set.seed(213)
    output <- matchit(
      Disease ~ Age,
      data = dat,
      method = "nearest",
      distance = "glm"
    )
    
  } else {
    set.seed(213)
    output <- matchit(
      Disease ~ Age + Sex,
      data = dat,
      method = "nearest",
      distance = "glm"
    )
  }
  
  
  match.data(output)
  
}

# Function to run differential expression using limma
do_limma_disease <-
  function(data_wide, 
           metadata,
           disease,
           controls,
           correct = T,
           cutoff = 0) {
    
    # Select current disease
    dat <-
      data_wide %>% 
      inner_join(metadata %>% 
                   select(DAid, Sex, Age, BMI, Disease), by = "DAid") %>% 
      rename(Group = Disease) %>% 
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control")) 
    
    
    if(correct == T) {
     dat <- 
       dat |> 
       filter(!is.na(Sex),
              !is.na(Age))
    } else {
     
      dat <- dat
    }
    
    # Design a model
    if(correct == T) {
      
      if(disease %in% c(male_diseases, female_diseases)) {
        design <- model.matrix(~0 + as.factor(dat$Group) + dat$Age) 
        colnames(design) <- c("control", "case", "Age")
      } else if(disease %in% pediatric_diseases & controls == "Healthy") {
        design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex)) 
        colnames(design) <- c("control", "case",  "Sex") 
      } else {
        design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age) 
        colnames(design) <- c("control", "case",  "Sex", "Age") 
      }
      
    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)
    
    # Fit linear model to each protein assay
    dat_fit <- 
      dat %>% 
      select(-Sex, -Age, -BMI, -Group)  %>% 
      column_to_rownames("DAid") %>% 
      t()
    
    fit <- lmFit(dat_fit, design = design, maxit = 100000) #method = "robust", 
    
    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)
    
    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit, robust = T)
    
    # Extract DE results
    DE_results <-
      topTable(ebays_fit,
               n = nrow(ebays_fit$p.value), 
               adjust.method = "fdr", 
               confint = TRUE)
    
    DE_res <- 
      DE_results %>% 
      as_tibble(rownames = "Assay") %>% 
      mutate(Disease = disease,
             sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
                             adj.P.Val < 0.05 & logFC > cutoff ~ "significant up", 
                             T ~ "not significant"),
             Control = controls)
    
    return(DE_res)
  }

# Function
do_limma_discrete <- function(join_data,
                              variable = "Disease",
                              case,
                              control,
                              correct = c("Sex", "Age"),
                              correct_type = c("factor", "numeric"),
                              only_female = NULL,
                              only_male = NULL,
                              pval_lim = 0.05,
                              logfc_lim = 0) {
  Variable <- rlang::sym(variable)
  
  if (variable == "Disease") {
    sex_specific <- FALSE
    # Filter for Sex if disease is Sex specific
    if(!is.null(only_female) & case %in% only_female) {
      join_data <- join_data |>
        dplyr::filter(Sex == "F")
      sex_specific <- TRUE
    } else {
      join_data <- join_data
    }
    
    if(!is.null(only_male) & case %in% only_male) {
      join_data <- join_data |>
        dplyr::filter(Sex == "M")
      sex_specific <- TRUE
    } else {
      join_data <- join_data
    }
  }
  
  nrows_before <- nrow(join_data)
  
  join_data <- join_data |>
    dplyr::filter(!dplyr::if_any(dplyr::all_of(c(variable, correct)), is.na)) |>  # Remove NAs from columns in formula
    dplyr::filter(!!Variable %in% c(case, control)) |>
    dplyr::mutate(!!Variable := ifelse(!!Variable == case, "1_Case", "0_Control"))
  
  nrows_after <- nrow(join_data)
  if (nrows_before != nrows_after){
    warning(paste0(nrows_before - nrows_after,
                   " rows were removed because they contain NAs in ",
                   variable,
                   " or ",
                   paste(correct, collapse = ", "),
                   "!"))
  }
  
  # Design a model
  formula <- paste("~0 + as.factor(", variable, ")")
  
  if (!is.null(correct)) {
    for (i in 1:length(correct)) {
      if (correct_type[i] == "factor") {
        if (correct[i] == "Sex" && sex_specific == TRUE) {
          join_data <- join_data |>
            dplyr::select(-Sex)
          next
        } else {
          cofactor = paste("as.factor(", correct[i], ")")
        }
      } else {
        cofactor = correct[i]
      }
      formula <- paste(formula, "+", cofactor)
    }
  }
  
  if (c("Sex") %in% correct && sex_specific == TRUE) {
    correct <- correct[!correct == 'Sex']
  }
  
  design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
  cols <- c("control", "case", correct)
  cols <- cols[!is.null(cols)]
  colnames(design) <- paste(cols)
  contrast <- limma::makeContrasts(Diff = case - control, levels = design)
  
  # Fit linear model to each protein assay
  data_fit <- join_data |>
    dplyr::select(-dplyr::any_of(c(variable, "Sex", correct))) |>
    tibble::column_to_rownames("DAid") |>
    t()
  
  fit <- limma::lmFit(data_fit, design = design, method = "robust", maxit = 10000)
  contrast_fit <- limma::contrasts.fit(fit, contrast)
  
  # Apply empirical Bayes smoothing to the SE
  ebays_fit <- limma::eBayes(contrast_fit)
  
  # Extract DE results
  de_results <- limma::topTable(ebays_fit,
                                n = nrow(ebays_fit$p.value),
                                adjust.method = "fdr",
                                confint = TRUE)
  
  de_res <- de_results |>
    tibble::as_tibble(rownames = "Assay") |>
    dplyr::mutate(!!Variable := case) |>
    dplyr::mutate(sig = dplyr::case_when(
      adj.P.Val < pval_lim & logFC < -logfc_lim ~ "significant down",
      adj.P.Val < pval_lim & logFC > logfc_lim ~ "significant up",
      T ~ "not significant")
    ) |>
    dplyr::arrange(adj.P.Val)
  
  return(de_res)
}


do_limma_continuous <- function(join_data,
                                   variable,
                                   correct = c("Sex"),
                                   correct_type = c("factor"),
                                   pval_lim = 0.05,
                                   logfc_lim = 0) {
  nrows_before <- nrow(join_data)
  
  join_data <- join_data |>
    dplyr::filter(!dplyr::if_any(dplyr::all_of(c(variable, correct)), is.na))  # Remove NAs from columns in formula
  
  nrows_after <- nrow(join_data)
  if (nrows_before != nrows_after){
    warning(paste0(nrows_before - nrows_after,
                   " rows were removed because they contain NAs in ",
                   variable,
                   " or ",
                   paste(correct, collapse = ", "),
                   "!"))
  }
  
  # Design a model
  # formula <- paste("~0 +" , variable)
  
  # if (!is.null(correct)) {
  #   for (i in 1:length(correct)) {
  #     if (correct_type[i] == "factor") {
  #       cofactor = paste("as.factor(", correct[i], ")")
  #     } else {
  #       cofactor = correct[i]
  #     }
  #     formula <- paste(formula, "+", cofactor)
  #   }
  # }
  #design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
  
  # Design a model
  formula <- paste("~" , variable)  # Include intercept
  
  if (!is.null(correct)) {
    for (i in 1:length(correct)) {
      if (correct_type[i] == "factor") {
        cofactor = paste("as.factor(", correct[i], ")")
      } else {
        cofactor = correct[i]
      }
      formula <- paste(formula, "+", cofactor)
    }
  }
  design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
  
  
  # Fit linear model to each protein assay
  data_fit <- 
    join_data |>
    dplyr::select(-dplyr::any_of(c(variable, correct))) |>
    tibble::column_to_rownames("DAid") |>
    t()
  
  fit <- limma::lmFit(data_fit, design = design, method = "robust", maxit = 10000)
  
  
  # Apply empirical Bayes smoothing to the SE
  ebays_fit <- limma::eBayes(fit)
  
  # Extract DE results
  de_results <- limma::topTable(ebays_fit,
                                coef = variable,  # Extract results for Age specifically
                                n = nrow(ebays_fit$p.value),
                                adjust.method = "fdr",
                                confint = TRUE)
  
  de_res <- de_results |>
    tibble::as_tibble(rownames = "Assay") |>
    dplyr::rename(logFC = colnames(de_results)[1]) |>
    dplyr::mutate(sig = dplyr::case_when(
      adj.P.Val < pval_lim & logFC < -logfc_lim ~ "significant down",
      adj.P.Val < pval_lim & logFC > logfc_lim ~ "significant up",
      T ~ "not significant")
    ) |>
    dplyr::arrange(adj.P.Val)
  
  return(de_res)
}

# Function to generate data splits for ML
generate_split <- function(data, 
                           proportion = 0.7,
                           seed = 213,
                           variable_stratify) {
  
  set.seed(seed)
  
  data_split <-
    data |> 
    initial_split(prop = proportion, strata = !!sym(variable_stratify))
  
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  return(list("data_split" = data_split, 
              "data_train" = data_train,
              "data_test" = data_test))
  
}

# Function to run logistic regression
do_log_reg <- function(protein, 
                       disease_case,
                       disease_control,
                       disease_train,
                       disease_test) {
  
  # Filter and preprocess data for the current protein
  train_data <- 
    disease_train |> 
    filter(Disease %in% c(disease_case, disease_control)) |> 
    mutate(Disease = ifelse(Disease %in% disease_control, "Control", "Case")) |>
    select(DAid, Disease, protein) |> 
    mutate(Disease = factor(Disease))
  
  test_data <- 
    disease_test |> 
    filter(Disease %in% c(disease_case, disease_control)) |> 
    mutate(Disease = ifelse(Disease %in% disease_control, "Control", "Case")) |>
    select(DAid, Disease, protein) |> 
    mutate(Disease = factor(Disease))
  
  # Define the recipe and model specification
  recipe <- 
    recipe(Disease ~ ., data = train_data) %>%
    update_role(DAid, new_role = "id") # Update role for the sample ID
  
  model_spec <- 
    logistic_reg() %>%
    set_engine("glm")
  
  # Create the workflow
  workflow <-
    workflow() %>%
    add_recipe(recipe) %>%
    add_model(model_spec)
  
  # Fit the model on the training data
  set.seed(213)
  fit <- 
    workflow %>%
    fit(data = train_data)
  
  # Predict on the test data
  predictions <- 
    predict(fit, test_data, type = "prob") %>%
    bind_cols(test_data |> 
                select(DAid, Disease))
  
  # Evaluate the model using ROC AUC
  performance <- roc_auc(predictions, truth = Disease, `.pred_Case`)
  
  roc_curve <- roc_curve(predictions, truth = Disease, `.pred_Case`) 
  
  predictions_discrete <- 
    predict(fit, test_data) |> 
    bind_cols(test_data |> 
                select(DAid, Disease))
  
  cm <- conf_mat(predictions_discrete, truth = Disease, `.pred_class`) #|> autoplot(type = "heatmap")
  
  
  return(list(AUC = performance$.estimate,
              predictions = predictions,
              roc_curve = roc_curve,
              cm = cm))
}
  
# Function to perform over-representation analysis
do_ORA <-
  function(gene_associations,
           database,
           universe,
           n_clusters = 5,
           minGSSize = 10,
           maxGSSize = Inf,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) {
    require(clusterProfiler)
    require(multidplyr)
    
    if(n_clusters != 1) {
      worker_cluster <- new_cluster(n = n_clusters)
      cluster_library(worker_cluster, c("dplyr",
                                        "tidyverse"))
      cluster_copy(worker_cluster, c("enricher",
                                     "universe",
                                     "database",
                                     "minGSSize",
                                     "maxGSSize",
                                     "pvalueCutoff",
                                     "qvalueCutoff" ))
      
      pre_out <- 
        gene_associations %>%
        group_by(partition) %>%
        partition(worker_cluster) 
    } else {
      pre_out <- 
        gene_associations %>%
        group_by(partition)
    }
    
    outdata <-
      pre_out %>% 
      do({
        g_data <- .
        pull(g_data, gene) %>%
          enricher(universe = universe,
                   TERM2GENE = database, 
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   qvalueCutoff = qvalueCutoff) %>%
          as_tibble()
      }) %>%
      ungroup() %>%
      collect()
    
    if(n_clusters != 1) rm(worker_cluster)
    outdata
  }


extract_roc_multiclass <- function(predictions) {
  roc_data <- 
    predictions |> 
    filter(seed == 1)
  
  dat_true <-
    roc_data |> 
    select(DAid, Disease) |>  
    mutate(value = 1) |> 
    spread(Disease,value, fill= 0) |> 
    column_to_rownames("DAid") 
  
  dat_true_final <- 
    dat_true |> 
    set_names(paste(names(dat_true), "_true", sep = "")) |> 
    rownames_to_column("DAid")
  
  dat_prob <- 
    roc_data |> 
    select(DAid, starts_with(".pred")) |> 
    rename_all(~stringr::str_replace_all(.,".pred_","")) |> 
    column_to_rownames("DAid") 
  
  dat_prob_final <- 
    dat_prob |>  
    set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) |> 
    rownames_to_column("DAid")
  
  final_df <- 
    dat_true_final %>% 
    left_join(dat_prob_final, by = "DAid") %>% 
    select(-DAid)
  
  roc_res <- multi_roc(final_df, force_diag=T)
  plot_roc_df <- plot_roc_data(roc_res)
  
  return(plot_roc_df)
}


log_reg_multiclass <- function(protein, disease) {
  
  res_lr_classes <-
    
    map(c(classes, "Healthy"), function(class) {
      control_diseases <-
        resource_meta |>
        filter(Class == class) |>
        distinct(Disease) |>
        filter(Disease != disease) |>
        pull()
      
      do_log_reg (
        protein = protein,
        disease_case = disease,
        disease_control = control_diseases,
        disease_train = disease_train,
        disease_test = disease_test
      )
      
    }) |>
    set_names(c(classes, "Healthy"))
  
  
  roc_class_df <-
    map_df(c(classes, "Healthy"), function(class) {
      res <-  res_lr_classes[[class]]
      
      res$roc_curve |>
        mutate(AUC = res$AUC) |>
        mutate(
          Class = class,
          Disease = disease,
          Protein = protein,
          Combo = paste(disease, protein, sep = "_")
        )
    })
  
}

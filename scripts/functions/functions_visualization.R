#### Title: Functions for data visualization 
#### Author: María Bueno Álvez
#### Description: script collecting functions to plot data
#### Last edited : 12/08/2024

# Visualization packages
library(ggrepel)
library(ggbeeswarm)
library(patchwork)
library(ggplotify)
library(pheatmap)
library(tidygraph)
library(ggraph)
library(embed)
library(ggridges)
library(eulerr)
library(UpSetR)
library(tidytext)

# Function to generate volcano plot from differential expression results                                                                                                                                  
plot_volcano <- function(de_results, cutoff = 0) {
  
  labels <- 
    de_results |> 
    top_n(n = 10, wt = -log10(adj.P.Val)) 
  
  volcano_plot <- 
    de_results |> 
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = sig, label = Assay)) +
    geom_point(size = 1, alpha = 0.4, show.legend = F) + 
    geom_text_repel(data = labels, size = 2, show.legend = F) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = -cutoff, linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = cutoff, linetype = 'dashed', color = "darkgrey") +
    scale_color_manual(values = pal_de) +
    theme_hpa() +   
    theme(axis.text = element_text(size = 8),
          legend.position = "top") 
  
  return(volcano_plot)
}

# Function to generate confusion matrix plot                                                                                                                                
plot_cm <- function(confusion_matrix,
                    percentage = F) {
  
  ann_row <- 
    resource_meta |> 
    distinct(Class, Disease) |> 
    filter(Disease %in% c(include_diseases, "Healthy")) |> 
    rename(`True class` = Class) |> 
    column_to_rownames("Disease") 
  
  ann_col <- 
    resource_meta |> 
    distinct(Class, Disease) |> 
    filter(Disease %in% c(include_diseases, "Healthy")) |> 
    rename(`Predicted class` = Class) |> 
    column_to_rownames("Disease") 
  
  
  cm_dat <- 
    confusion_matrix$table |> 
    as.data.frame() |>
    group_by(Truth) |> 
    mutate(Truth = factor(Truth, levels = c(disease_class_order_ml, "Healthy")),
           Prediction = factor(Prediction, levels = c(disease_class_order_ml, "Healthy"))) 
  
  if(percentage == F) {
    dat <- 
      cm_dat |> 
      arrange(Truth, Prediction) |> 
      pivot_wider(names_from = Truth, 
                  values_from = Freq) |>
      column_to_rownames("Prediction") |> 
      t() |> 
      as.data.frame() 
    
    labels <- 
      cm_dat |> 
      mutate(Freq = ifelse(Freq > 0, Freq, ""),
             Freq = as.character(Freq)) |> 
      arrange(Truth, Prediction) |> 
      pivot_wider(names_from = Truth, 
                  values_from = Freq) |>
      column_to_rownames("Prediction") |> 
      t() |> 
      as.data.frame() 
    
    title <-  "Confusion matrix (n)"
    
    
  } else {
    
    dat <- 
      cm_dat |> 
      mutate(Freq = round((Freq/sum(Freq)) * 100, 0)) |> 
      arrange(Truth, Prediction) |> 
      pivot_wider(names_from = Truth, 
                  values_from = Freq) |>
      column_to_rownames("Prediction") |> 
      t() |> 
      as.data.frame() 
    
    labels <- 
      cm_dat |>
      group_by(Truth) |> 
      mutate(Freq = round((Freq/sum(Freq)) * 100, 0),
             Freq = ifelse(Freq > 10, Freq, ""),
             Freq = as.character(Freq)) |> 
      arrange(Truth, Prediction) |> 
      pivot_wider(names_from = Truth, 
                  values_from = Freq) |>
      column_to_rownames("Prediction") |> 
      t() |> 
      as.data.frame() 
    
    title <- "Confusion matrix (%)"
  }

  dat |> 
    pheatmap(cluster_rows = F, 
             annotation_row = ann_row,
             annotation_col = ann_col,
             cellwidth = 9,
             cellheight = 9, 
             annotation_colors = list("True class" = pal_class, "Predicted class" = pal_class),
             color = c("white", pal_heat),
             display_numbers = labels,
             cluster_cols = F) |> 
    as.ggplot() +
    coord_fixed() +
    theme(plot.title = element_text(face = "bold", size = rel(1))) +
    ggtitle(title) 
}

plot_top_proteins_seed <- function(protein_importance, 
                                   plot_color = "grey20",
                                   n = 25) {
  
  top_proteins <- 
    protein_importance|> 
    group_by(term) |> 
    summarise(abs_estimate = mean(abs(estimate)),
              n = n_distinct(seed)) |> 
    arrange(-n) |> 
    head(n)
  
  protein_importance |> 
    filter(term %in% top_proteins$term) |> 
    mutate(term = factor(term, levels = rev(top_proteins$term))) |> 
    ggplot(aes(fct_reorder(term, abs(estimate)), abs(estimate))) +
    geom_quasirandom(size = 0.5, color = plot_color) +
    geom_boxplot(alpha = 0.5, outlier.color = NA, fill = plot_color) +
    coord_flip() +
    theme_hpa() +
    xlab("") 
  
}

plot_boxplot <- function(proteins, 
                         type) {
  
  if(type == "sex") {
    
    resource_data |> 
      filter(Assay %in% proteins) |> 
      select(DAid, Assay, NPX) |> 
      left_join(resource_meta |> 
                  select(DAid, Sex), by = "DAid") |> 
      filter(!is.na(Sex)) |> 
      mutate(Assay = factor(Assay, levels = proteins)) |>
      ggplot(aes(Sex, NPX, fill = Sex, color = Sex)) +
      geom_quasirandom() +
      geom_boxplot(alpha = 0.5, outlier.color = NA, color = "black") +
      facet_wrap(~Assay, nrow = 1, scales = "free_y") +
      scale_color_manual(values = pal_sex) +
      scale_fill_manual(values = pal_sex) +
      theme_hpa() 
    
  } else if (type == "age") {
    
    resource_data |> 
      filter(Assay %in% proteins) |> 
      select(DAid, Assay, NPX) |> 
      left_join(resource_meta |> 
                  select(DAid, Age), by = "DAid") |> 
      filter(!is.na(Age)) |> 
      mutate(Assay = factor(Assay, levels = proteins)) |>
      ggplot(aes(Age, NPX)) +
      geom_point(size = 0.5) +
      facet_wrap(~Assay, nrow = 1, scales = "free_y") +
      theme_hpa() 
    
  } else if (type == "bmi") {
    
    resource_data |> 
      filter(Assay %in% proteins) |> 
      select(DAid, Assay, NPX) |> 
      left_join(resource_meta |> 
                  select(DAid, BMI), by = "DAid") |> 
      filter(!is.na(BMI),
             BMI != 0) |> 
      mutate(Assay = factor(Assay, levels = proteins)) |>
      ggplot(aes(BMI, NPX)) +
      geom_point(size = 0.5) +
      facet_wrap(~Assay, nrow = 1, scales = "free_y") +
      theme_hpa() 
    
  } else {
    cat("Invalid type")
  }
}


plot_heatmap <- function(protein_importance, 
                         type) {
  
  top_proteins <- 
    protein_importance|> 
    count(term) |> 
    filter(n == 100) |> 
    pull(term)
  
  if(type == "sex") {
    
    annotation <- 
      resource_meta |> 
      filter(!is.na(Sex)) |> 
      select(DAid, Sex) |> 
      column_to_rownames("DAid")
    
    resource_data |> 
      filter(DAid %in% rownames(annotation),
             Assay %in% top_proteins) |> 
      mutate(NPX = scales::rescale(NPX, to = c(0,1))) |> 
      select(DAid, NPX, Assay) |> 
      pivot_wider(names_from = Assay, 
                  values_from = NPX) |> 
      column_to_rownames("DAid") |> 
      pheatmap(annotation_row = annotation,
               show_rownames = F,
               color = pal_heat,
               clustering_method = "ward.D2") |> 
      as.ggplot()
    
  } else if (type == "age") {
    
    annotation <- 
      resource_meta |> 
      filter(!is.na(Age)) |> 
      select(DAid, Age) |> 
      column_to_rownames("DAid")
    
    resource_data |> 
      filter(DAid %in% rownames(annotation),
             Assay %in% top_proteins) |> 
      mutate(NPX = scales::rescale(NPX, to = c(0,1))) |> 
      select(DAid, NPX, Assay) |> 
      pivot_wider(names_from = Assay, 
                  values_from = NPX) |> 
      column_to_rownames("DAid") |> 
      pheatmap(annotation_row = annotation,
               show_rownames = F,
               color = pal_heat,
               clustering_method = "ward.D2") |> 
      as.ggplot() 
    
  } else if (type == "bmi") {
    
    annotation <- 
      resource_meta |> 
      filter(!is.na(BMI),
             BMI != 0) |> 
      select(DAid, BMI) |> 
      column_to_rownames("DAid")
    
    resource_data |> 
      filter(DAid %in% rownames(annotation),
             Assay %in% top_proteins) |> 
      mutate(NPX = scales::rescale(NPX, to = c(0,1))) |> 
      select(DAid, NPX, Assay) |> 
      pivot_wider(names_from = Assay, 
                  values_from = NPX) |> 
      column_to_rownames("DAid") |> 
      as.data.frame() |> 
      pheatmap(annotation_row = annotation,
               show_rownames = F,
               color = pal_heat,
               clustering_method = "ward.D2") |> 
      as.ggplot() 
    
  } else {
    cat("Invalid type")
  }
}

plot_n_importance <- function(importance_data) {
  importance_data |>
    group_by(class, seed) |> 
    summarize(n = n_distinct(term)) |> 
    left_join(resource_meta |> distinct(Disease, Class), by = c("class" = "Disease")) |>
    mutate(class = factor(class, levels = disease_class_order)) |> 
    ggplot(aes(class, n, color = Class, fill = Class)) +
    geom_quasirandom(size = 0.5) +
    geom_boxplot(width = 0.5, outlier.color = NA, color = "grey20", alpha = 0.5) +
    scale_color_manual(values = pal_class) + 
    scale_fill_manual(values = pal_class) + 
    theme_hpa(angled = T)
}

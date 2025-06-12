# other comparisons


## Compared to healthy turkey

```{r}
# Prepare data & metadata
patients_de <- 
  da_meta |> 
  mutate(Sex = case_when(Gender == 0 ~ "Female",
                         Gender == 1 ~ "Male")) |> 
  select(DAid, Disease, Age, Sex, BMI) |> 
  filter(Disease %in% c("Pancreatic cancer", "Healthy Turkey")) #|> # 941 -> 931 (not many exclusions)
#filter(!is.na(Sex),
#      !Age < 10,
#     !is.na(Age))

wide_data <- 
  da_data |> 
  #  select(-Disease) |> 
  filter(DAid %in% patients_de$DAid) |> 
  select(DAid, Assay, NPX) |> 
  pivot_wider(values_from = "NPX", names_from = "Assay")

# Run differential expression analysis
de_healthy_turkey<- 
  de_limma_disease(data_wide = wide_data, 
                   metadata = patients_de,
                   disease = "Pancreatic cancer",
                   correct = F) |> 
  mutate(sig = case_when(adj.P.Val < 0.05 & logFC < -1 ~ "significant down",
                         adj.P.Val < 0.05 & logFC > 1 ~ "significant up", 
                         T ~ "not significant"))

# Generate volcano plot
volcano_healthy <- plot_volcano(de_healthy_turkey)

de_healthy |> 
  select(Assay, logFC) |> 
  left_join(de_healthy_turkey |> 
              select(Assay, logFC_turkey = logFC)) |> 
  ggplot(aes(logFC, logFC_turkey)) +
  geom_point() +
  geom_text_repel(aes(label = Assay)) +
  theme_hpa()

```


## Compared to Wellness

```{r}
# Prepare data & metadata
patients_de <- 
  da_meta |> 
  mutate(Sex = case_when(Gender == 0 ~ "Female",
                         Gender == 1 ~ "Male")) |> 
  select(DAid, Disease, Age, Sex, BMI) |> 
  filter(Disease %in% c("Pancreatic cancer", "Wellness")) #|> # 941 -> 931 (not many exclusions)
#filter(!is.na(Sex),
#      !Age < 10,
#     !is.na(Age))

wide_data <- 
  da_data |> 
  #  select(-Disease) |> 
  filter(DAid %in% patients_de$DAid) |> 
  select(DAid, Assay, NPX) |> 
  pivot_wider(values_from = "NPX", names_from = "Assay")

# Run differential expression analysis
de_healthy_wellness <- 
  de_limma_disease(data_wide = wide_data, 
                   metadata = patients_de,
                   disease = "Pancreatic cancer",
                   correct = F) |> 
  mutate(sig = case_when(adj.P.Val < 0.05 & logFC < -1 ~ "significant down",
                         adj.P.Val < 0.05 & logFC > 1 ~ "significant up", 
                         T ~ "not significant"))

# Generate volcano plot
plot_volcano(de_healthy_wellness)

de_healthy |> 
  select(Assay, logFC) |> 
  left_join(de_healthy_wellness |> 
              select(Assay, logFC_wellness = logFC)) |> 
  ggplot(aes(logFC, logFC_wellness)) +
  geom_point() +
  geom_text_repel(aes(label = Assay)) +
  theme_hpa()

```


## Compared to EpiHealth

```{r}
# Prepare data & metadata
patients_de <- 
  da_meta |> 
  mutate(Sex = case_when(Gender == 0 ~ "Female",
                         Gender == 1 ~ "Male")) |> 
  select(DAid, Disease, Age, Sex, BMI) |> 
  filter(Disease %in% c("Pancreatic cancer", "EpiHealth")) #|> # 941 -> 931 (not many exclusions)
#filter(!is.na(Sex),
#      !Age < 10,
#     !is.na(Age))

wide_data <- 
  da_data |> 
  #  select(-Disease) |> 
  filter(DAid %in% patients_de$DAid) |> 
  select(DAid, Assay, NPX) |> 
  pivot_wider(values_from = "NPX", names_from = "Assay")

# Run differential expression analysis
de_healthy_epi <- 
  de_limma_disease(data_wide = wide_data, 
                   metadata = patients_de,
                   disease = "Pancreatic cancer",
                   correct = F) |> 
  mutate(sig = case_when(adj.P.Val < 0.05 & logFC < -1 ~ "significant down",
                         adj.P.Val < 0.05 & logFC > 1 ~ "significant up", 
                         T ~ "not significant"))

# Generate volcano plot
plot_volcano(de_healthy_epi)

de_healthy |> 
  select(Assay, logFC) |> 
  left_join(de_healthy_epi |> 
              select(Assay, logFC_epi = logFC)) |> 
  ggplot(aes(logFC, logFC_epi)) +
  geom_point() +
  geom_text_repel(aes(label = Assay)) +
  theme_hpa()

plot_volcano(de_ukb)


de_healthy |> 
  select(Assay, logFC) |> 
  left_join(de_ukb |> 
              select(Assay, logFC_ukb = logFC)) |> 
  ggplot(aes(logFC, logFC_ukb)) +
  geom_point() +
  geom_text_repel(aes(label = Assay)) +
  theme_hpa()


```

## UKB data

```{r ukb data}
ukb_healthy_data <- read_table("data/ukb/UKBB_healthy_Olink.txt")
ukb_healthy_metadata <-  read_table("data/UKB/UKBB_healthy_metatable.txt")
ukb_pancreatic_data <- read_table("data/ukb/pancreas_cancer_olink.txt")
ukb_pancreatic_metadata <- read_table("data/ukb/pancreas_cancer_meta.txt")

# Formatting
ukb_meta <- 
  ukb_healthy_metadata |> 
  select(id = eid, Sex, Age = Age_sample_collect) |> 
  mutate(Disease = "Healthy") |> 
  bind_rows(ukb_pancreatic_metadata |> 
              select(id = eid, Sex, Age = Age_sample_collect, cancer_diagnosis = Cancer_age_diagnosis) |> 
              mutate(Disease = "Pancreatic cancer") )

ukb_data <- 
  ukb_healthy_data |> 
  separate(meaning, c("Assay", "Description"), sep = ";") |> 
  select(id = eid, Assay, NPX = result) |> 
  mutate(Disease = "Healthy") |> 
  bind_rows(ukb_pancreatic_data |> 
              separate(meaning, c("Assay", "Description"), sep = ";") |> 
              select(id = eid, Assay, NPX = result) |> 
              mutate(Disease = "Pancreatic cancer"))

# Check duplicated healthy samples
dublicated_healthy_ids <-
  ukb_data |> 
  group_by(id, Assay) |> 
  summarise(n = n()) |> 
  filter(n>1) |> 
  distinct(id)

ukb_meta <- 
  ukb_meta |> 
  filter(!id %in% dublicated_healthy_ids$id) 

ukb_data <- 
  ukb_data |> 
  filter(!id %in% dublicated_healthy_ids$id) 

ukb_meta_ext <- 
  ukb_meta |> 
  distinct(id, .keep_all = T) |> 
  mutate(Difference = cancer_diagnosis - Age,
         Group = case_when(Difference > 7 ~ "> 7 years before",
                           Difference < 7 & Difference >= 5 ~ "5-7 years before",
                           Difference < 5 & Difference >= 3 ~ "3-5 years before",
                           Difference < 3 & Difference >= 1 ~ "1-3 years before",
                           Difference < 1 & Difference >= -1 ~ "1 year before/after",
                           Difference < -1 ~ "> 1 year after",
                           Disease == "Healthy" ~ "Healthy"),
         Group = factor(Group, levels = names(pal_ukb))) 
```

# General




# Signature against UKB

## glmnet

```{r}
ukb_data_model <-
  ukb_data |>
  select(DAid = id, Disease, Assay, NPX) |> 
  pivot_wider(names_from = "Assay", values_from = "NPX") |> 
  mutate(DAid = as.character(DAid))

final_wf <- 
  final_rf_fit %>%
  extract_workflow()

final_wf |> tidy() |> arrange(-abs(estimate)) |> filter(abs(estimate) > 0)

ukb_glmnet_fit <- 
  predict(final_wf, ukb_data_model, type = "prob") |> 
  mutate(True_class = ukb_data_model$Disease,
         DAid = ukb_data_model$DAid)

ukb_glmnet_fit |> 
  left_join(ukb_meta |> 
              select(-Disease) |> 
              mutate(Difference = round(cancer_diagnosis - Age),
                     DAid = as.character(id)), by = "DAid") |> 
  arrange(-`.pred_Pancreatic cancer`)

ukb_glmnet_fit |> 
  left_join(ukb_meta |> 
              select(-Disease) |> 
              mutate(Difference = round(cancer_diagnosis - Age),
                     DAid = as.character(id)), by = "DAid") |> 
  ggplot(aes(True_class, `.pred_Pancreatic cancer`, color = Difference)) +
  geom_quasirandom() +
  geom_violin(fill = NA) +
  scale_color_viridis_c(begin = 1, end = 0) +
  theme_hpa()

ggsave(savepath("prob_da_model_ukb.png"), h = 5, w = 5) 


ukb_glmnet_fit |> 
  left_join(ukb_meta |> 
              select(-Disease) |> 
              mutate(Difference = round(cancer_diagnosis - Age),
                     DAid = as.character(id)), by = "DAid") |> 
  mutate(Group = case_when(Difference > 7 ~ "> 7 years before",
                           Difference < 7 & Difference >= 5 ~ "5-7 years before",
                           Difference < 5 & Difference >= 3 ~ "3-5 years before",
                           Difference < 3 & Difference >= 1 ~ "1-3 years before",
                           Difference < 1 & Difference >= -1 ~ "1 year before/after",
                           Difference < -1 ~ "> 1 year after",
                           True_class == "Healthy" ~ "Healthy"),
         Group = factor(Group, levels = names(pal_ukb))) |> 
  ggplot(aes(Group, `.pred_Pancreatic cancer`, color = Difference)) +
  geom_quasirandom() +
  geom_violin(fill = NA, ) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black") +
  scale_color_viridis_c(begin = 1, end = 0) +
  theme_hpa() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(savepath("prob_da_model_ukb_ext.png"), h = 6, w = 8) 


# rf


```

# Multiclassification

```{r}

```


# Other cancers as controls

```{r}

```




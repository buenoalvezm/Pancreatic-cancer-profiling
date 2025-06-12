#### Title: Themes & palettes
#### Author: María Bueno Álvez
#### Description: script collecting themes & palettes
#### Last edited : 12/08/2024

# Packages for themes & palettes
library(RColorBrewer)

# Palettes
pal_de <-
  c("not significant" = "#D3D3D3",
    "significant up" = "#FF7176",
    "significant down" = "#92C9DA")

pal_age_sex_bmi <- 
  c("Age" = "#75C8AE", #459395
    "Sex" = "#EB7C6A", 
    "BMI" = "#F7B84D",
    "Disease" = "#635C94")


pal_anova <- c("Age" = "#75C8AE",    
               "Sex" = "#EB7C6A", 
               "BMI" = "#F7B84D",    
               "Disease" = "#EEE2D1")    


pal_ukb <- 
  c("Healthy" = "grey",
    "> 7 years before" = "#4D7654", 
    "5-7 years before" = "#748B5F", 
    "3-5 years before" = "#E3D6A0", 
    "1-3 years before" = "#C78240", 
    "1 year before/after" = "#A42F2E", 
    "> 1 year after" = "#510402")

pal_ukb <- 
  c("Healthy" = "grey",
    "> 7 years before" = "#4D7654", 
    "5-7 years before" = "#748B5F", 
    "3-5 years before" = "#E3D6A0", 
    "1-3 years before" = "#C78240", 
    "1 year before/after" = "#A42F2E", 
    "> 1 year after" = "#510402")

pal_ukb_2 <- 
  c("Healthy" = "grey",
    "> 7 years before" = "#4D7654", 
    "5-7 years before" = "#748B5F", 
    "3-5 years before" = "#E3D6A0", 
    "1-3 years before" = "#C78240", 
    "1 year before/after" = "#A42F2E", 
    "1-3 years after" = "#802020", 
    "> 3 years after" = "#510402")


pal_heat <- colorRampPalette(c("#FFF3E0",  "#FFCC80", "#B7410E"))(100)

pal_corr <- colorRampPalette(rev(brewer.pal(n=10, name="RdBu")))(200)

rd_bu_continuous <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
pal_cor_2 <- rd_bu_continuous[15:85]

pal_sex <-   c("F"= "#8a72be", "M" = "#A9D0EF")

getPalette3 = colorRampPalette(brewer.pal(8, "Set2"))
pal_class<-getPalette3(8)
names(pal_class)<-c("Psychiatric","Cardiovascular","Cancer","Autoimmune","Pediatric","Infection","Metabolic","Healthy") 
class_order <- c("Healthy", "Cardiovascular","Metabolic","Cancer","Psychiatric","Autoimmune","Infection","Pediatric")


pal_secreted <- c("Secreted to blood" = "#B30000", 
                  "Secreted in brain" = "#FFDD00", 
                  "Secreted to digestive system" = "#1280C4", 
                  "Secreted in male reproductive system" = "#95D4F5", 
                  "Secreted in female reproductive system" = "#F8BDD7", 
                  "Secreted to extracellular matrix"  = "#7F6A9C", 
                  "Secreted in other tissues" = "#FFD480", 
                  "Secreted - unknown location" = "#A1A8AA", 
                  "Intracellular and membrane" = "#F9A266", 
                  "Unknown" = "grey80")

pal_specificity <- c( "Tissue enriched" = "#e41a1c",
                      "Group enriched" = "#FF9D00",
                      "Tissue enhanced" = "#984ea3",
                      "Low tissue specificity" = "grey40",
                      "not detected " = "grey")
  

consensus_colors <- 
  read_tsv("../DA/data/hpa/colors_consensus.tsv")

consensus_colors <-
  consensus_colors %>% 
  mutate(color = ifelse(color == "#DE6C7D", "#B38C6D", color),
         color = ifelse(color == "#A1A8AA", "#DE6C7D", color)) 

pal_tissue <- 
  consensus_colors |> 
  filter(dataset_id == "tissue") |> 
  select(-1) |> 
  bind_rows(tibble(sample = c("amygdala", "appendix", "basal ganglia",
                              "cerebellum","cerebral cortex","colon",
                              "duodenum","hippocampal formation", "hypothalamus",
                              "lymph node","midbrain","rectum"   ,
                              "small intestine","spinal cord","spleen",
                              "thymus","tonsil"), 
                   color = c("#FFDD00","#DE6C7D","#FFDD00",
                             "#FFDD00","#FFDD00", "#1280C4",
                             "#1280C4", "#FFDD00","#FFDD00",
                             "#DE6C7D","#FFDD00", "#1280C4",
                             "#1280C4","#FFDD00", "#DE6C7D",
                             "#DE6C7D", "#DE6C7D"))) |> 
  deframe()

pal_tissue <- c(pal_tissue, "Not specific" = "grey80")

# Levels
female_diseases <- c("Breast cancer", "Breast ductal carcinoma in situ", "Cervical cancer", "Endometrial cancer", "Ovarian cancer")
male_diseases <-  c("Prostate cancer", "Abdominal aortic aneurysm")
pediatric_diseases <- c("Pediatric CNS tumor",
                        "Pediatric bone tumor",
                        "Pediatric retinoblastoma",
                        "Pediatric lymphoma", 
                        "Pediatric neuroblastoma",
                        "Pediatric sarcoma",
                        "Pediatric kidney tumor" , 
                        "Pediatric diffuse astrocytic and oligodendro. tumor",
                        "Pediatric long COVID",
                        "Pediatric systemic inflammatory disease")
# HPA theme
theme_hpa <- 
  function(angled = F, axis_x = T, axis_y = T, facet_title = T) {
    t <- 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = rel(1), hjust = 0.5),
        plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.8)),
        legend.key.size= unit(0.7, "cm"),
        legend.title = element_text(size=rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="grey90",fill="grey90"),
        strip.text = element_text(face="bold")
      )
    
    if(angled) {
      t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    }
    
    if(axis_x == F) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())
    } 
    
    if(axis_y == F) {
      t <- t +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank())
    }
    if(facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }



# Other data

#pancreatic_meta <- readxl::read_excel("data/clinical_metadata/Pancreatic/UPApatients.xlsx")
# Read in data (MEDECA)
# data_medeca <- import_df("../MEDECA/data/processed/final_data/final_olink_data.csv")  
# metadata_medeca <- import_df("data/processed/final_data/combined_metadata.csv")  
# meta_medeca <- import_df("data/processed/final_data/medeca_metadata.csv")  
# meta_allvos <- import_df("data/processed/final_data/allvos_metadata.csv")  
# Read in data (Phase 2)
#data_phase_2 <- read_tsv(file = "../Human-disease-blood-atlas/data/final_data/data_phase2_batch4_curated_20250425.tsv") 
#manifest <- import_df("../Human-disease-blood-atlas/data/samples_2025-05-12.xlsx")

# Phase 2
meta_phase2 <- 
  manifest |> 
  filter(DAid %in% data_phase_2$DAid)

dat_p2 <- hd_initialize(dat = data_phase_2,
                        metadata = meta_phase2)

qc_p2 <- hd_qc_summary(dat_p2,
                       variable = "Disease")
# MEDECA

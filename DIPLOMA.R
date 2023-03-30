install.packages('dplyr')
library(dplyr)
FMI = read.csv("FMI variants.csv")
frame = read.csv('FMI variants.csv')
MAVE_VAMP = read.csv('MAVE_VAMP_PTEN_scores.csv')
OncoKB_CKB = read.csv('OncoKB-CKB PTEN annotation.csv')
FMI_samples = read.csv("FMI samples ID's.csv")

MAVE_VAMP = MAVE_VAMP %>%    
  mutate(MAVE_status = ifelse(MAVE_score < -1.1, "LoF", "WT"), VAMP_status = ifelse(VAMP_score < 0.7, "LoF", "WT"))  


frame = frame %>% 
  mutate(MAVE_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "LoF")$amino_acid_change, "LoF", NA), MAVE_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "WT")$amino_acid_change, "WT", MAVE_status))

frame = frame %>%    
  mutate(VAMP_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "LoF")$amino_acid_change, "LoF", NA), VAMP_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "WT")$amino_acid_change, "WT", VAMP_status))


frame = frame %>% 
  mutate(oncoKB = ifelse(gene == 'PTEN' & amino_acid_change %in% filter(OncoKB_CKB, oncoKB == 'lof')$amino_acid_change, 'LoF', NA), oncoKB = ifelse(gene == 'PTEN' & amino_acid_change %in% filter(OncoKB_CKB, oncoKB == 'wt')$amino_acid_change, 'WT', oncoKB))


frame = frame %>% 
  mutate(CKB = ifelse(gene == 'PTEN' & amino_acid_change %in% filter(OncoKB_CKB, CKB == 'lof')$amino_acid_change, 'LoF', NA), CKB)


FMI_samples = frame %>%    group_by(sample_id) %>%   summarise(PTEN_presence_VAMP_LoF = sum(gene == "PTEN" & VAMP_status == "LoF"),             PTEN_presence_VAMP_WT = sum(gene == "PTEN" & VAMP_status == "WT"),             PTEN_presence_MAVE_LoF = sum(gene == "PTEN" & VAMP_status == "LoF"),             PTEN_presence_MAVE_WT = sum(gene == "PTEN" & VAMP_status == "WT"), PTEN_presence_OncoKB_Lof = sum(gene == "PTEN" & oncoKB == "LoF"), PTEN_presence_CKB_LoF = sum(gene == "PTEN" & CKB == "WT"), 
                                                               PTEN_presence_ANY_LoF = sum(VAMP_status == "LoF", MAVE_status == "LoF", oncoKB == "LoF", CKB == "LoF" & gene == "PTEN"),
                                                               PTEN_presence_ANY_WT = sum(VAMP_status == "WT", MAVE_status == "WT" & gene == "PTEN"),
                                                               APC_presence = sum(gene == "APC"),         
                                                               SMAD4_presence = sum(gene == "SMAD4"),             
                                                               TP53_presence = sum(gene == "TP53"),
                                                               KRAS_presence = sum(gene == "KRAS"),
                                                               NRAS_presence = sum(gene == "NRAS"))
                                                               

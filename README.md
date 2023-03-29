# KFU_Anya_diploma
install.packages('dplyr')
library(dplyr)
FMI = read.csv("FMI variants.csv")
frame = read.csv('FMI variants.csv')
MAVE_VAMP = read.csv('MAVE_VAMP_PTEN_scores.csv')
OncoKB_CKB = read.csv('OncoKB-CKB PTEN annotation.csv')
FMI.samples = read.csv("FMI samples ID's.csv")

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


PTEN_presence <- frame %>% 
  group_by(MAVE_status, VAMP_status, oncoKB, CKB, sample_id, gene) %>% 
  summarise(MAVE_status == 'LoF' |
              VAMP_status == 'LoF' |
              oncoKB == 'LoF' |
              CKB == 'Lof' &
              gene == 'PTEN')

# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c("tidyverse", "data.table", "readxl", "reshape2", "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename = dplyr::rename
select = dplyr::select
filter = dplyr::filter
group_by = dplyr::group_by
mutate = dplyr::mutate

FMI <- fread("FMI variants.csv") %>% as_tibble()
frame <- fread("FMI variants.csv") %>% as_tibble()
MAVE_VAMP <- fread("MAVE_VAMP_PTEN_scores.csv") %>% as_tibble()
OncoKB_CKB <-
  fread("OncoKB-CKB PTEN annotation.csv") %>% as_tibble()
FMI_samples <- fread("FMI samples ID's.csv") %>% as_tibble()
FMI_samples = filter(FMI_samples, msi == 'MT-L')

MAVE_VAMP <- MAVE_VAMP %>%
  mutate(
    MAVE_status = ifelse(MAVE_score < -1.1, "LoF", "WT"),
    VAMP_status = ifelse(VAMP_score < 0.7, "LoF", "WT"),
    oncoKB = ifelse(
      amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    oncoKB = ifelse(
      amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "WT")$amino_acid_change,
      "WT",
      oncoKB
    ),
    CKB = ifelse(
      amino_acid_change %in% filter(OncoKB_CKB, CKB == "LoF")$amino_acid_change,
      "LoF",
      NA
    )
  )


frame <- frame %>%
  mutate(
    MAVE_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    MAVE_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "WT")$amino_acid_change,
      "WT",
      MAVE_status
    )
  )

frame <- frame %>%
  mutate(
    VAMP_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    VAMP_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "WT")$amino_acid_change,
      "WT",
      VAMP_status
    )
  )


frame <- frame %>%
  mutate(
    oncoKB = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    oncoKB = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "WT")$amino_acid_change,
      "WT",
      oncoKB
    )
  )


frame <- frame %>%
  mutate(CKB = ifelse(
    gene == "PTEN" &
      amino_acid_change %in% filter(OncoKB_CKB, CKB == "LoF")$amino_acid_change,
    "LoF",
    NA
  ))


# Add samples with 0 mutations from FMI samples ID's.csv using left_join/right_join
# NA -> 0 using replace_na 

FMI_samples <- frame %>%
  group_by(sample_id) %>%
  summarise(
    PTEN_presence_VAMP_LoF = sign(sum(gene == "PTEN" &
                                   VAMP_status == "LoF")),
    PTEN_presence_VAMP_WT = sign(sum(gene == "PTEN" &
                                  VAMP_status == "WT")),
    PTEN_presence_MAVE_LoF = sign(sum(gene == "PTEN" &
                                   VAMP_status == "LoF")),
    PTEN_presence_MAVE_WT = sign(sum(gene == "PTEN" &
                                  VAMP_status == "WT")),
    PTEN_presence_OncoKB_LoF = sign(sum(gene == "PTEN" &
                                     oncoKB == "LoF")),
    PTEN_presence_CKB_LoF = sign(sum(gene == "PTEN" &
                                  CKB == "LoF")),
    PTEN_presence_ANY_LoF = sign(sum(gene == "PTEN" &
                                  VAMP_status == "LoF", MAVE_status == "LoF")),
    PTEN_presence_ANY_WT = sign(sum(gene == "PTEN" &
                                 VAMP_status == "WT", MAVE_status == "WT")),
    MAVE_VAMP_LoF = sign(sum(gene == "PTEN" &
                          MAVE_status == "LoF" &
                          VAMP_status == "LoF")),
    MAVE_VAMP_WT = sign(sum(gene == "PTEN" &
                         MAVE_status == "WT" &
                         VAMP_status == "WT")),
    ALL_LoF = sign(sum(
      MAVE_status == "LoF" &
        VAMP_status == "LoF" &
        oncoKB == "LoF" &
        CKB == "LoF"
    )),
    APC_presence = sign(sum(gene == "APC")),
    SMAD4_presence = sign(sum(gene == "SMAD4")),
    TP53_presence = sign(sum(gene == "TP53")),
    KRAS_presence = sign(sum(gene == "KRAS")),
    NRAS_presence = sign(sum(gene == "NRAS"))) %>% 
  right_join(y = FMI_samples, by = c("sample_id")) %>% 
  mutate(across(.cols = everything(), .fns = function(x) replace_na(x, 0)))

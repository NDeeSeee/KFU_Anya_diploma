# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c("tidyverse", "data.table", "readxl", "reshape2", "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate

FMI <- fread("FMI variants.csv") %>% as_tibble()

DB_consensus <-
  fread("Database_consensus_LoF_PTEN.csv") %>% as_tibble()

FMI_samples <- fread("FMI samples ID's.csv") %>% as_tibble()

FMI_samples <- filter(FMI_samples, msi == "MT-L")

MAVE_VAMP <- fread("MAVE_VAMP_PTEN_scores.csv") %>%
  as_tibble() %>%
  mutate(
    MAVE_status = ifelse(MAVE_score < -1.1, "LoF", "WT"),
    VAMP_status = ifelse(VAMP_score < 0.7, "LoF", "WT"),
    DB_status = ifelse(
      amino_acid_change %in% filter(DB_consensus, DB_consensus_assigned_effect == "LoF")$aa_change,
      "LoF",
      NA
    ),
    DB_status = ifelse(
      amino_acid_change %in% filter(DB_consensus, DB_consensus_assigned_effect == "WT")$aa_change,
      "WT",
      DB_status
    )
  )


frame <- fread("FMI variants.csv") %>% 
  as_tibble() %>% 
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
    ),
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
    ),
    DB_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, DB_status == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    DB_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, DB_status == "WT")$amino_acid_change,
      "WT",
      DB_status
    ),
    DB_MV_consensus = case_when(DB_status == "LoF" ~ "LoF",
                                DB_status == "WT" ~ "WT",
                                is.na(DB_status) & (VAMP_status == "LoF" | MAVE_status == "LoF") ~ "LoF",
                                is.na(DB_status) & (VAMP_status == "WT" & MAVE_status == "WT") ~ "WT"),
    ultimate_consensus = case_when(DB_MV_consensus == "LoF" ~ "LoF",
                                   DB_MV_consensus == "WT" ~ "WT",
                                   gene == "PTEN" & (variantClass_v2 == "truncation" | codingType_SV == "splice") ~ "LoF")
  )


# Add samples with 0 mutations from FMI samples ID's.csv using left_join/right_join
# NA -> 0 using replace_na

frame_grouped <- frame %>%
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
    PTEN_presence_DB_LoF = sign(sum(gene == "PTEN" &
                                      DB_status == "LoF")),
    PTEN_presence_DB_WT = sign(sum(gene == "PTEN" &
                                      DB_status == "WT")),
    PTEN_presence_DB_MV_LoF = sign(sum(gene == "PTEN" &
                                      DB_MV_consensus == "LoF")),
    PTEN_presence_DB_MV_WT = sign(sum(gene == "PTEN" &
                                      DB_MV_consensus == "WT")),
    PTEN_presence_ultimate_LoF = sign(sum(gene == "PTEN" &
                                      ultimate_consensus == "LoF")),
    PTEN_presence_ultimate_WT = sign(sum(gene == "PTEN" &
                                      ultimate_consensus == "WT")),
  
    APC_presence = sign(sum(gene == "APC")),
    SMAD4_presence = sign(sum(gene == "SMAD4")),
    TP53_presence = sign(sum(gene == "TP53")),
    KRAS_presence = sign(sum(gene == "KRAS")),
    NRAS_presence = sign(sum(gene == "NRAS")),
    PIK3CA_presence = sign(sum(gene == "PIK3CA"))
  ) %>%
  right_join(y = FMI_samples, by = c("sample_id")) %>%
  mutate(across(
    .cols = everything(),
    .fns = function(x)
      replace_na(x, 0)
  ))

# PTEN column names
pten_types = names(frame_grouped)[which(str_detect(names(frame_grouped), "PTEN"))]

# Other genes
gene_names = names(frame_grouped)[which(!str_detect(names(frame_grouped), "PTEN") & str_detect(names(frame_grouped), "presence"))]

# Start
gene_name_vector = pten_type_vector = odds_ratio_vector = upper_conf_vector = lower_conf_vector = p_value_vector = c()

for (gene in gene_names) {
  for (pten_type in pten_types) {
    if (all(dim(table(select(frame_grouped, pten_type, gene))) == 2)) {
      f_test_results = fisher.test(table(select(frame_grouped, pten_type, gene)))
      gene_name_vector = c(gene_name_vector, gene)
      pten_type_vector = c(pten_type_vector, pten_type)
      odds_ratio_vector = c(odds_ratio_vector, f_test_results$estimate)
      upper_conf_vector = c(upper_conf_vector, f_test_results$conf.int[2])
      lower_conf_vector = c(lower_conf_vector, f_test_results$conf.int[1])
      p_value_vector = c(p_value_vector, f_test_results$p.value)
    }
  }
}

final_df = tibble(gene_name = gene_name_vector,
                  pten_type = pten_type_vector,
                  odds_ratio = odds_ratio_vector,
                  upper_conf = upper_conf_vector,
                  lower_conf = lower_conf_vector,
                  p_value = p_value_vector) %>% 
  mutate(odds_ratio = log2(odds_ratio),
         upper_conf = log2(upper_conf),
         lower_conf = log2(lower_conf))

for (gene in gene_names) {
  final_df_temp = filter(final_df, gene_name == gene) %>% 
    mutate(gene_name = str_remove_all(gene_name, "_presence"),
           pten_type = str_remove_all(pten_type, "PTEN_presence_"),
           strata = ifelse(str_detect(pten_type, "LoF"), "LoF", "WT"),
           strata = as.factor(strata))
  
  ggplot(final_df_temp, aes(x = pten_type, y = odds_ratio, col = strata)) +
    geom_hline(yintercept = 0, col = "gray35") +
    geom_point(size = 3, position = position_dodge(0.4)) +
    geom_errorbar(
      aes(ymin = lower_conf, ymax = upper_conf, col = strata),
      width = .2,
      position = position_dodge(0.4),
      linewidth = 1
    ) +
    theme_minimal() +
    theme(legend.title = element_blank(), text = element_text(angle = 25))
  
  ggsave(filename = paste0(gene, ".png"),
         dpi = 400,
         width = 5,
         height = 5,
         units = "in")
}

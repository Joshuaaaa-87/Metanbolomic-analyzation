

# ============================================================
# Filter KEGG p005 results: Best + Expanded
# Keep KEGG Metabolism category only
# ============================================================

library(tidyverse)
library(stringr)
library(readr)

#----1. Define input and output folders----

best_input_dir <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data/4_KEGG_Enrichment_BestKEGG_15ppm_top10"

expanded_input_dir <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data/4_KEGG_Enrichment_ExpandedKEGG_15ppm_top10"

output_dir <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/Johnathan/0_R/4.1_KEGG_final"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


#----2. Define files----

files <- c(
  Best = file.path(
    best_input_dir,
    "8_Concise_KEGG_pathway_table_log2FC1_p005_Best_KEGG.csv"
  ),
  Expanded = file.path(
    expanded_input_dir,
    "8_Concise_KEGG_pathway_table_log2FC1_p005_Expanded_KEGG.csv"
  )
)


#----3. Check whether files exist----

print(files)

if (!all(file.exists(files))) {
  stop(
    "Some input files do not exist. Please check the file names or folder path:\n",
    paste(files[!file.exists(files)], collapse = "\n")
  )
}


#----4. Read KEGG files----

kegg_raw <- purrr::imap_dfr(files, function(file, mapping_name) {
  readr::read_csv(file, show_col_types = FALSE) %>%
    dplyr::mutate(Mapping = mapping_name, .before = 1)
})


#----5. Classify KEGG pathway by KEGG map ID----

class_by_kegg_id <- function(map_id) {
  id_num <- suppressWarnings(as.integer(map_id))
  
  dplyr::case_when(
    !is.na(id_num) & id_num >= 0    & id_num < 1300 ~ "09100 Metabolism",
    !is.na(id_num) & id_num >= 2000 & id_num < 3000 ~ "09130 Environmental Information Processing",
    !is.na(id_num) & id_num >= 3000 & id_num < 4000 ~ "09120 Genetic Information Processing",
    !is.na(id_num) & id_num >= 4000 & id_num < 5000 ~ "09140/09150 Cellular Processes or Organismal Systems",
    !is.na(id_num) & id_num >= 5000 & id_num < 6000 ~ "09160 Human Diseases",
    !is.na(id_num) & id_num >= 7000 & id_num < 8000 ~ "09180 Drug Development",
    TRUE ~ "Other / Unknown"
  )
}


#----6. Annotate KEGG category----

kegg_annotated <- kegg_raw %>%
  dplyr::mutate(
    map_id = stringr::str_extract(Pathway_ID, "\\d{5}$"),
    KEGG_top_class = class_by_kegg_id(map_id),
    Is_Metabolism = KEGG_top_class == "09100 Metabolism"
  )


#----7. Keep metabolism-only pathways----

pathway_p_cutoff <- 0.05
pathway_fdr_cutoff <- 0.10

kegg_metabolism_only <- kegg_annotated %>%
  dplyr::filter(Is_Metabolism)

kegg_metabolism_only_pathwayP005 <- kegg_metabolism_only %>%
  dplyr::filter(P_value < pathway_p_cutoff) %>%
  dplyr::arrange(Mapping, Condition, P_value)

kegg_metabolism_only_pathwayFDR01 <- kegg_metabolism_only %>%
  dplyr::filter(FDR < pathway_fdr_cutoff) %>%
  dplyr::arrange(Mapping, Condition, FDR)


#----8. Export results----

readr::write_csv(
  kegg_annotated,
  file.path(output_dir, "p005_Best_and_Expanded_KEGG_annotated_all.csv")
)

readr::write_csv(
  kegg_metabolism_only,
  file.path(output_dir, "p005_Best_and_Expanded_KEGG_metabolism_only_all.csv")
)

readr::write_csv(
  kegg_metabolism_only_pathwayP005,
  file.path(output_dir, "p005_Best_and_Expanded_KEGG_metabolism_only_pathwayP005.csv")
)

readr::write_csv(
  kegg_metabolism_only_pathwayFDR01,
  file.path(output_dir, "p005_Best_and_Expanded_KEGG_metabolism_only_pathwayFDR01.csv")
)


#----9. Summary table----

summary_table <- kegg_annotated %>%
  dplyr::group_by(Mapping, Condition) %>%
  dplyr::summarise(
    Total_pathways = dplyr::n(),
    Metabolism_pathways = sum(Is_Metabolism, na.rm = TRUE),
    Metabolism_pathway_P005 = sum(Is_Metabolism & P_value < pathway_p_cutoff, na.rm = TRUE),
    Metabolism_pathway_FDR01 = sum(Is_Metabolism & FDR < pathway_fdr_cutoff, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  summary_table,
  file.path(output_dir, "p005_Best_and_Expanded_KEGG_metabolism_summary.csv")
)

print(summary_table)

message("Done. Output files saved to: ", output_dir)

#----8. Format p-value columns for export----

format_sci <- function(x, digits = 2) {
  ifelse(
    is.na(x),
    NA_character_,
    toupper(formatC(as.numeric(x), format = "e", digits = digits))
  )
}

format_pvalue_columns <- function(df, digits = 2) {
  
  pvalue_cols <- names(df)[
    stringr::str_detect(
      names(df),
      stringr::regex("P_value|p_value|pvalue|p\\.value|FDR|q_value|qvalue|padj|adjust", 
                     ignore_case = TRUE)
    )
  ]
  
  df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(pvalue_cols),
        ~ format_sci(.x, digits = digits)
      )
    )
}


#----9. Export Best and Expanded separately----

export_by_mapping <- function(df, file_stub) {
  
  df %>%
    dplyr::group_split(Mapping) %>%
    purrr::walk(function(sub_df) {
      
      mapping_name <- unique(sub_df$Mapping)
      
      export_df <- sub_df %>%
        dplyr::select(-Mapping) %>%
        format_pvalue_columns(digits = 2)
      
      readr::write_csv(
        export_df,
        file.path(
          output_dir,
          paste0("p005_", mapping_name, "_", file_stub, ".csv")
        )
      )
    })
}


export_by_mapping(
  kegg_annotated,
  "KEGG_annotated_all"
)

export_by_mapping(
  kegg_metabolism_only,
  "KEGG_metabolism_only_all"
)

export_by_mapping(
  kegg_metabolism_only_pathwayP005,
  "KEGG_metabolism_only_pathwayP005"
)

export_by_mapping(
  kegg_metabolism_only_pathwayFDR01,
  "KEGG_metabolism_only_pathwayFDR01"
)

#----10. Summary table----

summary_table <- kegg_annotated %>%
  dplyr::group_by(Mapping, Condition) %>%
  dplyr::summarise(
    Total_pathways = dplyr::n(),
    Metabolism_pathways = sum(Is_Metabolism, na.rm = TRUE),
    Metabolism_pathway_P005 = sum(Is_Metabolism & P_value < pathway_p_cutoff, na.rm = TRUE),
    Metabolism_pathway_FDR01 = sum(Is_Metabolism & FDR < pathway_fdr_cutoff, na.rm = TRUE),
    .groups = "drop"
  )

summary_table %>%
  dplyr::group_split(Mapping) %>%
  purrr::walk(function(sub_df) {
    
    mapping_name <- unique(sub_df$Mapping)
    
    export_df <- sub_df %>%
      dplyr::select(-Mapping)
    
    readr::write_csv(
      export_df,
      file.path(
        output_dir,
        paste0("p005_", mapping_name, "_KEGG_metabolism_summary.csv")
      )
    )
  })

print(summary_table)

message("Done. Output files saved to: ", output_dir)


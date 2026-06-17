#===================================================================
# Script: 6_Up-regulated Metabolites, Keep All Classified, Trim "Others" to N=70
#===================================================================

#---- 0. Load packages ----
library(tidyverse)
library(pheatmap)

#---- 1. Set paths ----
base_path <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data"

data_path <- file.path(base_path, "5_Metabolite Annotated")
output_path <- file.path(base_path, "6_HeatMap")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

#---- 2. Download KEGG Pathway Mappings ----
message("正在下載 KEGG Pathway 關聯資料庫...")

kegg_link_raw <- readLines("https://rest.kegg.jp/link/pathway/cpd", warn = FALSE)
kegg_link_df <- tibble(raw = kegg_link_raw) %>%
  filter(str_detect(raw, "\t")) %>%
  separate(raw, into = c("Compound", "Pathway"), sep = "\t", extra = "drop") %>%
  mutate(KEGG_ID = str_replace(Compound, "cpd:", ""), Pathway_ID = str_replace(Pathway, "path:map", ""))

kegg_path_raw <- readLines("https://rest.kegg.jp/list/pathway", warn = FALSE)
kegg_path_df <- tibble(raw = kegg_path_raw) %>%
  filter(str_detect(raw, "\t")) %>%
  separate(raw, into = c("Pathway", "Pathway_Name"), sep = "\t", extra = "drop") %>%
  mutate(Pathway_ID = str_replace(Pathway, "path:map", ""))

cpd_pathways <- kegg_link_df %>%
  left_join(kegg_path_df, by = "Pathway_ID") %>%
  filter(!is.na(Pathway_Name)) %>%
  group_by(KEGG_ID) %>%
  summarise(Pathways = paste(Pathway_Name, collapse = " | "), .groups = "drop")

pathway_dict <- setNames(cpd_pathways$Pathways, cpd_pathways$KEGG_ID)

#---- 3. Custom Function: 廣域生化分類大師 ----
classify_metabolite_broad <- function(kegg_str, name_str) {
  ids <- str_squish(unlist(str_split(kegg_str, ";")))
  paths <- paste(na.omit(pathway_dict[ids]), collapse = " | ")
  
  combined_text <- str_to_lower(paste(paths, name_str))
  
  if (str_detect(combined_text, "citrate cycle|tca cycle|tricarboxylic acid|citrate|isocitrate|succinate|fumarate|malate|oxaloacetate|ketoglutarate|aconitate")) {
    return("TCA cycle")
  } else if (str_detect(combined_text, "amino sugar|nucleotide sugar|udp-|gdp-|cdp-|adp-|cmp-|glucosamine|galactosamine|mannosamine|neuramin|sialic|hexosamine|n-acetylglucosamine|n-acetylgalactosamine|n-acetylmannosamine")) {
    return("Nucleotide-Sugar")
  } else if (str_detect(combined_text, "purine|pyrimidine|nucleotide metabolism|adenosine|guanosine|uridine|cytidine|inosine|xanthine|hypoxanthine|adenine|guanine|uracil|cytosine|thymine|atp|adp|amp|utp|udp|ump|gtp|gdp|gmp|ctp|cdp|cmp|nad\\+|nadp\\+|flavin|nicotinamide")) {
    return("Nucleotide")
  } else if (str_detect(combined_text, "amino acid|glutathione|gsh|gssg|arginine|proline|alanine|aspartate|glutamate|glycine|serine|threonine|cysteine|methionine|valine|leucine|isoleucine|lysine|histidine|phenylalanine|tyrosine|tryptophan|polyamines|putrescine|spermidine|spermine|peptide|ornithine|citrulline|homocysteine|creatine|taurine")) {
    return("Amino Acid metabolism")
  } else if (str_detect(combined_text, "lipid|fatty acid|steroid|sphing|glycero|arachid|linole|choline|ethanolamine|ceramide|carnitine|palmitoyl|oleoyl|stearoyl|acyl|prostaglandin|eicosanoid|cholesterol|bile acid")) {
    return("Lipid")
  } else if (str_detect(combined_text, "carbohydrate|glycolysis|pentose phosphate|galactose|starch|sucrose|fructose|mannose|glucose|ribose|lactose|maltose|trehalose|sorbitol|mannitol|glyceraldehyde|pyruvate|lactate|glycerate|gluconate|glucuronate|saccharide|glycan|inositol")) {
    return("Sugar")
  } else {
    return("others")
  }
}

#---- 4. Custom Function: 專業代謝物縮寫與去括號處理 ----
simplify_metabolite_names_pro <- function(names_vec) {
  primary_names <- sapply(names_vec, function(x) {
    parts <- str_squish(str_split(x, ";")[[1]])
    paste(head(parts, 2), collapse = " / ")
  })
  
  clean_names <- primary_names %>%
    str_replace_all("\\(.*?\\)", "") %>% 
    str_replace_all("(?i)\\bAdenosine 5'-triphosphate\\b", "ATP") %>%
    str_replace_all("(?i)\\bAdenosine 5'-diphosphate\\b", "ADP") %>%
    str_replace_all("(?i)\\bAdenosine 5'-monophosphate\\b", "AMP") %>%
    str_replace_all("(?i)\\bUridine 5'-triphosphate\\b", "UTP") %>%
    str_replace_all("(?i)\\bUridine 5'-diphosphate\\b", "UDP") %>%
    str_replace_all("(?i)\\bUridine 5'-monophosphate\\b", "UMP") %>%
    str_replace_all("(?i)\\bGuanosine 5'-triphosphate\\b", "GTP") %>%
    str_replace_all("(?i)\\bGuanosine 5'-diphosphate\\b", "GDP") %>%
    str_replace_all("(?i)\\bGuanosine 5'-monophosphate\\b", "GMP") %>%
    str_replace_all("(?i)\\bCytidine 5'-triphosphate\\b", "CTP") %>%
    str_replace_all("(?i)\\bCytidine 5'-diphosphate\\b", "CDP") %>%
    str_replace_all("(?i)\\bCytidine 5'-monophosphate\\b", "CMP") %>%
    str_replace_all("(?i)\\bNicotinamide adenine dinucleotide phosphate\\b", "NADP+") %>%
    str_replace_all("(?i)\\bNicotinamide adenine dinucleotide\\b", "NAD+") %>%
    str_replace_all("(?i)\\bCoenzyme A\\b", "CoA") %>%
    str_replace_all("(?i)Acetyl-CoA", "Ac-CoA") %>%
    str_replace_all("(?i)\\bUDP-glucose\\b", "UDP-Glc") %>%
    str_replace_all("(?i)UDP-N-acetyl-glucosamine", "UDP-GlcNAc") %>%
    str_replace_all("(?i)UDP-N-acetyl-galactosamine", "UDP-GalNAc") %>%
    str_replace_all("(?i)N-Acetyl-glucosamine", "GlcNAc") %>%
    str_replace_all("(?i)N-Acetyl-galactosamine", "GalNAc") %>%
    str_replace_all("(?i)N-Acetyl-mannosamine", "ManNAc") %>%     
    str_replace_all("(?i)N-Acetyl-neuraminate", "Neu5Ac") %>%
    str_replace_all("(?i)\\bPhosphatidylcholine\\b", "PC") %>%
    str_replace_all("(?i)\\bPhosphatidylethanolamine\\b", "PE") %>%
    str_replace_all("(?i)\\bPhosphatidylserine\\b", "PS") %>%
    str_replace_all("(?i)\\bPhosphatidylinositol\\b", "PI") %>%
    str_replace_all("(?i)\\bSphingomyelin\\b", "SM") %>%
    str_replace_all("(?i)\\bCeramide\\b", "Cer") %>%
    str_replace_all("(?i)\\btriphosphate\\b", "TP") %>%
    str_replace_all("(?i)\\bdiphosphate\\b", "DP") %>%
    str_replace_all("(?i)\\bphosphate\\b", "P") %>%
    str_replace_all("(?i)N-Acetyl-", "NAc-") %>%
    str_replace_all("(?i)\\bacid\\b", "") %>%
    str_replace_all("(?i)alpha-D-|beta-D-|alpha-|beta-|D-|L-|sn-", "") %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("-$", "") %>%
    str_replace_all(" / $", "") %>%
    str_squish()
  
  return(unname(clean_names))
}

#---- 5. 設定要處理的文件列表 ----
file_list <- list(
  POS = "1_POS_Final_Measurement_with_KEGG_Names.csv",
  NEG = "1_NEG_Final_Measurement_with_KEGG_Names.csv"
)

#---- 6. 主循環：分類保送、others 淘汰機制 ----
for (mode_name in names(file_list)) {
  
  file_name <- file_list[[mode_name]]
  full_data <- read_csv(file.path(data_path, file_name), show_col_types = FALSE)
  
  # 步驟 A: 基本清洗與全 0 過濾
  processed_df <- full_data %>%
    filter(!is.na(Primary_Name) & Primary_Name != "Unknown") %>%
    rename(
      `Ferroptosis-1` = `RSL-1`,
      `Ferroptosis-2` = `RSL-2`,
      `Ferroptosis-3` = `RSL-3`
    ) %>%
    filter((`DMSO-1` + `DMSO-2` + `DMSO-3` + `Ferroptosis-1` + `Ferroptosis-2` + `Ferroptosis-3`) > 0)
  
  # 步驟 B: 統計計算 (Log2FC > 1, p < 0.05)
  stats_df <- processed_df %>%
    mutate(across(c(`DMSO-1`, `DMSO-2`, `DMSO-3`, `Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`), ~ .x + 1)) %>%
    rowwise() %>%
    mutate(
      Mean_DMSO = mean(c(`DMSO-1`, `DMSO-2`, `DMSO-3`)),
      Mean_Ferro = mean(c(`Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`)),
      Log2FC = log2(Mean_Ferro / Mean_DMSO),
      P_value = tryCatch(
        t.test(c(`Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`), c(`DMSO-1`, `DMSO-2`, `DMSO-3`))$p.value, 
        error = function(e) NA
      )
    ) %>%
    ungroup() %>%
    filter(P_value < 0.05 & Log2FC > 1)
  
  if(nrow(stats_df) == 0) {
    next
  }
  
  # 步驟 C: 廣域分類與瘦身
  annotated_df <- stats_df %>%
    rowwise() %>%
    mutate(Category = classify_metabolite_broad(KEGG_ID, Primary_Name)) %>%
    ungroup() %>%
    mutate(Short_Name = simplify_metabolite_names_pro(Primary_Name)) %>%
    group_by(Short_Name) %>%
    mutate(Unique_Name = if(n() > 1) paste0(Short_Name, "_", row_number()) else Short_Name) %>%
    ungroup()
  
  # 💡 步驟 D: [精華] VIP 保送與 Others 淘汰機制 (修改為 70 個)
  # 1. 抓出所有「有分類」的名單 (保證一個都不砍)
  df_known <- annotated_df %>% filter(Category != "others")
  n_known <- nrow(df_known)
  
  # 2. 💡 總名額改成 70！計算總共還缺幾個名額才滿 70
  target_total <- 70
  n_others_needed <- max(0, target_total - n_known)
  
  # 3. 處理 others：按 Log2FC 由高到低排，只取缺額數量
  df_others <- annotated_df %>% 
    filter(Category == "others") %>% 
    arrange(desc(Log2FC)) %>% 
    head(n_others_needed)
  
  # 4. 把保送名單和倖存的 others 合併，並強制按類別排列
  final_df <- bind_rows(df_known, df_others) %>%
    arrange(Category, desc(Log2FC))
  
  # 步驟 E: 準備繪圖矩陣與側條註解
  mat <- final_df %>%
    select(Unique_Name, `DMSO-1`, `DMSO-2`, `DMSO-3`, `Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`) %>%
    column_to_rownames("Unique_Name") %>%
    as.matrix() %>%
    log2()
  
  annotation_row <- final_df %>%
    select(Unique_Name, Category) %>%
    column_to_rownames("Unique_Name")
  
  ann_colors <- list(
    Category = c(
      "Nucleotide" = "#E41A1C",
      "Nucleotide-Sugar" = "#377EB8",
      "Sugar" = "#4DAF4A",
      "TCA cycle" = "#984EA3",
      "Amino Acid metabolism" = "#FF7F00",
      "Lipid" = "#FFFF33",
      "others" = "#A65628"
    )
  )
  
  # 準備 Gaps 分隔線
  category_counts <- table(final_df$Category)[unique(final_df$Category)]
  gaps_pos <- cumsum(category_counts)
  gaps_pos <- gaps_pos[-length(gaps_pos)]
  
  # 步驟 F: 繪製完美分群 Heatmap (固定 cell 高低，不會被拉伸)
  pheatmap(
    mat,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = FALSE,        
    gaps_row = gaps_pos,         
    annotation_row = annotation_row, 
    annotation_colors = ann_colors,  
    show_rownames = TRUE,
    show_colnames = FALSE,       
    border_color = NA,           
    fontsize_row = 10,            
    cellwidth = 25,              
    cellheight = 10,             
    color = colorRampPalette(c("blue", "white", "red"))(100), 
    main = paste0("Up-regulated Metabolites Classified (", mode_name, "): Ferroptosis vs DMSO"),
    filename = file.path(output_path, paste0("6_Heatmap_", mode_name, "_Target70_SmartTrim.png"))
  )
  
  write_csv(
    final_df, 
    file.path(output_path, paste0("6_UpRegulated_Metabolites_", mode_name, "_Target70_SmartTrim.csv"))
  )
}

message("大哥/大姊，目標已上調至 70 個！一樣保送所有分類名單，剩下的再從 Others 去補，快去 6_HeatMap 資料夾看看 6_ 開頭的新檔案吧！")
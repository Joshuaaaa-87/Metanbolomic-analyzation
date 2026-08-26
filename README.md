# Metanbolomic-analyzation

本專案是一套以 **R / R Markdown** 撰寫的非標靶代謝體資料分析流程，主要用於比較不同細胞死亡條件與 DMSO 對照組之間的代謝差異。分析內容涵蓋原始資料整理、背景訊號過濾、代謝物 ID 整合、KEGG 註解、PCA、UMAP、熱圖，以及 KEGG pathway enrichment analysis。

> 本 repository 為研究中的分析腳本集合，不是可直接安裝的 R package。執行前需準備原始 QI 匯出資料、修改腳本內的本機路徑，並依照下方順序執行。

## 分析流程

```text
QI 匯出資料（POS / NEG）
        │
        ▼
原始資料整理與 5× Blank filter
        │
        ▼
ChemSpider CSID 對應 ChEBI / HMDB / PubChem
        │
        ▼
Mass error ±15 ppm + 每個 feature 最多 10 個候選註解
        │
        ▼
ChEBI / PubChem / HMDB 轉換為 KEGG Compound ID
        │
        ▼
產生 KEGG-mapped POS / NEG measurement tables
        │
        ├── PCA
        ├── UMAP
        ├── 代謝物名稱註解與 Heatmap
        └── KEGG pathway enrichment
                ├── Single KEGG（保守主分析）
                ├── Best KEGG per feature（探索分析）
                └── Expanded KEGG（敏感度分析 / QC）
```

## 實驗分組

每組包含 3 個 biological replicates：

| 欄位前綴 | 實驗條件 | 分析標籤 |
|---|---|---|
| `DMSO` | Vehicle control | control |
| `STS` | STS | intrinsic apoptosis |
| `Fas` | Fas | extrinsic apoptosis |
| `RSL` | RSL | ferroptosis |
| `Nec` | Nec | necroptosis |
| `Blank` | Blank sample | 背景訊號過濾用，不納入 PCA / UMAP |

measurement table 經整理後必須包含以下欄位：

```text
Compound,
DMSO-1, DMSO-2, DMSO-3,
STS-1, STS-2, STS-3,
Fas-1, Fas-2, Fas-3,
RSL-1, RSL-2, RSL-3,
Nec-1, Nec-2, Nec-3,
Blank
```

## Repository 內容

| 檔案 | 功能 |
|---|---|
| `20260605_Raw_data_Tidy.Rmd` | 讀取 QI 原始資料、整理 measurement table、執行 5× Blank filter、抓取 ChemSpider external IDs、套用 mass-error 與 top-10 candidate filters，並建立 KEGG conversion 所需的 ID 清單。 |
| `20260608_SCID&PubChem_to_KEGG.Rmd` | 將 ChEBI、PubChem 與 HMDB ID 轉為 KEGG Compound ID，整合 mapping 結果，再對回 POS / NEG measurement tables。 |
| `20260608-1_UMAP.Rmd` | 合併 POS 與 NEG features，執行 UMAP 並輸出座標、圖檔與 RDS model。 |
| `20260608-2_PCA.Rmd` | 對 KEGG-mapped features 執行 PCA，輸出 explained variance、scores、loadings、Top 50 loading features 與圖檔。 |
| `20260608-3_KEGG_pathway_enrichment_death_vs_DMSO_15ppm_top10.Rmd` | 僅使用 `Single_KEGG` features 進行保守版 pathway enrichment，建議作為主分析。 |
| `20260608-3B(Mutiple)_KEGG_pathway_enrichment_death_vs_DMSO_15ppm_top10.Rmd.R` | 納入 Single 與 Multiple KEGG features，但每個 feature 只保留一個 KEGG ID，作探索分析。 |
| `20260608-3C(Expend)_KEGG_pathway_enrichment_death_vs_DMSO_15ppm_top10.Rmd.R` | 將 Multiple KEGG IDs 全部展開，作敏感度分析或 QC；可能放大 pathway hits。 |
| `20260609-3D_KEGG_Data整理.Rmd.R` | 合併 Best 與 Expanded enrichment 結果，依 KEGG map ID 分類，並篩選 metabolism pathways。 |
| `202606010-5_Metabolite Annotated.Rmd.R` | 由 KEGG REST 取得 compound primary names，加入最終 POS / NEG measurement tables。 |
| `202606010-6_HeatMap.Rmd.R` | 比較 ferroptosis 與 DMSO，篩選顯著差異代謝物並繪製 POS / NEG heatmaps。 |
| `202606010-6_HeatMap_CategoryVersion.Rmd.R` | 將上升代謝物分為 TCA cycle、nucleotide、amino acid、lipid、sugar 等類別，並產生分類 heatmap。 |

## 系統需求

- R 4.x
- 建議使用 RStudio 逐段執行 R Markdown
- 可連線至 ChemSpider、KEGG REST API 與 MetaboAnalyst
- 原始資料與分析輸出需要足夠的本機儲存空間

### R packages

```r
install.packages(c(
  "tidyverse",
  "pheatmap",
  "umap",
  "rvest",
  "matrixStats",
  "FactoMineR",
  "factoextra",
  "ggcorrplot",
  "httr2",
  "usethis",
  "webchem",
  "BiocManager"
))

BiocManager::install("Rdisop")
```

實際使用的主要功能如下：

- `tidyverse`：資料整理、統計與作圖
- `rvest`、`httr2`：讀取 ChemSpider 網頁資料
- `umap`：UMAP dimensionality reduction
- `factoextra`：PCA 輔助分析
- `pheatmap`：heatmap 繪製
- `Rdisop`、`matrixStats`：代謝體資料處理輔助功能

## 輸入資料

請將 QI 匯出的四個 CSV 檔放在同一個資料夾：

```text
0_QI report/
├── POS_Compound Identification.csv
├── NEG_Compound Identification.csv
├── POS_Compound Measurement.csv
└── NEG_Compound Measurement.csv
```

原始資料整理腳本目前假設：

1. measurement CSV 的前兩列不是正式 feature data，因此會被移除。
2. `Compound` 位於第 1 欄，樣本 intensity 位於第 14–29 欄。
3. identification table 至少含有 `Compound`、`Compound ID`、`CSID`、`Mass Error (ppm)`、`Score`、`Fragmentation Score` 與 `Isotope Similarity` 等欄位。
4. identification 與 measurement table 可透過 `Compound` 欄位連接。

若 QI 匯出的版型不同，請先修改 `20260605_Raw_data_Tidy.Rmd` 中的欄位選取與命名設定。

## 使用方式

### 1. 下載 repository

```bash
git clone https://github.com/Joshuaaaa-87/Metanbolomic-analyzation.git
cd Metanbolomic-analyzation
```

### 2. 修改資料路徑

腳本目前使用作者電腦上的絕對路徑，例如：

```r
base_path <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data"
```

執行前請搜尋並修改各檔案中的 `base_path`、`project_path`、`raw_path`、`raw_data_path`、`input_path` 與 `output_path`。

建議使用以下目錄結構：

```text
project/
├── 0_Raw Data/
│   ├── 0_QI report/
│   └── 1_Data(reanalysis)_Celldeath_Sup_Each_Filtered/
├── 2_PCA_KEGG_mapped_15ppm_top10/
├── 3_UMAP_Filtered_Measurement_POS_NEG/
├── 4_KEGG_Enrichment_15ppm_top10/
├── 4_KEGG_Enrichment_BestKEGG_15ppm_top10/
├── 4_KEGG_Enrichment_ExpandedKEGG_15ppm_top10/
├── 5_Metabolite Annotated/
└── 6_HeatMap/
```

### 3. 整理原始資料

在 RStudio 中開啟並由上到下執行：

```text
20260605_Raw_data_Tidy.Rmd
```

主要處理步驟：

1. 計算各組平均 intensity。
2. 當任一實驗條件的平均 intensity 大於 `5 × Blank` 時保留該 feature。
3. 以 measurement 中保留的 `Compound` 過濾 identification table。
4. 使用 CSID 讀取 ChemSpider 頁面中的 ChEBI、HMDB 與 PubChem ID。
5. 僅保留 mass error 位於 ±15 ppm 的註解。
6. 依序參考 mass error、fragmentation score、score、isotope similarity 與 external ID 完整度，每個 `Mode + Compound` 最多保留 10 個 candidate rows。
7. 每個 `Mode + Compound + ID type` 最多保留 10 個 external IDs。

ChemSpider 批次查詢預設每筆暫停 4–8 秒，並每 50 筆寫入 checkpoint。資料量大時可能需要數小時，請勿移除延遲，以免對外部服務送出過量請求。

### 4. 將 external IDs 轉為 KEGG IDs

接續同一個 R session，執行：

```text
20260608_SCID&PubChem_to_KEGG.Rmd
```

此腳本會：

- 透過 KEGG REST 將 ChEBI 與 PubChem ID 轉為 KEGG Compound ID。
- 匯出 HMDB ID 清單供 MetaboAnalyst 使用。
- 合併 ChEBI、PubChem、HMDB direct mapping 與 rescue mapping。
- 產生 `Single_KEGG` / `Multiple_KEGG` 狀態。
- 將 KEGG annotation 合併回 measurement tables。

HMDB 轉換需要手動完成：

1. 將 `11_HMDB_ID_list_for_KEGG_conversion_15ppm_top10.csv` 上傳至 [MetaboAnalyst ID Conversion](https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml)。
2. Original ID type 選擇 **HMDB ID**。
3. Target ID type 選擇 **KEGG ID**。
4. 將結果下載並命名為 `14_HMDB_to_KEGG_MetaboAnalyst_raw.csv`。
5. 把檔案放回 `output_path` 後，再繼續執行後續區塊。

後續分析的核心輸入為：

```text
17_NEG_measurement_with_KEGG_mapped_only_15ppm_top10.csv
17_POS_measurement_with_KEGG_mapped_only_15ppm_top10.csv
```

### 5. 執行 PCA 與 UMAP

```text
20260608-1_UMAP.Rmd
20260608-2_PCA.Rmd
```

UMAP 流程：

- 使用通過 5× Blank filter 的 POS + NEG features。
- 執行 `log2(x + 1)`、移除 zero-variance features，再進行 scaling。
- 參數為 `n_neighbors = 3`、`min_dist = 0.15`、Euclidean distance、seed = 1。
- 輸出 UMAP score、PNG 圖與 RDS object。

PCA 流程：

- 使用 KEGG-mapped POS + NEG features。
- 執行 `log2(x + 1)` 並移除 zero-variance features。
- 使用 `prcomp(center = TRUE, scale. = TRUE)`。
- 輸出 PC1/PC2 與 PC2/PC3 圖、scores、explained variance、loadings、Top 50 loading features 與 RDS object。

### 6. 執行 KEGG pathway enrichment

三種版本可依分析目的分別執行：

| 版本 | 使用資料 | 建議用途 |
|---|---|---|
| Single KEGG | 只保留一個 feature 對應一個 KEGG ID 的資料 | 最保守，建議作為正式主分析 |
| Best KEGG | Single + Multiple features；Multiple feature 僅取一個 KEGG ID | 擴大 coverage 的探索分析 |
| Expanded KEGG | 將每個 Multiple feature 的 KEGG IDs 全部展開 | 敏感度分析或 QC，不建議單獨作為主要結論 |

對每個 metabolite feature，腳本分別比較 STS、Fas、RSL、Nec 與 DMSO：

- `log2FC = log2((Condition mean + 1) / (DMSO mean + 1))`
- 使用 two-sample t-test 計算 p-value
- 使用 Benjamini–Hochberg 方法校正 FDR
- 探索性門檻：`log2FC >= 1` 且 `p < 0.05`
- 較嚴格門檻：`log2FC >= 1` 且 `FDR < 0.1`

pathway over-representation analysis 使用 hypergeometric test，並以本實驗中實際偵測且成功對應 KEGG 的 metabolites 作為 background，而非使用整個 KEGG database。

### 7. 整理 enrichment 結果

執行：

```text
20260609-3D_KEGG_Data整理.Rmd.R
```

此步驟會合併 Best 與 Expanded KEGG enrichment 結果，依 KEGG map ID 指派 top-level category，並輸出 metabolism-only pathways、`P < 0.05`、`FDR < 0.1` 與 summary tables。

### 8. 取得代謝物名稱並繪製 heatmap

依序執行：

```text
202606010-5_Metabolite Annotated.Rmd.R
202606010-6_HeatMap.Rmd.R
202606010-6_HeatMap_CategoryVersion.Rmd.R
```

一般 heatmap 比較 ferroptosis 與 DMSO，保留 `p < 0.05` 且 `|log2FC| > 1` 的 metabolites。分類版只保留上升 metabolites（`log2FC > 1`），優先保留可歸入已知生化類別的項目，再以 log2FC 最高的 `others` 補至目標 70 個。

## 主要輸出

完整流程會產生大量中間 QC tables。最重要的結果包括：

- 5× Blank filtered POS / NEG measurement 與 identification tables
- ChEBI、HMDB、PubChem 與 KEGG mapping tables
- 最終 KEGG-mapped POS / NEG measurement tables
- PCA scores、loadings、Top 50 features、PNG 與 RDS
- UMAP scores、sample annotation、PNG 與 RDS
- differential metabolite tables
- high-up metabolite lists
- KEGG compound-to-pathway mappings
- 各條件 pathway enrichment、Top 20 pathways 與 concise tables
- metabolism-only KEGG pathway tables
- POS / NEG significant-metabolite heatmaps

## 分析注意事項

1. 本 repository 未包含原始 QI data 或已產生的分析結果，因此無法只下載程式碼就重現完整分析。
2. 多個腳本使用作者電腦的絕對路徑，執行前必須修改。
3. 部分 R Markdown 區塊依賴前一步建立在 R workspace 中的 objects，建議依序在同一個 R session 執行。
4. ChemSpider、KEGG 與 MetaboAnalyst 的服務內容、回傳格式或 rate limit 可能改變，若查詢失敗需檢查外部服務狀態。
5. `Multiple_KEGG` 代表註解歧義，不代表一個 LC-MS feature 已被確認為多個代謝物。
6. Best KEGG 版本可提升 coverage，但所選 ID 不一定是真實 metabolite identity；正式解讀時應與 Single KEGG 結果並列比較。
7. Expanded KEGG 版本可能重複計算同一 feature 的 intensity 訊號，使 pathway hit 膨脹，應限於 QC 或敏感度分析。
8. 每組目前為 `n = 3`。UMAP ellipse 僅供視覺化，不應解讀為正式信賴區間；逐 feature t-test 的統計力亦有限。
9. 熱圖中的 broad biochemical categories 同時依賴 KEGG pathway 名稱與關鍵字規則，屬於視覺化用的廣義分類。
10. 此專案目前未提供 package lockfile、session information 或自動化測試；在正式發表前建議記錄 R 與 package versions，並逐步核對中間輸出。

## 資料庫與線上資源

- [KEGG](https://www.kegg.jp/)
- [KEGG REST API](https://rest.kegg.jp/)
- [ChemSpider](https://www.chemspider.com/)
- [HMDB](https://www.hmdb.ca/)
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- [ChEBI](https://www.ebi.ac.uk/chebi/)
- [MetaboAnalyst](https://www.metaboanalyst.ca/)

## License

本 repository 目前未提供 `LICENSE` 檔案。若計畫修改、散布或用於其他研究，請先向 repository 作者確認授權方式。

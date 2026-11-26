# ==============================================================================
# 文件说明：PCN（胰腺囊性肿瘤）预测指标 binary RR 数据清洗脚本
# 作用：读取 bprf_data/ 目录下各指标的 xlsx 文件，标准化清洗后
#       输出供 binary risk.r（MR-BRT 建模）使用的 CSV 文件
# 输入：/home/hopehappy/PCN_BPRF/bprf_data2/*.csv
# 输出：{out_dir}/{indicator}_all_extraction.csv（每个指标一个文件）
#
# 数据列说明：
#   pmid           - PubMed ID，唯一研究标识（等同于原流程的 nid）
#   or_rr          - 效应量（OR 或 RR）
#   lower/upper    - 95% CI
#   cov_adj        - 是否调整协变量（0=未调整，1=已调整）
#   cv_outcome     - 结局测量偏倚（0=低偏倚，1=高偏倚）
#                    注：部分文件中该列含非标准编码（行号/文本），已统一截断为 0/1
#   cv_studydesign - 研究设计偏倚（1=回顾性/较差设计）
#   cv_multicentre - 是否多中心（1=多中心研究）
#   cv_aging       - 老龄化相关偏倚
#   cv_mri/ct/mixmod/eus - 影像学模式（仅影像类指标有此列，其余填 0）
#
# 特殊处理说明：
#   - row_type: 仅保留 "standard" 和 "positive"，排除有序分类/合并类别等
#   - cv_outcome: 非数值或 >1 的值均截断为 1（保守假设存在偏倚）
#   - 空 xlsx 文件（0数据行）自动跳过
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 加载依赖包
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(tidyverse)

# ==============================================================================
# 第一部分：全局参数
# ==============================================================================

data_dir <- "/home/hopehappy/PCN_BPRF/bprf_data2"    # 原始数据目录
out_dir  <- "/home/hopehappy/PCN_BPRF/cleaned_data"  # 清洗结果输出目录

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 保留的 row_type 类型：
#   standard    - 标准二元比较（指标阳性 vs 阴性）
#   positive    - 正向配对（部分连续变量指标的阳性组，如 MPD_Dilation）
# 排除：ordered_category（有序分类分析）、negative_pair（反向配对）、merged_*（合并阈值分析）
valid_row_types <- c("standard", "positive")

# 影像学协变量列名（仅影像类指标存在；临床/血清类指标缺失此列，后续填 0）
# cv_mri    : 仅使用 MRI 评估
# cv_ct     : 仅使用 CT 评估
# cv_mixmod : 使用混合影像模态
# cv_eus    : 仅使用超声内镜（EUS）
imaging_cv_cols <- c("cv_mri", "cv_ct", "cv_mixed", "cv_eus")

# 所有指标均含的基础偏倚协变量
base_cv_cols <- c(
  "cv_studydesign",
  "cv_multicentre",
  "cv_aging",
  "cv_reader_quality",
  "cv_adj",
  "cv_population",
  "cv_surgical_criteria_reported"
)

new_cv_cols <- c(
  "cv_disease_IPMN",
  "cv_disease_MCN",
  "cv_subgroup_imaging",
  "cv_subgroup_pathology"
)

ipmn_pmids <- c(
  "25822687", "26055134", "26689209", "26905937", "26942155",
  "26959611", "27283593", "27607098", "27829085", "27882420",
  "28030856", "28079542", "28088598", "28122068", "28144718",
  "28593459", "28697135", "29258117", "29730906", "29961172",
  "30170819", "30278929", "30361110", "30591368", "30672797",
  "31117111", "31388508", "31478800", "31611130", "31745906",
  "31798784", "32008558", "32407968", "32584315", "32655269",
  "33114028", "33145018", "33250091", "33385536", "33409798",
  "33578680", "33750453", "33767983", "34009303", "34112528",
  "34474037", "35275259", "35703792", "35781178", "36029231",
  "36197166", "36402960", "36675810", "36716959", "36761937",
  "36881188", "37078135", "37434479", "39151506", "39418295",
  "39258374", "40620507", "41409549", "41820087", "41392225",
  "42878968", "41720652", "42001163"
)

mcn_pmids <- c(
  "27760255", "28936758", "32011526", "32187927",
  "36822873", "37346158", "41905867"
)

mixed_disease_pmids <- c("40856421")

subgroup_pmids <- c(
  "25822687", "26055134", "26689209", "26905937",
  "27607098", "27829085", "28079542", "28122068",
  "28088598", "28593459", "29730906", "29961172",
  "30591368", "31388508", "31478800", "31611130",
  "31745906", "32584315", "32655269", "33145018",
  "33385536", "33578680", "34009303", "35275259",
  "36197166", "36402960", "36675810", "37078135",
  "37434479", "40620507"
)

imaging_definition_pmids <- c(
  "28088598", "31478800", "40620507", "28079542",
  "36402960", "29961172", "33145018", "37434479"
)

# 后续人工复核确认的 PMID 限定：
# 1) 以下 PMID 在所有指标中全局剔除；
# 2) PMID 34009303 仅允许用于复合 HRS 指标。
manual_excluded_pmids <- c("33915626", "34059873", "42277584")
hrs_only_pmid <- "34009303"
hrs_indicator_keys <- c(
  "hrs",
  "hrs_ref0",
  "high_risk_stigmata",
  "high_risk_stigmata_present",
  "high_risk_stigmata_hrs_present"
)

# ==============================================================================
# 第二部分：辅助函数
# ==============================================================================

# 将协变量列安全截断为 0/1 整型
# 处理规则：
#   - 数值 0    → 0
#   - 数值 > 0  → 1（含行号、评分等非标准编码，保守赋值为"存在偏倚"）
#   - 非数值/NA → 1（保守假设存在偏倚）
to_binary_cv <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))   # 非数值（如中文文本）强制转为 NA
  result <- ifelse(is.na(x_num), 1L,         # NA → 1（保守）
             ifelse(x_num == 0, 0L, 1L))     # 0 → 0，其余 → 1
  return(as.integer(result))
}

# ==============================================================================
# 第三部分：扫描并读取所有 xlsx 文件
# ==============================================================================

csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
message("发现 csv 文件数：", length(csv_files))

dt_list <- lapply(csv_files, function(f) {
  df <- tryCatch(
    fread(file = f, colClasses = "character", encoding = "UTF-8"),
    error = function(e) {
      message("读取失败：", basename(f), " - ", e$message)
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) {
    message("跳过空文件：", basename(f))
    return(NULL)
  }
  df <- as.data.table(df)
  df[, source_indicator := tools::file_path_sans_ext(basename(f))]
  return(df)
})

non_empty <- !sapply(dt_list, is.null)
message("有数据的文件数：", sum(non_empty))

# 纵向合并（fill=TRUE 处理各指标列差异）
dt_all <- rbindlist(dt_list[non_empty], use.names = TRUE, fill = TRUE)
message("合并后总行数：", nrow(dt_all))

# ==============================================================================
# 第四部分：基础过滤
# ==============================================================================

# 仅保留标准二元比较行
n_before <- nrow(dt_all)
dt_all <- dt_all[row_type %in% valid_row_types]
message("row_type 过滤后：", nrow(dt_all), " 行（删除 ", n_before - nrow(dt_all), " 行）")

# ==============================================================================
# 第五部分：数值列类型转换
# （因读取时全部为文本，此处统一转换）
# ==============================================================================

num_cols <- c("or_rr", "lower", "upper", "n_eh", "n_el", "n_event", "n", "year", "cov_adj")

for (col in intersect(num_cols, names(dt_all))) {
  dt_all[, (col) := suppressWarnings(as.numeric(get(col)))]
}

dt_all[, pmid := suppressWarnings(as.integer(pmid))]

# ==============================================================================
# 第六部分：效应量有效性过滤
# ==============================================================================

n_before <- nrow(dt_all)
dt_all <- dt_all[!is.na(or_rr) & !is.na(lower) & !is.na(upper)]
dt_all <- dt_all[or_rr > 0 & lower > 0 & upper > 0]   # 对数变换要求正数
message("效应量过滤后：", nrow(dt_all), " 行（删除 ", n_before - nrow(dt_all), " 行）")

# ==============================================================================
# 第七部分：对数变换效应量
# MR-BRT 在对数空间建模，以满足正态性假设
# ==============================================================================

dt_all[, ln_effect := log(or_rr)]
# 用 95% CI 反推对数空间标准误：SE_log ≈ (log(upper) - log(lower)) / 3.92
dt_all[, ln_se := (log(upper) - log(lower)) / 3.92]

# 诊断
cat("\n=== ln_se 分布 ===\n"); print(summary(dt_all$ln_se))
cat("=== ln_effect 分布 ===\n"); print(summary(dt_all$ln_effect))

# 删除标准误异常行
n_before <- nrow(dt_all)
dt_all <- dt_all[!is.na(ln_se) & ln_se > 0]
if (n_before > nrow(dt_all))
  message("删除 SE 异常行：", n_before - nrow(dt_all), " 行")

# ==============================================================================
# 第八部分：偏倚风险协变量处理
# ==============================================================================

# --- 基础协变量（所有指标均有）---

# cv_outcome：部分文件含非标准编码（行号/中文文本），统一截断为 0/1
# Lymphadenopathy: 值为 1/3/5/7/8/10/12/13/14/15 → 均截断为 1（存在偏倚）
# MPD_Dilation:    部分行值为 "定义不一致"/"缺乏定义数据" → NA → 截断为 1
# Standardize renamed/added covariates before binary conversion.
if ("cv_mixmod" %in% names(dt_all) && !"cv_mixed" %in% names(dt_all)) {
  dt_all[, cv_mixed := cv_mixmod]
}

if ("cv_multicenter" %in% names(dt_all) && !"cv_multicentre" %in% names(dt_all)) {
  dt_all[, cv_multicentre := cv_multicenter]
}

if (!"cov_adj" %in% names(dt_all)) {
  dt_all[, cov_adj := 0]
}
dt_all[, cov_adj := as.integer(replace(cov_adj, is.na(cov_adj), 0L))]
dt_all[, cv_adj := cov_adj]

if ("cv_outcome" %in% names(dt_all)) {
  dt_all[, cv_outcome := NULL]
}

for (col in base_cv_cols) {
  if (col %in% names(dt_all)) {
    dt_all[, (col) := to_binary_cv(get(col))]
  } else {
    dt_all[, (col) := 1L]   # 缺列：保守赋值为存在偏倚
  }
}

# --- 影像学协变量（仅影像类指标有）---
# 临床/血清类指标（如 CA19-9、CEA、Diabetes）无此列，缺失填 0（不适用该偏倚来源）
for (col in imaging_cv_cols) {
  if (col %in% names(dt_all)) {
    dt_all[, (col) := to_binary_cv(get(col))]
  } else {
    dt_all[, (col) := 0L]
  }
}

dt_all[, pmid := trimws(as.character(pmid))]

# --- 后续 PMID 人工复核规则 ---
# 先完成全局剔除，再进行疾病范围与亚组协变量赋值。
dt_all <- dt_all[!pmid %in% manual_excluded_pmids]

# PMID 34009303 只保留复合 HRS 指标的记录。
# 指标键仅用于此 PMID 限定，不改变其他指标的清洗逻辑。
indicator_key <- tolower(trimws(as.character(dt_all$source_indicator)))
indicator_key <- gsub("[^a-z0-9]+", "_", indicator_key)
indicator_key <- gsub("^_+|_+$", "", indicator_key)
dt_all <- dt_all[
  pmid != hrs_only_pmid | indicator_key %in% hrs_indicator_keys
]

matched_disease_pmid <- dt_all$pmid %in% c(ipmn_pmids, mcn_pmids, mixed_disease_pmids)
unmatched_disease_pmids <- sort(unique(dt_all$pmid[!matched_disease_pmid]))

if (length(unmatched_disease_pmids) > 0) {
  stop(
    "Unmatched PMID(s) for disease-scope covariates: ",
    paste(unmatched_disease_pmids, collapse = ", ")
  )
}

dt_all[, cv_disease_IPMN := fifelse(pmid %in% ipmn_pmids, 1L, 0L)]
dt_all[, cv_disease_MCN := fifelse(pmid %in% mcn_pmids, 1L, 0L)]

pmid_chr <- trimws(as.character(dt_all$pmid))
in_subgroup <- pmid_chr %in% subgroup_pmids
in_imaging_definition <- pmid_chr %in% imaging_definition_pmids

dt_all[, cv_subgroup_imaging := as.integer(in_subgroup & in_imaging_definition)]
dt_all[, cv_subgroup_pathology := as.integer(in_subgroup & !in_imaging_definition)]

for (col in new_cv_cols) {
  dt_all[is.na(get(col)), (col) := 0L]
  dt_all[, (col) := as.integer(get(col))]
}

all_cv_cols <- c(base_cv_cols, imaging_cv_cols, new_cv_cols)

# ==============================================================================
# 第九部分：协变量调整程度哑变量
# cov_adj 为二元（0=未调整任何协变量，1=至少调整了一个）
# cv_adj_L0：1=完全未调整（最差质量），0=有调整
# ==============================================================================

dt_all[, cv_adj_L0 := ifelse(cov_adj == 0L, 1L, 0L)]

# ==============================================================================
# 第十部分：添加 MR-BRT 必需的参考/暴露标识
# ==============================================================================

dt_all[, ref := 0L]   # 参考组（指标阴性）
dt_all[, alt := 1L]   # 暴露组（指标阳性）

# ==============================================================================
# 第十一部分：最终列选取与分指标保存
# ==============================================================================

final_cols <- c(
  # --- 研究标识 ---
  "pmid",                      # 研究 ID（对应原流程 nid）
  "author", "country", "year", # 研究元信息
  # --- 指标描述 ---
  "indicators",                # 指标全称（含阈值定义）
  "type_ind",                  # 指标类型（Imaging / Clinical / Tumor marker）
  "major_group",               # 大类分组（01_Mural_Nodule 等）
  "imaging",                   # 影像学检查方式
  # --- 效应量 ---
  "or_rr", "lower", "upper",  # 原始 OR/RR 及 95% CI
  "or_source",                 # 效应量来源（reported / calculated_from_2x2）
  "ln_effect", "ln_se",        # 对数效应量及标准误（MR-BRT 主要输入）
  # --- 样本量 ---
  "n_eh", "n_el",              # 高/低风险组事件数
  "n_event", "n",              # 总事件数及总样本量
  # --- 协变量 ---
  "cov_adj",                   # 原始调整协变量数（0/1）
  all_cv_cols,                 # 所有偏倚风险协变量（cv_ 前缀）
  "cv_adj_L0",                 # 调整程度哑变量
  # --- MR-BRT 标识 ---
  "ref", "alt",                # binary 分析的参考/暴露标识
  # --- 元数据 ---
  "row_type",                  # 行类型（standard / positive）
  "source_indicator"           # 来源文件名（即指标名称）
)

final_cols <- intersect(final_cols, names(dt_all))
dt_out <- dt_all[, ..final_cols]

# 诊断 cv_ 列分布
cat("\n=== 各 cv_ 协变量分布（均值，越高说明该偏倚越普遍）===\n")
cv_summary <- dt_out[, lapply(.SD, mean, na.rm = TRUE), .SDcols = all_cv_cols]
print(cv_summary)

# --- 分指标保存 ---
cat("\n=== 各指标清洗后行数 ===\n")
indicators <- dt_out[, unique(source_indicator)]

for (ind in sort(indicators)) {
  dt_ind <- dt_out[source_indicator == ind]

  # 诊断：无方差的 cv_ 列（CovFinder 会自动剔除，此处仅提示）
  cv_present <- intersect(all_cv_cols, names(dt_ind))
  zero_var <- cv_present[sapply(cv_present, function(col) {
    vals <- dt_ind[[col]]
    length(unique(vals[!is.na(vals)])) <= 1
  })]
  if (length(zero_var) > 0)
    message("  [提示] ", ind, " 无方差 cv_ 列（binary risk.r 中将自动剔除）: ",
            paste(zero_var, collapse = ", "))

  out_path <- file.path(out_dir, paste0(ind, "_all_extraction.csv"))
  write.csv(dt_ind, out_path, row.names = FALSE)
  cat(sprintf("  %-40s → %2d 行 → %s\n", ind, nrow(dt_ind), basename(out_path)))
}

message("\n清洗完成。共 ", nrow(dt_out), " 行，", length(indicators), " 个指标。")
message("输出目录：", out_dir)

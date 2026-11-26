library(readxl)
library(dplyr)
library(stringr)

read_sheet_to_prep_input <- function(filepath, sheet_name) {
  
  # 1. 原始读取，保留前两列
  raw <- readxl::read_excel(filepath, sheet = sheet_name, col_names = FALSE, .name_repair = "minimal")
  raw <- raw[, 1:2]
  colnames(raw) <- c("col1", "col2")
  
  # 提前获取原表头，判断是否为 Cumulative（累计概率）
  y_header <- as.character(raw$col2[1])
  is_cumulative <- str_detect(str_to_lower(y_header), "cum|incidence|risk")
  
  # 2. 安全转换为数字，划分区块
  raw <- raw %>%
    mutate(
      num1 = suppressWarnings(as.numeric(col1)),
      num2 = suppressWarnings(as.numeric(col2)),
      is_num = !is.na(num1) & !is.na(num2) 
    ) %>% 
    mutate(group_id = cumsum(!is_num))
  
  num_data <- raw %>% filter(is_num)
  
  if (nrow(num_data) == 0) {
    stop(paste("在", sheet_name, "中未检测到任何有效的数字数据！"))
  }
  
  valid_groups <- unique(num_data$group_id)
  
  # ---------------------------------------------------------
  # 3. 提取坐标（永远是第一块纯数字区）
  # ---------------------------------------------------------
  group_coords <- valid_groups[1]
  coords_raw <- num_data %>% filter(group_id == group_coords) %>% mutate(time = num1, y_val = num2)
  
  # 自动处理百分数与生存率转换
  if (max(coords_raw$y_val, na.rm = TRUE) > 1) {
    coords_raw$y_val <- coords_raw$y_val / 100
  }
  if (is_cumulative) {
    coords_raw$survival <- 1 - coords_raw$y_val
  } else {
    coords_raw$survival <- coords_raw$y_val
  }
  
  # 清洗单调性与补齐起点
  df_coords <- coords_raw %>% select(time, survival) %>% arrange(time)
  if (min(df_coords$time, na.rm = TRUE) > 0) {
    df_coords <- bind_rows(data.frame(time = 0, survival = 1), df_coords) %>% arrange(time)
  }
  df_coords$survival[1] <- 1
  df_coords$survival <- cummin(df_coords$survival) 
  
  # ---------------------------------------------------------
  # 4. 提取风险人数表（灵活兼容机制）
  # ---------------------------------------------------------
  # 如果包含2块及以上的数字区，最后一块就是风险表
  if (length(valid_groups) >= 2) {
    group_risk <- valid_groups[length(valid_groups)]
    df_risk <- num_data %>% 
      filter(group_id == group_risk) %>%
      select(trisk = num1, nrisk = num2) %>%
      arrange(trisk)
    
    # 截断对齐
    max_coord_time <- max(df_coords$time, na.rm = TRUE)
    df_risk <- df_risk %>% filter(trisk <= max_coord_time)
    
    totalpts_val <- df_risk$nrisk[1]
    
  } else {
    # 如果只有1块数字区（缺失风险表），生成空表
    df_risk <- data.frame(trisk = numeric(), nrisk = numeric())
    totalpts_val <- NA # 标记为缺失
  }
  
  # ---------------------------------------------------------
  # 5. 返回结果
  # ---------------------------------------------------------
  return(list(
    sheet = sheet_name,
    dat = as.data.frame(df_coords), 
    trisk = df_risk$trisk,          
    nrisk = df_risk$nrisk,          
    totalpts = totalpts_val     
  ))
}
###主分析
library(IPDfromKM)

filepath <- "Data Include_new.xlsx"
all_sheets <- readxl::excel_sheets(filepath)
target_sheets <- all_sheets[-1] # 跳过第一个说明 sheet

prep_results_list <- list()

for (sheet_name in target_sheets) {
  cat("----------------------------------\n")
  cat("正在读取:", sheet_name, "...\n")
  
  # 1. 提取数据
  input_data <- read_sheet_to_prep_input(filepath, sheet_name)
  
  # 2. 检查是否缺失风险表 (即 totalpts 是否为 NA)
  if (is.na(input_data$totalpts)) {
    warning(paste("跳过 preprocess 预处理：", sheet_name, "缺少风险人数表(N_risk)，算法无法获取初始人数(totalpts)。"))
    next # 跳过这一个，继续跑下一个 Sheet
  }
  
  # 3. 数据完整，进行预处理
  cat("数据提取成功，开始 preprocess...\n")
  prep_res <- preprocess(
    dat = input_data$dat, 
    trisk = input_data$trisk, 
    nrisk = input_data$nrisk, 
    totalpts = input_data$totalpts, 
    maxy = 1
  )
  
  prep_results_list[[sheet_name]] <- prep_res
  cat(sheet_name, "预处理完成！\n")
}

cat("\n==================================\n")
cat("批量跑批结束！成功处理", length(prep_results_list), "条曲线。\n")



# -------------------------------------------------------------------------
# 1. 加载必要的 R 包 (如果尚未加载)
# -------------------------------------------------------------------------
basic_info <- data.frame(
  sheet_name = c(
    "37530526",
    "31449158",
    "35511697",
    "39253809",
    "27390303_1cm",
    "27390303_2cm",
    "27390303_3cm",
    "39581796_1cm",
    "39581796_2cm",
    "39581796_3cm",
    "41188169",
    "33309223_小于15mm",
    "33309223_大于15mm",
    "39793956_1cm",
    "39793956_2_3cm",
    "42852691"
  ),
  study_id = c(
    "37530526",
    "31449158",
    "35511697",
    "39253809",
    "27390303",
    "27390303",
    "27390303",
    "39581796",
    "39581796",
    "39581796",
    "41188169",
    "33309223",
    "33309223",
    "39793956",
    "39793956",
    "42852691"
  ),
  analysis_set = c(
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "primary",
    "sensitivity",
    "sensitivity",
    "sensitivity",
    "sensitivity",
    "sensitivity",
    "primary"
  ),
  size_group = c(
    "Overall",
    "Overall",
    "Overall",
    "Overall",
    "<1 cm",
    "1-2 cm",
    "2-3 cm",
    "<1 cm",
    "1-2 cm",
    "2-3 cm",
    "Overall",
    "<15 mm",
    ">15 mm",
    "<1 cm",
    "2-3 cm",
    "Overall"
  ),
  stringsAsFactors = FALSE
)

library(IPDfromKM)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(broom)

# 设定输出文件夹，确保它们存在
out_dir <- "IPDfromKM_WF_HRS_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "01_IPDfromKM_native_plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "03_survival_curves"), showWarnings = FALSE)

# -------------------------------------------------------------------------
# 2. 批量个体数据重建 (getIPD) 与合并
# -------------------------------------------------------------------------
cat("\n=== 开始批量运行 getIPD 算法重建个体数据 ===\n")

ipd_results_list <- list()

# 遍历已经成功通过 preprocess 的结果
for (sheet_name in names(prep_results_list)) {
  
  cat("正在重建:", sheet_name, "...\n")
  
  prep_data <- prep_results_list[[sheet_name]]
  
  # 获取该 sheet 的元数据 (从你事先准备好的 basic_info 表中匹配)
  meta <- basic_info %>% filter(sheet_name == !!sheet_name)
  
  # 运行改良版迭代 KM 算法 (这里假设你没有确切的总事件数，设为 NULL 盲估)
  # 注意：如果你在 basic_info 里有每个组的确切事件数，请替换 NULL 为对应变量
  est <- tryCatch({
    getIPD(prep = prep_data, armind = 1, tot.events = NULL)
  }, error = function(e) {
    warning(paste("getIPD 失败:", sheet_name, "-", e$message))
    return(NULL)
  })
  
  if (is.null(est)) next
  
  # 提取生成的 IPD 数据
  ipd_df <- as.data.frame(est$IPD)
  # 统一列名以防 IPDfromKM 版本差异
  time_col <- intersect(c("time", "Time", "t", "T"), names(ipd_df))[1]
  status_col <- intersect(c("status", "Status", "event", "Event", "d"), names(ipd_df))[1]
  
  ipd_clean <- data.frame(
    time = as.numeric(ipd_df[[time_col]]),
    status = as.integer(ipd_df[[status_col]] == 1)
  )
  
  # 贴上该队列的专属标签，用于后续合并与亚组分析
  if(nrow(meta) > 0){
    ipd_clean$sheet_name <- sheet_name
    ipd_clean$study_id <- meta$study_id[1]
    ipd_clean$analysis_set <- meta$analysis_set[1]
    ipd_clean$size_group <- meta$size_group[1]
  } else {
    warning(paste("未在 basic_info 找到", sheet_name, "的元数据，标签赋为空。"))
    ipd_clean$sheet_name <- sheet_name
    ipd_clean$study_id <- "Unknown"
    ipd_clean$analysis_set <- "Unknown"
    ipd_clean$size_group <- "Unknown"
  }
  
  ipd_results_list[[sheet_name]] <- ipd_clean
  
  # ================= 原生校准绘图 (plotIPD) =================
  native_pdf <- file.path(out_dir, "01_IPDfromKM_native_plots", paste0(sheet_name, "_plotIPD.pdf"))
  tryCatch({
    pdf(native_pdf, width = 8, height = 6)
    plotIPD(est)
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
  })
}

# 垂直合并所有的子队列，形成终极主干数据集
final_master_ipd <- bind_rows(ipd_results_list)

cat("\n所有曲线重建并合并完成！主干数据集共有", nrow(final_master_ipd), "名虚拟患者。\n")
write.csv(final_master_ipd, file.path(out_dir, "Reconstructed_Master_IPD.csv"), row.names = FALSE)


# -------------------------------------------------------------------------
# 3. 生存分析与 K-M 曲线重绘
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 3. 生存分析与 K-M 曲线重绘
# -------------------------------------------------------------------------
cat("\n=== 开始进行生存分析与绘图 ===\n")

# -- 分析 1: 整体队列的无 WF/HRS 生存期分析 --
cat("正在生成主队列生存曲线...\n")
primary_ipd <- final_master_ipd %>%
  filter(analysis_set == "primary")

fit_primary <- survfit(Surv(time, status) ~ 1, data = primary_ipd)

p_primary <- ggsurvplot(
  fit_primary,
  data = primary_ipd, # 注意：此处建议直接传 primary_ipd 保持一致
  conf.int = TRUE,
  risk.table = TRUE,
  censor = TRUE,
  xlab = "Time (months)",
  ylab = "WF/HRS-free survival",
  title = "Overall Reconstructed Survival",
  ggtheme = theme_bw()
)

# 保存为 JPG 格式
ggsave(file.path(out_dir, "03_survival_curves", "overall_survival.jpg"), 
       print(p_primary), width = 8, height = 6, dpi = 300)

# 保存为 PDF 格式
ggsave(file.path(out_dir, "03_survival_curves", "overall_survival.pdf"), 
       print(p_primary), width = 8, height = 6)


# -- 分析 2: 按尺寸 (Size) 的亚组分析 --
subgroup_data <- final_master_ipd %>%
  filter(
    analysis_set == "primary",
    size_group %in% c("<1 cm", "1-2 cm", "2-3 cm")
  ) %>%
  mutate(size_group = factor(size_group, levels = c("<1 cm", "1-2 cm", "2-3 cm")))

if(nrow(subgroup_data) > 0) {
  cat("正在进行按尺寸分组的亚组分析...\n")
  
  fit_size <- survfit(Surv(time, status) ~ size_group, data = subgroup_data)
  
  p_size <- ggsurvplot(
    fit_size,
    data = subgroup_data,
    pval = TRUE,              # 显示 Log-rank 检验 p 值
    conf.int = TRUE,
    risk.table = TRUE,
    legend.title = "Cyst Size",
    xlab = "Time (months)",
    ylab = "WF/HRS-free survival",
    title = "Survival by Cyst Size",
    ggtheme = theme_bw()
  )
  
  # 保存为 JPG 格式
  ggsave(file.path(out_dir, "03_survival_curves", "subgroup_size_survival.jpg"), 
         print(p_size), width = 9, height = 6.5, dpi = 300)
  
  # 保存为 PDF 格式
  ggsave(file.path(out_dir, "03_survival_curves", "subgroup_size_survival.pdf"), 
         print(p_size), width = 9, height = 6.5)
  
  # -------------------------------------------------------------------------
  # 4. Cox 回归模型与比例风险 (PH) 假设检验
  # -------------------------------------------------------------------------
  cat("\n=== 开始拟合 Cox 模型并检验 PH 假设 ===\n")
  
  # 拟合 Cox 模型，以 <1 cm 为参照组 (因为上面 factor 的 level 设了 <1 cm 为第一位)
  cox_model <- coxph(
    Surv(time, status) ~ size_group + strata(study_id),
    data = subgroup_data
  )
  
  # 提取并导出 HR (Hazard Ratio) 表格
  cox_summary <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, HR = estimate, conf.low, conf.high, p.value)
  
  cat("Cox 模型 HR 结果:\n")
  print(cox_summary)
  write.csv(cox_summary, file.path(out_dir, "Cox_HR_results.csv"), row.names = FALSE)
  
  # 执行 Schoenfeld 残差检验 (PH 假设检验)
  ph_test <- cox.zph(cox_model)
  cat("\nPH 假设检验结果 (P值 < 0.05 意味着违反了等比例风险假设):\n")
  print(ph_test)
  
  # 导出 PH 检验表格
  ph_table <- as.data.frame(ph_test$table) %>% tibble::rownames_to_column("Variable")
  write.csv(ph_table, file.path(out_dir, "PH_Assumption_Test.csv"), row.names = FALSE)
  
  # 绘制 Schoenfeld 残差图 (直观检查随时间变化的 HR 趋势)
  pdf(file.path(out_dir, "03_survival_curves", "PH_Schoenfeld_Residuals.pdf"), width = 8, height = 5)
  plot(ph_test, main="Schoenfeld Residuals (PH Test)")
  dev.off()
  global_ph_p <- ph_test$table["GLOBAL", "p"]
  
  if (!is.na(global_ph_p) && global_ph_p < 0.05) {
    subgroup_data <- subgroup_data %>%
      mutate(
        size_1_2 = ifelse(size_group == "1-2 cm", 1, 0),
        size_2_3 = ifelse(size_group == "2-3 cm", 1, 0)
      )
    
    tvc_model <- coxph(
      Surv(time, status) ~ size_1_2 + size_2_3 +
        tt(size_1_2) + tt(size_2_3) + strata(study_id),
      data = subgroup_data,
      tt = function(x, t, ...) x * log(pmax(t, 0.01))
    )
    
    print(summary(tvc_model))
  }
} else {
  warning("未能找到用于亚组分析的尺寸数据标签 (<1 cm, 1-2 cm, 2-3 cm)，跳过亚组分析。")
}

cat("\n==================================\n")
cat("全部分析完成！结果已保存在", out_dir, "文件夹中。\n")


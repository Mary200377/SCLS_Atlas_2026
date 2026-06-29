# ============================================================================
# АННОТАЦИЯ КЛЕТОЧНЫХ ТИПОВ SCLC
# На основе объекта sclc_epithelial.rds
# ============================================================================

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(harmony)
library(future)

# Увеличиваем лимит памяти для future
options(future.globals.maxSize = 8000 * 1024^2)
plan("multisession", workers = 4)

set.seed(1234)
setwd("~/Magister_work/")

# ============================================================================
# 1. ЗАГРУЗКА ДАННЫХ
# ============================================================================

epithelial_cells <- readRDS("./results/sclc_epithelial.rds")

cat("Клеток:", ncol(epithelial_cells), "\n")
cat("Генов:", nrow(epithelial_cells), "\n")
cat("Кластеров SCLC:", length(unique(epithelial_cells$SCLC_clusters)), "\n")

# ============================================================================
# 2. РАСЧЕТ СКОРОВ ДЛЯ СУБТИПОВ SCLC
# ============================================================================

# Переключаемся на SCT assay (где хранятся нормализованные данные)
DefaultAssay(epithelial_cells) <- "SCT"

# Маркеры субтипов SCLC (из литературы: George et al. 2015, Gay et al. 2021)
sclc_subtype_markers <- list(
  SCLC_A = c("ASCL1", "DLL3", "INSM1", "SCG2", "SCG3", "CHGA", "SYP", 
             "NCAM1", "BEX1", "HES6", "SOX2", "POU3F2"),
  SCLC_N = c("NEUROD1", "NEUROG1", "NEUROG2", "ASCL2", "NKX2-1", "TTF1",
             "BCL2", "MYCL", "REST", "SOX4", "PAX8"),
  SCLC_P = c("POU2F3", "POU2F1", "POU2F2", "TRPM5", "GNAT3", "PLCB2",
             "CLDN4", "IL25", "DCLK1", "GUCA2B"),
  SCLC_Y = c("YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
             "VGLL4", "FOSL1", "FOSL2", "JUNB", "JUND")
)

# Фильтруем маркеры, которые есть в данных
found_markers <- list()
for (subtype in names(sclc_subtype_markers)) {
  found <- intersect(sclc_subtype_markers[[subtype]], rownames(epithelial_cells))
  if (length(found) > 0) {
    found_markers[[subtype]] <- found
    cat(sprintf("  %s: %d/%d генов найдено\n", subtype, length(found), 
                length(sclc_subtype_markers[[subtype]])))
  } else {
    cat(sprintf(" %s: гены не найдены\n", subtype))
  }
}

# Рассчитываем скоры для каждого субтипа
for (subtype in names(found_markers)) {
  epithelial_cells <- AddModuleScore(
    epithelial_cells,
    features = list(found_markers[[subtype]]),
    name = paste0(subtype, "_score"),
    ctrl = 50  # количество контрольных генов
  )
  cat(sprintf(" Скор %s рассчитан\n", subtype))
}

# ============================================================================
# 3. АННОТАЦИЯ КЛАСТЕРОВ
# ============================================================================

# Для каждого кластера определяем доминирующий субтип
cluster_annotations <- data.frame(
  cluster = unique(epithelial_cells$SCLC_clusters),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(cluster_annotations)) {
  cluster_id <- cluster_annotations$cluster[i]
  cluster_mask <- epithelial_cells$SCLC_clusters == cluster_id
  
  # Вычисляем средний скор для каждого субтипа в кластере
  scores <- sapply(names(found_markers), function(subtype) {
    score_col <- paste0(subtype, "_score1")
    mean(epithelial_cells[[score_col]][cluster_mask])
  })
  
  # Определяем доминирующий субтип
  dominant_subtype <- names(which.max(scores))
  max_score <- max(scores)
  
  cluster_annotations$dominant_subtype[i] <- dominant_subtype
  cluster_annotations$max_score[i] <- max_score
  
  cat(sprintf("  Кластер %s: %s (скор: %.3f)\n", 
              cluster_id, dominant_subtype, max_score))
}

# Добавляем аннотации в метаданные
epithelial_cells$SCLC_subtype <- cluster_annotations$dominant_subtype[
  match(epithelial_cells$SCLC_clusters, cluster_annotations$cluster)
]

# ============================================================================
# 4. ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ
# ============================================================================


# UMAP по субтипам
p_subtype <- DimPlot(epithelial_cells, group.by = "SCLC_subtype",
                     label = TRUE, label.size = 4, pt.size = 0.5) +
  ggtitle("SCLC Subtypes (A/N/P/Y)")

# UMAP по скорам
p_scores <- FeaturePlot(epithelial_cells, 
                        features = c("SCLC_A_score1", "SCLC_N_score1", 
                                    "SCLC_P_score1", "SCLC_Y_score1"),
                        ncol = 2, pt.size = 0.5)

# DotPlot маркеров
dotplot <- DotPlot(epithelial_cells, 
                   features = found_markers,
                   group.by = "SCLC_subtype") +
  RotatedAxis() +
  ggtitle("SCLC Subtype Markers")

# Сохранение визуализаций
png("./figures/sclc_subtype_umap.png", width = 12, height = 10, units = "in", res = 300)
print(p_subtype)
dev.off()

png("./figures/sclc_subtype_scores.png", width = 14, height = 12, units = "in", res = 300)
print(p_scores)
dev.off()

png("./figures/sclc_subtype_dotplot.png", width = 14, height = 10, units = "in", res = 300)
print(dotplot)
dev.off()


# ============================================================================
# 5. СТАТИСТИКА ПО СУБТИПАМ
# ============================================================================


subtype_stats <- epithelial_cells@meta.data %>%
  group_by(SCLC_subtype) %>%
  summarise(
    n_cells = n(),
    percent = n() / nrow(epithelial_cells@meta.data) * 100
  ) %>%
  arrange(desc(n_cells))

print(subtype_stats)

# Статистика по кластерам
cluster_stats <- epithelial_cells@meta.data %>%
  group_by(SCLC_clusters, SCLC_subtype) %>%
  summarise(n_cells = n()) %>%
  ungroup()

write.csv(cluster_stats, "./results/sclc_cluster_stats.csv", row.names = FALSE)

# ============================================================================
# 6. ДОПОЛНИТЕЛЬНЫЙ АНАЛИЗ: EPCAM ЭКСПРЕССИЯ
# ============================================================================


# Согласно Zhang et al. 2022, ASCL1+ клетки имеют высокую экспрессию EPCAM
if ("EPCAM" %in% rownames(epithelial_cells)) {
  epcam_by_subtype <- epithelial_cells@meta.data %>%
    group_by(SCLC_subtype) %>%
    summarise(
      mean_epcam = mean(EPCAM, na.rm = TRUE),
      median_epcam = median(EPCAM, na.rm = TRUE)
    )
  
  print(epcam_by_subtype)
  
  # Визуализация
  p_epcam <- VlnPlot(epithelial_cells, features = "EPCAM",
                     group.by = "SCLC_subtype", pt.size = 0) +
    ggtitle("EPCAM Expression by SCLC Subtype")
  
  png("./figures/sclc_epcam_violin.png", width = 12, height = 8, units = "in", res = 300)
  print(p_epcam)
  dev.off()
  
}

# ============================================================================
# 7. СОХРАНЕНИЕ АННОТИРОВАННОГО ОБЪЕКТА
# ============================================================================

# Сохраняем аннотированный объект
saveRDS(epithelial_cells, "./results/sclc_epithelial_annotated.rds")

# Сохраняем метаданные
write.csv(epithelial_cells@meta.data, "./results/sclc_epithelial_metadata.csv", 
          row.names = TRUE)

# Сохраняем аннотации кластеров
write.csv(cluster_annotations, "./results/sclc_cluster_annotations.csv", 
          row.names = FALSE)

# ============================================================================
# 8. ДОПОЛНИТЕЛЬНАЯ ВИЗУАЛИЗАЦИЯ: HEATMAP МАРКЕРОВ
# ============================================================================


# Выбираем топ маркеры для каждого субтипа
top_markers <- list()
for (subtype in names(found_markers)) {
  top_markers[[subtype]] <- head(found_markers[[subtype]], 5)
}

# Heatmap
heatmap <- DoHeatmap(epithelial_cells, 
                     features = unlist(top_markers),
                     group.by = "SCLC_subtype",
                     size = 3) +
  ggtitle("Top SCLC Subtype Markers")

png("./figures/sclc_subtype_heatmap.png", width = 10, height = 12, units = "in", res = 300)
print(heatmap)
dev.off()

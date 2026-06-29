# ============================================================================
# SCLC scRNA-seq Analysis Pipeline (16 samples)
# ============================================================================

# Загрузка необходимых библиотек
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(harmony)
library(future)
library(cowplot)

# Установка количества потоков для параллельных вычислений
plan(multisession, workers = 8)
set.seed(1234)

# Рабочая директория
setwd("~/storage/Magister_work/")

# ============================================================================
# 0. ЗАГРУЗКА ДАННЫХ (16 образцов SCLC)
# ============================================================================

# Список образцов
sample_ids <- c(
  "SRR13342017", "SRR13342018", "SRR13342019", "SRR13342020",
  "SRR13342021", "SRR13342022", "SRR13342023", "SRR13342024",
  "SRR13417495", "SRR13417496", "SRR13417497", "SRR13417498",
  "SRR13417499", "SRR13417500", "SRR13417501", "SRR13417502"
)

# Создание списка Seurat объектов
seurat_list <- list()

for (sample in sample_ids) {
  cat("Загрузка образца:", sample, "\n")
  
  # Путь к данным Cell Ranger
  data_path <- paste0("./cellranger_out/", sample, "/outs/filtered_feature_bc_matrix/")
  
  # Чтение данных (10x Genomics формат)
  seurat_list[[sample]] <- Read10X(data.dir = data_path)
  
  # Создание Seurat объекта с метаданными
  seurat_list[[sample]] <- CreateSeuratObject(
    counts = seurat_list[[sample]],
    project = "SCLC",
    min.cells = 3,
    min.features = 200
  )
  
  # Добавление информации об образце
  seurat_list[[sample]]$sample <- sample
  seurat_list[[sample]]$batch <- sample
}


# ============================================================================
# 1. КОНТРОЛЬ КАЧЕСТВА (QC) И ФИЛЬТРАЦИЯ
# ============================================================================

# Расчет метрик QC для каждого образца
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^MT-", col.name = "percent.mt")
  seurat_list[[i]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^RP[SL]", col.name = "percent.ribo")
}

# Визуализация QC метрик
qc_plots <- list()
for (i in 1:length(seurat_list)) {
  sample_name <- names(seurat_list)[i]
  
  p1 <- VlnPlot(seurat_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0) + ggtitle(sample_name)
  qc_plots[[i]] <- p1
}

# Сохранение QC plots
png("./figures/01_QC_violin.png", width = 20, height = 15, units = "in", res = 300)
print(wrap_plots(qc_plots, ncol = 4))
dev.off()

# Фильтрация клеток (пороги из статьи Zhang et al.)
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- subset(
    seurat_list[[i]],
    subset = nFeature_RNA > 200 & 
             nFeature_RNA < 9000 & 
             percent.mt < 25
  )
  cat("Образец", names(seurat_list)[i], ":", ncol(seurat_list[[i]]), "клеток после фильтрации\n")
}

# ============================================================================
# 2. ИНТЕГРАЦИЯ ДАННЫХ (Harmony)
# ============================================================================

# Объединение всех образцов в один объект
sclc_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "SCLC_Integrated"
)

cat("Всего клеток после объединения:", ncol(sclc_merged), "\n")

# Стандартная нормализация для интеграции
DefaultAssay(sclc_merged) <- "RNA"
sclc_merged <- NormalizeData(sclc_merged, normalization.method = "LogNormalize", scale.factor = 10000)

# Поиск вариабельных генов
sclc_merged <- FindVariableFeatures(sclc_merged, selection.method = "vst", nfeatures = 3000)

# Масштабирование
all_genes <- rownames(sclc_merged)
sclc_merged <- ScaleData(sclc_merged, features = all_genes)

# PCA
sclc_merged <- RunPCA(sclc_merged, features = VariableFeatures(sclc_merged), npcs = 50)

# Интеграция с помощью Harmony (удаление батч-эффекта)
cat("Запуск Harmony для интеграции...\n")
sclc_merged <- RunHarmony(sclc_merged, group.by.vars = "batch", assay.use = "RNA", reduction = "pca")

# Сохранение интегрированных данных
saveRDS(sclc_merged, file = "./results/sclc_merged_harmony.rds")

# ============================================================================
# 3. КЛАСТЕРИЗАЦИЯ И UMAP НА ИНТЕГРИРОВАННЫХ ДАННЫХ
# ============================================================================

# Построение графа соседей на интегрированных данных (Harmony embeddings)
sclc_merged <- FindNeighbors(sclc_merged, reduction = "harmony", dims = 1:30)

# Кластеризация Leiden
sclc_merged <- FindClusters(sclc_merged, resolution = 0.5, algorithm = 4, cluster.name = "RNA_clusters")

# UMAP
sclc_merged <- RunUMAP(sclc_merged, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")

# Визуализация
p1 <- DimPlot(sclc_merged, reduction = "umap_harmony", group.by = "sample", 
              label = FALSE, pt.size = 0.5) + ggtitle("Batch (Sample)")
p2 <- DimPlot(sclc_merged, reduction = "umap_harmony", group.by = "RNA_clusters", 
              label = TRUE, label.size = 4, pt.size = 0.5) + ggtitle("Clusters")

png("./figures/02_integration_umap.png", width = 16, height = 8, units = "in", res = 300)
print(p1 + p2)
dev.off()

cat("Найдено кластеров:", length(unique(sclc_merged$RNA_clusters)), "\n")

# ============================================================================
# 4. SCTransform PIPELINE (Альтернативный подход)
# ============================================================================

# SCTransform для каждого образца отдельно
sclc_list_sct <- list()
for (i in 1:length(seurat_list)) {
  sample_name <- names(seurat_list)[i]
  cat("SCTransform для образца:", sample_name, "\n")
  
  DefaultAssay(seurat_list[[i]]) <- "RNA"
  seurat_list[[i]] <- SCTransform(
    seurat_list[[i]], 
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
  sclc_list_sct[[i]] <- seurat_list[[i]]
}

# Интеграция SCTransform объектов
sclc_sct_anchors <- FindIntegrationAnchors(
  object.list = sclc_list_sct,
  normalization.method = "SCT",
  anchor.features = 3000
)

sclc_sct_integrated <- IntegrateData(
  anchorset = sclc_sct_anchors,
  normalization.method = "SCT"
)

DefaultAssay(sclc_sct_integrated) <- "integrated"

# PCA на интегрированных данных
sclc_sct_integrated <- RunPCA(sclc_sct_integrated, npcs = 50, verbose = FALSE)

# Harmony интеграция
sclc_sct_integrated <- RunHarmony(sclc_sct_integrated, group.by.vars = "sample", reduction = "pca")

# Кластеризация и UMAP
sclc_sct_integrated <- FindNeighbors(sclc_sct_integrated, reduction = "harmony", dims = 1:30)
sclc_sct_integrated <- FindClusters(sclc_sct_integrated, resolution = 0.5, algorithm = 4, cluster.name = "SCT_clusters")
sclc_sct_integrated <- RunUMAP(sclc_sct_integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap_sct")

# Визуализация SCT
p3 <- DimPlot(sclc_sct_integrated, reduction = "umap_sct", group.by = "sample", 
              label = FALSE, pt.size = 0.5) + ggtitle("SCT Batch")
p4 <- DimPlot(sclc_sct_integrated, reduction = "umap_sct", group.by = "SCT_clusters", 
              label = TRUE, label.size = 4, pt.size = 0.5) + ggtitle("SCT Clusters")

png("./figures/03_sct_umap.png", width = 16, height = 8, units = "in", res = 300)
print(p3 + p4)
dev.off()

# Сохранение SCT объекта
saveRDS(sclc_sct_integrated, file = "./results/sclc_sct_integrated.rds")

# ============================================================================
# 5. АННОТАЦИЯ КЛЕТОЧНЫХ ТИПОВ
# ============================================================================

# Маркеры клеточных типов
cell_markers <- list(
  T_cells = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD19"),
  NK_cells = c("NCAM1", "NKG7", "KLRD1"),
  Myeloid = c("CD68", "CD14", "LYZ", "FCGR3A"),
  Fibroblasts = c("COL1A1", "COL1A2", "DCN", "LUM", "FAP"),
  Endothelial = c("PECAM1", "VWF", "CDH5"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "CDH1"),
  SCLC = c("ASCL1", "NEUROD1", "CHGA", "SYP", "INSM1", "SOX2")
)

# Визуализация маркеров
DefaultAssay(sclc_merged) <- "RNA"

# DotPlot
dotplot <- DotPlot(sclc_merged, features = cell_markers, group.by = "RNA_clusters") + 
  RotatedAxis() + ggtitle("Cell Type Markers")

png("./figures/04_dotplot_markers.png", width = 12, height = 10, units = "in", res = 300)
print(dotplot)
dev.off()

# FeaturePlot для ключевых маркеров
feature_plots <- FeaturePlot(
  sclc_merged, 
  features = c("EPCAM", "PTPRC", "COL1A1", "ASCL1", "NEUROD1"),
  reduction = "umap_harmony",
  ncol = 3
)

png("./figures/05_featureplot_markers.png", width = 15, height = 10, units = "in", res = 300)
print(feature_plots)
dev.off()

# ============================================================================
# 6. АНАЛИЗ СУБТИПОВ SCLC
# ============================================================================

# Фильтрация только эпителиальных/SCLC клеток
epithelial_cells <- subset(sclc_merged, subset = EPCAM > 1)

# Перекластеризация
epithelial_cells <- FindVariableFeatures(epithelial_cells, nfeatures = 2000)
epithelial_cells <- ScaleData(epithelial_cells)
epithelial_cells <- RunPCA(epithelial_cells, npcs = 30)
epithelial_cells <- RunHarmony(epithelial_cells, group.by.vars = "sample")
epithelial_cells <- FindNeighbors(epithelial_cells, reduction = "harmony", dims = 1:20)
epithelial_cells <- FindClusters(epithelial_cells, resolution = 0.6, cluster.name = "SCLC_clusters")
epithelial_cells <- RunUMAP(epithelial_cells, reduction = "harmony", dims = 1:20)

# Визуализация субтипов
p5 <- DimPlot(epithelial_cells, group.by = "SCLC_clusters", label = TRUE, pt.size = 0.5)
p6 <- FeaturePlot(epithelial_cells, features = c("ASCL1", "NEUROD1", "POU2F3", "YAP1"), ncol = 2)

png("./figures/06_sclc_subtypes.png", width = 16, height = 12, units = "in", res = 300)
print(p5 + p6)
dev.off()

# Сохранение
saveRDS(epithelial_cells, file = "./results/sclc_epithelial.rds")

# ============================================================================
# 7. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
# ============================================================================

# Метаданные всех клеток
metadata <- sclc_merged@meta.data
write.csv(metadata, file = "./results/cell_metadata.csv")

# Информация о кластерах
cluster_info <- data.frame(
  cluster = unique(sclc_merged$RNA_clusters),
  n_cells = table(sclc_merged$RNA_clusters)
)
write.csv(cluster_info, file = "./results/cluster_info.csv")

# ============================================================================
# ДОПОЛНИТЕЛЬНО: Cell Cycle Scoring
# ============================================================================

# Гены клеточного цикла
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Пересечение с нашими данными
s.genes <- intersect(s.genes, rownames(sclc_merged))
g2m.genes <- intersect(g2m.genes, rownames(sclc_merged))

sclc_merged <- CellCycleScoring(
  sclc_merged,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# Визуализация клеточного цикла
p7 <- VlnPlot(sclc_merged, features = c("S.Score", "G2M.Score"), group.by = "RNA_clusters", pt.size = 0)
p8 <- DimPlot(sclc_merged, reduction = "umap_harmony", group.by = "Phase", pt.size = 0.5)

png("./figures/07_cell_cycle.png", width = 16, height = 8, units = "in", res = 300)
print(p7 + p8)
dev.off()

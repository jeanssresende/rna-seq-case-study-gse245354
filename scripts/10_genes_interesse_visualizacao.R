# ==========================================================
# BLOCO 0: Carregar dados da Lesson 09
# ==========================================================

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Carregar a tabela integrada (formato longo)
expr_metadata <- read.csv(
  "data/processed/expr_genes_interesse_long.csv"
)

# Carregar a submatriz (formato matriz)
expr_genes_interesse <- read.csv(
  "data/processed/expr_genes_interesse_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

cat("✓ Dados carregados da Lesson 09\n")
cat("  expr_metadata:", dim(expr_metadata), "\n")
cat("  expr_genes_interesse:", dim(expr_genes_interesse), "\n\n")

# Verificar estrutura
cat("Primeiras linhas de expr_metadata:\n")
print(head(expr_metadata, 10))

# ==========================================================
# BLOCO 1: Exploração inicial
# ==========================================================

cat("\n=== EXPLORAÇÃO INICIAL ===\n")

# Genes únicos
genes_unicos <- unique(expr_metadata$Gene)
cat("Genes analisados:", length(genes_unicos), "\n")
print(genes_unicos)

# Resumo estatístico por gene
cat("\nResumo de expressão (CPM) por gene:\n")
resumo_genes <- expr_metadata %>%
  group_by(Gene) %>%
  summarise(
    Mean = mean(CPM),
    Median = median(CPM),
    SD = sd(CPM),
    Min = min(CPM),
    Max = max(CPM),
    .groups = 'drop'
  ) %>%
  arrange(desc(Mean))

print(resumo_genes)

# Salvar resumo
write.csv(resumo_genes,
          file = "results/lesson10_gene_expression_summary.csv",
          row.names = FALSE)

cat("\n✓ Resumo salvo em: results/lesson10_gene_expression_summary.csv\n")

# ==========================================================
# BLOCO 2: Distribuição geral de expressão
# ==========================================================

cat("\n=== DISTRIBUIÇÃO GERAL DE EXPRESSÃO ===\n")

# Histograma de todos os valores de CPM
hist_all <- ggplot(expr_metadata, aes(x = CPM)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of CPM Values (All 20 Genes)",
    x = "CPM (Counts Per Million)",
    y = "Frequency"
  )

print(hist_all)

ggsave("results/lesson10_01_histogram_all_cpm.png",
       hist_all,
       width = 7,
       height = 5,
       dpi = 300)

cat("✓ Histograma salvo em: results/lesson10_01_histogram_all_cpm.png\n")

# Density plot de todos os valores
density_all <- ggplot(expr_metadata, aes(x = CPM)) +
  geom_density(fill = "steelblue", alpha = 0.6, color = "black") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Density Distribution of CPM Values",
    x = "CPM",
    y = "Density"
  )

print(density_all)

ggsave("results/lesson10_02_density_all_cpm.png",
       density_all,
       width = 7,
       height = 5,
       dpi = 300)

cat("✓ Density plot salvo em: results/lesson10_02_density_all_cpm.png\n")

# ==========================================================
# BLOCO 3: Boxplot de todos os genes
# ==========================================================

cat("\n=== BOXPLOT DE EXPRESSÃO POR GENE ===\n")

# Reordenar genes por mediana de expressão
gene_order <- expr_metadata %>%
  group_by(Gene) %>%
  summarise(Median = median(CPM), .groups = 'drop') %>%
  arrange(desc(Median)) %>%
  pull(Gene)

expr_metadata$Gene <- factor(expr_metadata$Gene, levels = gene_order)

# Boxplot
boxplot_genes <- ggplot(expr_metadata, aes(x = Gene, y = CPM, fill = Gene)) +
  geom_boxplot(alpha = 0.7, outlier.size = 2) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Expression Distribution of 20 Target Genes",
    x = "Gene",
    y = "CPM (Counts Per Million)"
  )

print(boxplot_genes)

ggsave("results/lesson10_03_boxplot_genes.png",
       boxplot_genes,
       width = 10,
       height = 6,
       dpi = 300)

cat("✓ Boxplot salvo em: results/lesson10_03_boxplot_genes.png\n")

# ==========================================================
# BLOCO 4: Boxplot com jitter (amostra individual)
# ==========================================================

cat("\n=== BOXPLOT + JITTER ===\n")

boxplot_jitter <- ggplot(expr_metadata, aes(x = Gene, y = CPM, fill = Gene)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Expression Distribution with Individual Samples",
    x = "Gene",
    y = "CPM"
  )

print(boxplot_jitter)

ggsave("results/lesson10_04_boxplot_jitter.png",
       boxplot_jitter,
       width = 10,
       height = 6,
       dpi = 300)

cat("✓ Boxplot + jitter salvo em: results/lesson10_04_boxplot_jitter.png\n")

# ==========================================================
# BLOCO 5: Violin plot
# ==========================================================

cat("\n=== VIOLIN PLOT ===\n")

violin_genes <- ggplot(expr_metadata, aes(x = Gene, y = CPM, fill = Gene)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Expression Distribution (Violin Plot)",
    x = "Gene",
    y = "CPM"
  )

print(violin_genes)

ggsave("results/lesson10_05_violin_genes.png",
       violin_genes,
       width = 10,
       height = 6,
       dpi = 300)

cat("✓ Violin plot salvo em: results/lesson10_05_violin_genes.png\n")

# ==========================================================
# BLOCO 6: Facet por gene (distribuição por tratamento)
# ==========================================================

cat("\n=== FACET POR GENE - DISTRIBUIÇÃO POR TRATAMENTO ===\n")

facet_treatment <- ggplot(expr_metadata, aes(x = Treatment, y = CPM, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 5) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Expression by Treatment (20 Target Genes)",
    x = "Treatment",
    y = "CPM"
  )

print(facet_treatment)

ggsave("results/lesson10_06_facet_treatment.png",
       facet_treatment,
       width = 14,
       height = 10,
       dpi = 300)

cat("✓ Facet treatment salvo em: results/lesson10_06_facet_treatment.png\n")

# ==========================================================
# BLOCO 7: Facet por gene (distribuição por tempo)
# ==========================================================

cat("\n=== FACET POR GENE - DISTRIBUIÇÃO POR TEMPO ===\n")

facet_time <- ggplot(expr_metadata, aes(x = Time, y = CPM, fill = Time)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 5) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Expression by Time (20 Target Genes)",
    x = "Time",
    y = "CPM"
  )

print(facet_time)

ggsave("results/lesson10_07_facet_time.png",
       facet_time,
       width = 14,
       height = 10,
       dpi = 300)

cat("✓ Facet time salvo em: results/lesson10_07_facet_time.png\n")

# ==========================================================
# BLOCO 8: Facet por gene (distribuição por doador)
# ==========================================================

cat("\n=== FACET POR GENE - DISTRIBUIÇÃO POR DOADOR ===\n")

# Converter Donor para factor para garantir ordem correta
expr_metadata$Donor <- factor(expr_metadata$Donor, 
                              levels = c("1", "2", "3"))

facet_donor <- ggplot(expr_metadata, aes(x = Donor, y = CPM, fill = Donor)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 5) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    title = "Expression by Donor (20 Target Genes)",
    x = "Donor",
    y = "CPM"
  )

print(facet_donor)

ggsave("results/lesson10_08_facet_donor.png",
       facet_donor,
       width = 14,
       height = 10,
       dpi = 300)

cat("✓ Facet donor salvo em: results/lesson10_08_facet_donor.png\n")

# ==========================================================
# BLOCO 9: Heatmap dos 20 genes (sem transformação log)
# ==========================================================

cat("\n=== HEATMAP DOS 20 GENES ===\n")

library(pheatmap)

# Preparar matriz para heatmap
heatmap_matrix <- expr_genes_interesse

# Criar anotação de colunas
metadata_samples <- expr_metadata %>%
  select(MatrixSample, Treatment, Time, Donor) %>%
  distinct() %>%
  as.data.frame()

rownames(metadata_samples) <- metadata_samples$MatrixSample
metadata_samples$MatrixSample <- NULL

# Converter para factors com ordem correta
metadata_samples$Treatment <- factor(metadata_samples$Treatment,
                                     levels = c("Ctrl", "Hb2+", "Hb3+", "Hemin"))
metadata_samples$Time <- factor(metadata_samples$Time,
                                levels = c("30 min", "1.5 hr"))
metadata_samples$Donor <- factor(metadata_samples$Donor,
                                 levels = c("1", "2", "3"))

# Reordenar para corresponder às colunas da matriz
metadata_samples <- metadata_samples[colnames(heatmap_matrix), ]

# Reordenar metadata_samples para que Donor seja a última coluna
metadata_samples <- metadata_samples %>%
  select(Treatment, Time, Donor)

# Criar cores para anotação
# IMPORTANTE: os nomes das cores devem corresponder aos níveis dos factors
annotation_colors <- list(
  Treatment = c(
    "Ctrl" = "#1b9e77",
    "Hb2+" = "#d95f02",
    "Hb3+" = "#7570b3",
    "Hemin" = "#e7298a"
  ),
  Time = c(
    "30 min" = "#66c2a5",
    "1.5 hr" = "#fc8d62"
  ),
  Donor = c(
    "1" = "#8da0cb",
    "2" = "#e78ac3",
    "3" = "#a6d854"
  )
)

# Verificar correspondência
cat("Níveis de Donor em metadata_samples:\n")
print(levels(metadata_samples$Donor))

cat("\nNomes de cores para Donor:\n")
print(names(annotation_colors$Donor))

# Criar heatmap
png("results/lesson10_09_heatmap_genes.png",
    width = 10,
    height = 8,
    units = "in",
    res = 300)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = metadata_samples,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(
    rev(RColorBrewer::brewer.pal(9, "RdBu"))
  )(100),
  fontsize = 10,
  fontsize_col = 8,
  main = "Heatmap of 20 Target Genes (Row-Scaled CPM)",
  show_rownames = TRUE,
  show_colnames = FALSE
)

dev.off()

cat("✓ Heatmap salvo em: results/lesson10_09_heatmap_genes.png\n")

# ==========================================================
# BLOCO 10: Resumo visual (número de amostras por grupo)
# ==========================================================

cat("\n=== RESUMO DE AMOSTRAS POR GRUPO ===\n")

# Contar amostras por tratamento
summary_treatment <- expr_metadata %>%
  select(MatrixSample, Treatment) %>%
  distinct() %>%
  group_by(Treatment) %>%
  summarise(N = n(), .groups = 'drop')

cat("Amostras por tratamento:\n")
print(summary_treatment)

# Contar amostras por tempo
summary_time <- expr_metadata %>%
  select(MatrixSample, Time) %>%
  distinct() %>%
  group_by(Time) %>%
  summarise(N = n(), .groups = 'drop')

cat("\nAmostras por tempo:\n")
print(summary_time)

# Contar amostras por doador
summary_donor <- expr_metadata %>%
  select(MatrixSample, Donor) %>%
  distinct() %>%
  group_by(Donor) %>%
  summarise(N = n(), .groups = 'drop')

cat("\nAmostras por doador:\n")
print(summary_donor)

# ==========================================================
# BLOCO 11: Identificar padrões interessantes
# ==========================================================

cat("\n=== IDENTIFICAÇÃO DE PADRÕES INTERESSANTES ===\n")

# Genes com maior variabilidade entre tratamentos
var_por_gene_tratamento <- expr_metadata %>%
  group_by(Gene, Treatment) %>%
  summarise(Mean_CPM = mean(CPM), .groups = 'drop') %>%
  pivot_wider(names_from = Treatment, values_from = Mean_CPM) %>%
  mutate(
    Range = pmax(Ctrl, `Hb2+`, `Hb3+`, Hemin) - pmin(Ctrl, `Hb2+`, `Hb3+`, Hemin)
  ) %>%
  arrange(desc(Range))

cat("Genes com maior variação entre tratamentos:\n")
print(var_por_gene_tratamento)

# Genes com maior variabilidade entre tempos
var_por_gene_tempo <- expr_metadata %>%
  group_by(Gene, Time) %>%
  summarise(Mean_CPM = mean(CPM), .groups = 'drop') %>%
  pivot_wider(names_from = Time, values_from = Mean_CPM) %>%
  mutate(
    Diff = abs(`30 min` - `1.5 hr`)
  ) %>%
  arrange(desc(Diff))

cat("\nGenes com maior variação entre tempos:\n")
print(var_por_gene_tempo)

# Genes com maior variabilidade entre doadores
var_por_gene_donor <- expr_metadata %>%
  group_by(Gene, Donor) %>%
  summarise(Mean_CPM = mean(CPM), .groups = 'drop') %>%
  pivot_wider(names_from = Donor, values_from = Mean_CPM) %>%
  mutate(
    Range = pmax(`1`, `2`, `3`) - pmin(`1`, `2`, `3`)
  ) %>%
  arrange(desc(Range))

cat("\nGenes com maior variação entre doadores:\n")
print(var_por_gene_donor)

# Salvar tabelas de variação
write.csv(var_por_gene_tratamento,
          file = "results/lesson10_variation_by_treatment.csv",
          row.names = FALSE)

write.csv(var_por_gene_tempo,
          file = "results/lesson10_variation_by_time.csv",
          row.names = FALSE)

write.csv(var_por_gene_donor,
          file = "results/lesson10_variation_by_donor.csv",
          row.names = FALSE)

cat("\n✓ Tabelas de variação salvas em results/\n")

# ==========================================================
# BLOCO 12: Resumo final
# ==========================================================

cat("\n=== RESUMO DA LESSON 10 ===\n")

cat("\nVisualizações criadas:\n")
cat("  1. Histograma de distribuição geral (CPM)\n")
cat("  2. Density plot de distribuição geral\n")
cat("  3. Boxplot dos 20 genes\n")
cat("  4. Boxplot + jitter (amostras individuais)\n")
cat("  5. Violin plot dos 20 genes\n")
cat("  6. Facet por gene × tratamento\n")
cat("  7. Facet por gene × tempo\n")
cat("  8. Facet por gene × doador\n")
cat("  9. Heatmap dos 20 genes (row-scaled)\n")

cat("\nTabelas de análise criadas:\n")
cat("  - Gene expression summary (média, mediana, SD, range)\n")
cat("  - Variation by treatment\n")
cat("  - Variation by time\n")
cat("  - Variation by donor\n")

cat("\nObservações importantes:\n")
cat("  - Todos os 20 genes estão presentes\n")
cat("  - Todos os 24 tratamentos/tempo/doador estão representados\n")
cat("  - Heatmap mostra padrões de clustering por doador\n")
cat("  - Próxima etapa: análise comparativa por tratamento\n")

# ==========================================================
# FIM DO SCRIPT
# ==========================================================

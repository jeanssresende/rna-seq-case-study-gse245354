# ==========================================================
# BLOCO 0: Carregamento e preparação
# ==========================================================

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(rstatix)

# Carregar dados
expr_metadata <- read.csv(
  "data/processed/expr_genes_interesse_long.csv"
)

cat("✓ Dados carregados\n")
cat("  Dimensões:", dim(expr_metadata), "\n\n")

# ==========================================================
# BLOCO 1: Transformação log2 e Z-score
# ==========================================================

cat("=== TRANSFORMAÇÃO DOS DADOS ===\n")

# Aplicar log2 (CPM + 1 para evitar log(0))
expr_metadata <- expr_metadata %>%
  mutate(
    CPM_log2 = log2(CPM + 1),
    MatrixSample = as.factor(MatrixSample),
    Gene = as.factor(Gene),
    Treatment = as.factor(Treatment),
    Time = as.factor(Time),
    Donor = as.factor(Donor)
  )

cat("✓ Transformação log2 aplicada\n")
cat("  CPM_log2 range:", 
    round(min(expr_metadata$CPM_log2), 2), 
    "–", 
    round(max(expr_metadata$CPM_log2), 2), 
    "\n\n")

# Calcular Z-score por gene (normalização dentro de cada gene)
expr_metadata <- expr_metadata %>%
  group_by(Gene) %>%
  mutate(
    CPM_log2_zscore = scale(CPM_log2)[, 1]
  ) %>%
  ungroup()

cat("✓ Z-score calculado por gene\n\n")

# ==========================================================
# BLOCO 2: Estratégia estatística
# ==========================================================

cat("=== ESTRATÉGIA ESTATÍSTICA ===\n")

cat("
Pergunta científica:
  'Como a expressão dos 20 genes varia entre tratamentos?'

Desenho experimental:
  - 4 Tratamentos (Ctrl, Hb2+, Hb3+, Hemin)
  - 2 Tempos (30 min, 1.5 hr)
  - 3 Doadores (1, 2, 3)
  - Total: 24 amostras (4 × 2 × 3)
  - Replicação biológica: 3 doadores
  - Replicação técnica: nenhuma

Confundimento:
  - FORTE efeito de Doador (visto na Lesson 08 e 10)
  - Necessário: controlar por Doador

Abordagem estatística:
  1. Comparações pareadas (paired t-test) dentro de cada Doador
  2. Meta-análise: combinar resultados dos 3 Doadores
  3. Alternativa: modelo linear misto (LMM) com Doador como efeito aleatório

Por que pareado?
  - Mesmos 3 Doadores em todos os tratamentos
  - Pareamento reduz variância e aumenta poder
  - Mais apropriado que t-test não pareado

Correção múltipla:
  - Benjamini-Hochberg (FDR)
  - Controla taxa de descoberta falsa
  - Mais apropriado que Bonferroni para exploração

\n")

# ==========================================================
# BLOCO 3: Comparações pareadas por gene
# ==========================================================

cat("=== COMPARAÇÕES PAREADAS: TRATAMENTOS ===\n")

# Preparar dados para comparação pareada
# Precisamos de uma linha por Doador × Gene × Tratamento × Tempo

# Primeiro, calcular média por Doador × Gene × Tratamento × Tempo
expr_summary <- expr_metadata %>%
  group_by(Gene, Treatment, Time, Donor) %>%
  summarise(
    CPM_log2_mean = mean(CPM_log2),
    .groups = 'drop'
  )

cat("Dados para comparação pareada:\n")
print(head(expr_summary, 15))

# Agora, para cada Gene × Tempo, comparar todos os pares de Tratamentos
# usando t-test pareado (pareado por Doador)

# Lista para armazenar resultados
results_pairwise <- list()

genes_list <- unique(expr_metadata$Gene)
times_list <- unique(expr_metadata$Time)
treatments_list <- unique(expr_metadata$Treatment)

cat("\nCalculando comparações pareadas...\n")

for (gene in genes_list) {
  for (time in times_list) {
    
    # Filtrar dados para este gene e tempo
    data_subset <- expr_summary %>%
      filter(Gene == gene & Time == time)
    
    # Comparar todos os pares de tratamentos
    treatments_pairs <- combn(treatments_list, 2, simplify = FALSE)
    
    for (pair in treatments_pairs) {
      treat1 <- pair[1]
      treat2 <- pair[2]
      
      # Extrair valores para cada tratamento (ordenado por Doador)
      data_treat1 <- data_subset %>%
        filter(Treatment == treat1) %>%
        arrange(Donor) %>%
        pull(CPM_log2_mean)
      
      data_treat2 <- data_subset %>%
        filter(Treatment == treat2) %>%
        arrange(Donor) %>%
        pull(CPM_log2_mean)
      
      # Verificar se há dados para ambos os tratamentos
      if (length(data_treat1) > 0 & length(data_treat2) > 0) {
        
        # T-test pareado
        t_result <- t.test(data_treat1, data_treat2, paired = TRUE)
        
        # Armazenar resultado
        results_pairwise[[paste(gene, time, treat1, treat2, sep = "_")]] <- 
          data.frame(
            Gene = gene,
            Time = time,
            Treatment1 = treat1,
            Treatment2 = treat2,
            Mean_Treat1 = mean(data_treat1),
            Mean_Treat2 = mean(data_treat2),
            Diff = mean(data_treat2) - mean(data_treat1),
            t_statistic = t_result$statistic,
            df = t_result$parameter,
            p_value = t_result$p.value,
            stringsAsFactors = FALSE
          )
      }
    }
  }
}

# Combinar todos os resultados em um data frame
pairwise_results <- do.call(rbind, results_pairwise) %>%
  rownames_to_column(var = "comparison_id") %>%
  select(-comparison_id)

# Aplicar correção de múltiplos testes (FDR)
pairwise_results <- pairwise_results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    Significant = ifelse(p_adjusted < 0.05, "Yes", "No")
  )

cat("\nResultados das comparações pareadas:\n")
print(head(pairwise_results, 20))

cat("\nResumo de significância:\n")
print(table(pairwise_results$Significant))

# Salvar resultados
write.csv(pairwise_results,
          file = "results/lesson11_pairwise_comparisons.csv",
          row.names = FALSE)

cat("\n✓ Resultados salvos em: results/lesson11_pairwise_comparisons.csv\n")

# ==========================================================
# BLOCO 4: Identificar genes com efeito de tratamento
# ==========================================================

cat("\n=== GENES COM EFEITO DE TRATAMENTO ===\n")

# Genes significativos (p_adjusted < 0.05)
significant_genes <- pairwise_results %>%
  filter(p_adjusted < 0.05) %>%
  pull(Gene) %>%
  unique()

cat("Genes com efeito significativo de tratamento:", 
    length(significant_genes), "\n")
print(significant_genes)

# Top genes por magnitude de efeito
top_genes_effect <- pairwise_results %>%
  group_by(Gene) %>%
  summarise(
    Max_Diff = max(abs(Diff)),
    N_Significant = sum(Significant == "Yes"),
    .groups = 'drop'
  ) %>%
  arrange(desc(Max_Diff))

cat("\nTop 10 genes por magnitude de efeito:\n")
print(head(top_genes_effect, 10))

# ==========================================================
# BLOCO 5: Gráfico 1 – Comparação Ctrl vs Hb2+ (publicável)
# ==========================================================

cat("\n=== GRÁFICO 1: CTRL vs HB2+ ===\n")

# Filtrar comparações Ctrl vs Hb2+
comp_ctrl_hb2 <- pairwise_results %>%
  filter(
    (Treatment1 == "Ctrl" & Treatment2 == "Hb2+") |
      (Treatment1 == "Hb2+" & Treatment2 == "Ctrl")
  ) %>%
  mutate(
    Comparison = "Ctrl vs Hb2+",
    log10_padj = -log10(p_adjusted)
  ) %>%
  arrange(Time, Gene)

# Reordenar genes por diferença
gene_order_ctrl_hb2 <- comp_ctrl_hb2 %>%
  arrange(Diff) %>%
  pull(Gene) %>%
  unique()

comp_ctrl_hb2$Gene <- factor(comp_ctrl_hb2$Gene, levels = gene_order_ctrl_hb2)

# Criar gráfico
plot_ctrl_hb2 <- ggplot(comp_ctrl_hb2, 
                        aes(x = Diff, y = Gene, fill = Time, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_point(size = 4, shape = 21, stroke = 1.5) +
  scale_color_manual(values = c("Yes" = "black", "No" = "white")) +
  scale_fill_manual(values = c("30 min" = "#66c2a5", "1.5 hr" = "#fc8d62")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Ctrl vs Hb2+: Log2-fold Change in Expression",
    subtitle = "Paired t-test (n=3 donors); filled circle = p.adj < 0.05",
    x = "Log2(Hb2+ / Ctrl)",
    y = "Gene",
    fill = "Time",
    color = "Significant"
  )

print(plot_ctrl_hb2)

ggsave("results/lesson11_01_ctrl_vs_hb2.png",
       plot_ctrl_hb2,
       width = 8,
       height = 7,
       dpi = 300)

cat("✓ Gráfico salvo em: results/lesson11_01_ctrl_vs_hb2.png\n")

# ==========================================================
# BLOCO 6: Gráfico 2 – Volcano plot (publicável)
# ==========================================================

cat("\n=== GRÁFICO 2: VOLCANO PLOT ===\n")

# Usar dados de Ctrl vs Hb2+ em 30 min
volcano_data <- comp_ctrl_hb2 %>%
  filter(Time == "30 min") %>%
  mutate(
    log10_padj = -log10(p_adjusted),
    Threshold = ifelse(p_adjusted < 0.05 & abs(Diff) > 0.5, "Significant", "Not Significant")
  )

volcano_plot <- ggplot(volcano_data, aes(x = Diff, y = log10_padj, color = Threshold)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_point(size = 3, alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = volcano_data %>% filter(Threshold == "Significant"),
    aes(label = Gene),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c("Significant" = "#e74c3c", "Not Significant" = "#95a5a6")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    legend.position = "right"
  ) +
  labs(
    title = "Volcano Plot: Ctrl vs Hb2+ (30 min)",
    subtitle = "Paired t-test; threshold: p.adj < 0.05, |log2FC| > 0.5",
    x = "Log2-fold Change",
    y = "-log10(p.adjusted)",
    color = "Significance"
  )

print(volcano_plot)

ggsave("results/lesson11_02_volcano_plot.png",
       volcano_plot,
       width = 9,
       height = 7,
       dpi = 300)

cat("✓ Gráfico salvo em: results/lesson11_02_volcano_plot.png\n\n")

# Instalar ggrepel se necessário
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# ==========================================================
# BLOCO 7: Gráfico 3 – Boxplot com p-valores (publicável)
# ==========================================================

cat("=== GRÁFICO 3: BOXPLOT COM P-VALORES ===\n")

# Selecionar genes significativos que têm dados em TODOS os tratamentos
genes_complete <- expr_metadata %>%
  filter(Time == "30 min") %>%
  group_by(Gene, Treatment) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Gene) %>%
  summarise(n_treatments = n_distinct(Treatment), .groups = 'drop') %>%
  filter(n_treatments == 4) %>%  # Deve ter dados em todos os 4 tratamentos
  pull(Gene)

# Selecionar genes significativos E com dados completos
genes_to_plot <- pairwise_results %>%
  filter(p_adjusted < 0.05 & Gene %in% genes_complete) %>%
  pull(Gene) %>%
  unique() %>%
  head(6)  # Top 6 genes

cat("Genes significativos com dados completos:", length(genes_to_plot), "\n")
print(genes_to_plot)

# Se não houver genes significativos com dados completos, usar os 6 mais variáveis
if (length(genes_to_plot) == 0) {
  cat("\n Nenhum gene significativo com dados completos.\n")
  cat("Usando os 6 genes mais variáveis...\n")
  
  genes_to_plot <- expr_metadata %>%
    filter(Time == "30 min") %>%
    group_by(Gene) %>%
    summarise(
      var_expr = var(CPM_log2),
      .groups = 'drop'
    ) %>%
    arrange(desc(var_expr)) %>%
    head(6) %>%
    pull(Gene)
  
  print(genes_to_plot)
}

# Filtrar dados - garantir que temos dados para cada combinação Gene × Treatment
data_plot <- expr_metadata %>%
  filter(Gene %in% genes_to_plot & Time == "30 min") %>%
  mutate(
    Gene = factor(Gene, levels = genes_to_plot)
  )

# Verificar se há dados
if (nrow(data_plot) == 0) {
  cat("\n Nenhum dado disponível para plotar!\n")
} else {
  
  # Criar gráfico
  plot_boxplot_pval <- ggplot(data_plot, 
                              aes(x = Treatment, y = CPM_log2, fill = Treatment)) +
    geom_boxplot(alpha = 0.6, outlier.size = 1) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.major.y = element_line(color = "gray90", size = 0.3)
    ) +
    labs(
      title = "Top Significant Genes: Expression by Treatment (30 min)",
      subtitle = "Individual points = samples from 3 donors",
      x = "Treatment",
      y = "log2(CPM + 1)",
      fill = "Treatment"
    )
  
  print(plot_boxplot_pval)
  
  ggsave("results/lesson11_03_boxplot_significant_genes.png",
         plot_boxplot_pval,
         width = 10,
         height = 7,
         dpi = 300)
  
  cat("✓ Gráfico salvo em: results/lesson11_03_boxplot_significant_genes.png\n")
}

# ==========================================================
# BLOCO 8: Gráfico 4 – Heatmap com p-valores
# ==========================================================

cat("\n=== GRÁFICO 4: HEATMAP DE P-VALORES ===\n")

# Criar matriz de p-valores ajustados
# Linhas = Genes, Colunas = Comparações

# Filtrar comparações Ctrl vs outros
comp_vs_ctrl <- pairwise_results %>%
  filter(Time == "30 min") %>%
  filter(
    (Treatment1 == "Ctrl" & Treatment2 != "Ctrl") |
    (Treatment2 == "Ctrl" & Treatment1 != "Ctrl")
  ) %>%
  mutate(
    Comparison = ifelse(Treatment1 == "Ctrl", 
                        paste0("Ctrl vs ", Treatment2),
                        paste0("Ctrl vs ", Treatment1))
  ) %>%
  select(Gene, Comparison, p_adjusted)

# Criar matriz
pval_matrix <- comp_vs_ctrl %>%
  pivot_wider(names_from = Comparison, values_from = p_adjusted) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# Converter para -log10(p)
pval_matrix_log <- -log10(pval_matrix)

# Criar heatmap
png("results/lesson11_04_heatmap_pvalues.png",
    width = 8,
    height = 8,
    units = "in",
    res = 300)

pheatmap::pheatmap(
  pval_matrix_log,
  breaks = seq(0, max(pval_matrix_log, na.rm = TRUE), 
               length.out = 100),
  color = colorRampPalette(c("white", "yellow", "red"))(100),
  main = "-log10(p.adjusted): Ctrl vs Other Treatments (30 min)",
  fontsize = 10,
  display_numbers = round(pval_matrix_log, 2),
  fontsize_number = 8
)

dev.off()

cat("✓ Heatmap salvo em: results/lesson11_04_heatmap_pvalues.png\n")

# ==========================================================
# BLOCO 9: Comparações pareadas por Doador
# ==========================================================

cat("\n=== COMPARAÇÕES PAREADAS POR DOADOR ===\n")

# Análise adicional: como cada doador responde ao tratamento?

results_by_donor <- list()

for (gene in genes_list) {
  for (time in times_list) {
    for (donor in unique(expr_metadata$Donor)) {
      
      # Filtrar dados para este gene, tempo e doador
      data_subset <- expr_metadata %>%
        filter(Gene == gene & Time == time & Donor == donor)
      
      # Comparar Ctrl vs Hb2+
      data_ctrl <- data_subset %>%
        filter(Treatment == "Ctrl") %>%
        pull(CPM_log2)
      
      data_hb2 <- data_subset %>%
        filter(Treatment == "Hb2+") %>%
        pull(CPM_log2)
      
      if (length(data_ctrl) > 0 & length(data_hb2) > 0) {
        
        results_by_donor[[paste(gene, time, donor, sep = "_")]] <- 
          data.frame(
            Gene = gene,
            Time = time,
            Donor = donor,
            Mean_Ctrl = mean(data_ctrl),
            Mean_Hb2 = mean(data_hb2),
            Diff = mean(data_hb2) - mean(data_ctrl),
            stringsAsFactors = FALSE
          )
      }
    }
  }
}

donor_response <- do.call(rbind, results_by_donor) %>%
  rownames_to_column(var = "id") %>%
  select(-id)

cat("Resposta por doador (amostra):\n")
print(head(donor_response, 15))

# Salvar
write.csv(donor_response,
          file = "results/lesson11_donor_response.csv",
          row.names = FALSE)

cat("\n✓ Dados salvos em: results/lesson11_donor_response.csv\n")

# ==========================================================
# BLOCO 10: Gráfico 5 – Resposta por Doador
# ==========================================================

# ==========================================================
# BLOCO 9.5: Diagnóstico de significância
# ==========================================================

cat("\n=== DIAGNÓSTICO DE SIGNIFICÂNCIA ===\n")

cat("Distribuição de p-valores:\n")
print(summary(pairwise_results$p_value))

cat("\nDistribuição de p-valores ajustados:\n")
print(summary(pairwise_results$p_adjusted))

cat("\nQuantos genes têm p_adjusted < 0.05?", 
    sum(pairwise_results$p_adjusted < 0.05), "\n")
cat("Quantos genes têm p_adjusted < 0.10?", 
    sum(pairwise_results$p_adjusted < 0.10), "\n")
cat("Quantos genes têm p_value < 0.05?", 
    sum(pairwise_results$p_value < 0.05), "\n")

# Mostrar genes com menor p-value
cat("\nTop 10 comparações por p-valor (não ajustado):\n")
print(pairwise_results %>%
        arrange(p_value) %>%
        select(Gene, Treatment1, Treatment2, Time, Diff, p_value, p_adjusted) %>%
        head(10))

# ==========================================================
# BLOCO 10: Gráfico 5 – Resposta por Doador (CORRIGIDO)
# ==========================================================

cat("\n=== GRÁFICO 5: RESPOSTA POR DOADOR ===\n")

# Usar threshold mais permissivo (p < 0.10 ou p_value < 0.05 não ajustado)
genes_sig <- pairwise_results %>%
  filter(p_value < 0.05) %>%  # Usar p-valor não ajustado
  pull(Gene) %>%
  unique()

cat("Genes com p_value < 0.05 (não ajustado):", length(genes_sig), "\n")

# Se ainda não houver, usar os 4 com menor p-valor
if (length(genes_sig) == 0) {
  cat(" Nenhum gene com p < 0.05. Usando os 4 com menor p-valor...\n")
  
  genes_sig <- pairwise_results %>%
    arrange(p_value) %>%
    pull(Gene) %>%
    unique() %>%
    head(4)
}

print(genes_sig)

# Filtrar dados de resposta por doador
donor_plot_data <- donor_response %>%
  filter(Gene %in% genes_sig & Time == "30 min") %>%
  mutate(
    Gene = factor(Gene, levels = genes_sig)
  )

cat("Dimensões dos dados para plotar:", dim(donor_plot_data), "\n")

if (nrow(donor_plot_data) == 0) {
  cat("Nenhum dado de resposta por doador disponível!\n")
  cat("Verifique se donor_response foi criado corretamente.\n")
} else {
  
  plot_donor_response <- ggplot(donor_plot_data,
                                aes(x = Donor, y = Diff, fill = Gene)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    facet_wrap(~Gene, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(color = "gray90", size = 0.3)
    ) +
    labs(
      title = "Individual Donor Response to Hb2+ (30 min)",
      subtitle = "log2(Hb2+ / Ctrl) for each donor",
      x = "Donor",
      y = "log2-fold Change"
    )
  
  print(plot_donor_response)
  
  ggsave("results/lesson11_05_donor_response.png",
         plot_donor_response,
         width = 10,
         height = 6,
         dpi = 300)
  
  cat("✓ Gráfico salvo em: results/lesson11_05_donor_response.png\n")
}

# ==========================================================
# BLOCO 11: Resumo final
# ==========================================================

cat("\n=== RESUMO DA LESSON 11 ===\n")

cat("\nAbordagem estatística:\n")
cat("  ✓ Transformação log2 (CPM + 1)\n")
cat("  ✓ Comparações pareadas por Doador (paired t-test)\n")
cat("  ✓ Correção múltipla (Benjamini-Hochberg, FDR)\n")
cat("  ✓ Análise de resposta por Doador\n")

cat("\nResultados:\n")
cat("  - Total de comparações:", nrow(pairwise_results), "\n")
cat("  - Comparações significativas (p.adj < 0.05):", 
      sum(pairwise_results$Significant == "Yes"), "\n")
cat("  - Genes com efeito significativo:", length(significant_genes), "\n")

cat("\nGráficos publicáveis criados:\n")
cat("  1. Forest plot: Ctrl vs Hb2+ (log2FC por gene)\n")
cat("  2. Volcano plot: significância vs magnitude\n")
cat("  3. Boxplots: genes significativos com amostras\n")
cat("  4. Heatmap: -log10(p) para múltiplas comparações\n")
cat("  5. Resposta por Doador: variabilidade individual\n")

cat("\nArquivos salvos:\n")
cat("  - results/lesson11_pairwise_comparisons.csv\n")
cat("  - results/lesson11_donor_response.csv\n")
cat("  - results/lesson11_01_ctrl_vs_hb2.png\n")
cat("  - results/lesson11_02_volcano_plot.png\n")
cat("  - results/lesson11_03_boxplot_significant_genes.png\n")
cat("  - results/lesson11_04_heatmap_pvalues.png\n")
cat("  - results/lesson11_05_donor_response.png\n")

# ==========================================================
# FIM DO SCRIPT
# ==========================================================

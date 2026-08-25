# ==========================================================
# BLOCO 0: Preparação e carregamento dos dados brutos
# ==========================================================

# Limpar ambiente
rm(list = ls())

# Carregar bibliotecas
library(tidyverse)
library(dplyr)

# Carregar a matriz de expressão bruta (CPM)
counts <- read.delim(
  "data/raw/GSE245354/GSE245354_Meegan_Normalized_counts.txt.gz"
)

# Extrair apenas as colunas de expressão (começam na coluna 9)
expr <- counts[, 9:ncol(counts)]

# Atribuir nomes de genes (símbolos gênicos)
rownames(expr) <- counts$Gene.Symbol

# Remover genes sem símbolo
expr <- expr[rownames(expr) != "" & !is.na(rownames(expr)), ]

cat("✓ Matriz de expressão carregada\n")
cat("  Dimensões:", dim(expr)[1], "genes ×", dim(expr)[2], "amostras\n\n")

# Carregar metadados
metadata <- read.csv(
  "data/metadata/sample_metadata.csv",
  stringsAsFactors = FALSE
)

# Processar metadados (mesmo procedimento da Lesson 08)
metadata <- metadata %>%
  mutate(
    SampleID = str_extract(Description, "(?<=NP-)\\d+"),
    MatrixSample = paste0("X5958.NP.", SampleID, "_S1")
  ) %>%
  select(MatrixSample, Treatment, Time, Donor)

cat("✓ Metadados carregados\n")
cat("  Dimensões:", dim(metadata)[1], "amostras\n\n")

# Verificar correspondência entre nomes de colunas e metadados
cat("Primeiras amostras na matriz:\n")
print(head(colnames(expr)))

cat("\nPrimeiras amostras nos metadados:\n")
print(head(metadata$MatrixSample))

# ==========================================================
# BLOCO 1: Definir os 20 genes de interesse
# ==========================================================

genes_interesse <- c(
  "CDH5",
  "TJP1",
  "CTTN",
  "PTK2",
  "HMOX1",
  "NFE2L2",
  "GPX4",
  "SLC7A11",
  "ACSL4",
  "TFRC",
  "FTH1",
  "NCOA4",
  "AIFM2",
  "SOD2",
  "NQO1",
  "AKT1",
  "MAPK1",
  "MAPK3",
  "BMPR2",
  "TGFB1"
)

cat("\n=== DEFINIÇÃO DOS GENES DE INTERESSE ===\n")
cat("Total de genes a validar:", length(genes_interesse), "\n\n")
cat("Genes de interesse:\n")
print(genes_interesse)

# ==========================================================
# BLOCO 2: Verificar presença dos genes na matriz
# ==========================================================

genes_matriz <- rownames(expr)

cat("\n=== VERIFICAÇÃO DE PRESENÇA ===\n")
cat("Total de genes na matriz:", length(genes_matriz), "\n\n")

# Verificar quais genes estão presentes
genes_presentes <- genes_interesse %in% genes_matriz

cat("Resumo de presença:\n")
print(table(genes_presentes))

# Identificar genes encontrados e ausentes
genes_encontrados <- genes_interesse[genes_presentes]
genes_ausentes <- genes_interesse[!genes_presentes]

cat("\nGenes encontrados (", length(genes_encontrados), "):\n")
print(genes_encontrados)

if (length(genes_ausentes) > 0) {
  cat("\n ATENÇÃO: Genes NÃO encontrados (", length(genes_ausentes), "):\n")
  print(genes_ausentes)
  cat("\nEsses genes serão EXCLUÍDOS da análise.\n")
} else {
  cat("\n✓ Bom! Todos os 20 genes foram encontrados na matriz.\n")
}

# ==========================================================
# BLOCO 3: Verificar duplicatas de símbolos gênicos
# ==========================================================

cat("\n=== VERIFICAÇÃO DE DUPLICATAS ===\n")

# Verificar duplicatas na matriz completa
duplicatas_matriz <- genes_matriz[duplicated(genes_matriz)]

if (length(duplicatas_matriz) > 0) {
  cat("Símbolos de genes duplicados na matriz completa:\n")
  print(table(duplicatas_matriz))
} else {
  cat("✓ Nenhum símbolo duplicado na matriz completa.\n")
}

# Verificar especificamente nos genes de interesse
duplicatas_interesse <- genes_encontrados[duplicated(genes_encontrados)]

if (length(duplicatas_interesse) > 0) {
  cat("\n ATENÇÃO: Genes de interesse com símbolos duplicados:\n")
  print(duplicatas_interesse)
} else {
  cat("✓ Nenhum gene de interesse possui símbolo duplicado.\n")
}

# ==========================================================
# BLOCO 4: Verificar quantas linhas cada gene ocupa
# ==========================================================

cat("\n=== VERIFICAÇÃO DE UNICIDADE ===\n")

# Contar quantas linhas cada gene ocupa
contagem_genes <- table(genes_matriz[genes_matriz %in% genes_encontrados])

cat("Contagem de linhas por gene de interesse:\n")
print(contagem_genes)

# Verificar se algum gene ocupa mais de uma linha
genes_multiplas_linhas <- names(contagem_genes[contagem_genes > 1])

if (length(genes_multiplas_linhas) > 0) {
  cat("\n ATENÇÃO: Genes que ocupam múltiplas linhas:\n")
  print(genes_multiplas_linhas)
  cat("\nIsso pode indicar:\n")
  cat("  - Isoformas diferentes\n")
  cat("  - Transcritos alternativos\n")
  cat("  - Erros de anotação\n")
  cat("\nPor enquanto, usaremos apenas a PRIMEIRA ocorrência de cada gene.\n")
  
  # Estratégia: manter apenas a primeira ocorrência
  # (poderia ser refinada em aulas futuras)
  genes_encontrados_unicos <- genes_encontrados[!duplicated(genes_encontrados)]
  
} else {
  cat("✓ Cada gene ocupa exatamente uma linha.\n")
  genes_encontrados_unicos <- genes_encontrados
}

# ==========================================================
# BLOCO 5: Extrair submatriz com genes de interesse
# ==========================================================

cat("\n=== EXTRAÇÃO DOS GENES DE INTERESSE ===\n")

# Filtrar a matriz para conter apenas genes encontrados (versão única)
expr_genes_interesse <- expr[genes_encontrados_unicos, ]

cat("Dimensões da submatriz (genes × amostras):\n")
print(dim(expr_genes_interesse))

cat("\nPrimeiras linhas da submatriz:\n")
print(head(expr_genes_interesse[, 1:5]))

cat("\nÚltimas linhas da submatriz:\n")
print(tail(expr_genes_interesse[, 1:5]))

# Verificar valores faltantes
n_na <- sum(is.na(expr_genes_interesse))
cat("\nValores faltantes (NA) na submatriz:", n_na, "\n")

# Resumo estatístico
cat("\nResumo estatístico dos valores de expressão (CPM):\n")
print(summary(as.vector(as.matrix(expr_genes_interesse))))

# Verificar distribuição por gene
cat("\nMédia de expressão (CPM) por gene:\n")
media_por_gene <- rowMeans(expr_genes_interesse)
print(sort(media_por_gene, decreasing = TRUE))

# ==========================================================
# BLOCO 6: Converter para formato longo e integrar metadados
# ==========================================================

cat("\n=== INTEGRAÇÃO COM METADADOS ===\n")

# Converter a submatriz para formato longo
expr_long <- expr_genes_interesse %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "MatrixSample",
    values_to = "CPM"
  )

cat("Dimensões da tabela em formato longo (antes do join):\n")
print(dim(expr_long))

cat("\nPrimeiras linhas:\n")
print(head(expr_long, 10))

# Integrar com metadados
expr_metadata <- expr_long %>%
  left_join(
    metadata,
    by = "MatrixSample"
  )

cat("\nDimensões após integração com metadados:\n")
print(dim(expr_metadata))

cat("\nPrimeiras linhas da tabela integrada:\n")
print(head(expr_metadata, 10))

cat("\nEstrutura da tabela:\n")
str(expr_metadata)

# Verificar valores faltantes após o join
n_na_join <- sum(is.na(expr_metadata))
cat("\nValores faltantes após join:", n_na_join, "\n")

if (n_na_join > 0) {
  cat("ATENÇÃO: Há valores faltantes após o join!\n")
  cat("Possível motivo: nomes de amostras não correspondem.\n")
  
  # Verificar quais amostras não tiveram match
  amostras_sem_metadata <- expr_long %>%
    anti_join(metadata, by = "MatrixSample") %>%
    pull(MatrixSample) %>%
    unique()
  
  cat("Amostras sem metadados:\n")
  print(amostras_sem_metadata)
}

# ==========================================================
# BLOCO 7: Validação final
# ==========================================================

cat("\n=== VALIDAÇÃO FINAL ===\n")

# Verificar genes únicos
genes_unicos <- unique(expr_metadata$Gene)
cat("Genes únicos na tabela integrada:", length(genes_unicos), "\n")

# Verificar amostras únicas
amostras_unicas <- unique(expr_metadata$MatrixSample)
cat("Amostras únicas na tabela integrada:", length(amostras_unicas), "\n")

# Verificar tratamentos
cat("\nTratamentos presentes:\n")
print(table(expr_metadata$Treatment))

# Verificar tempos
cat("\nTempos presentes:\n")
print(table(expr_metadata$Time))

# Verificar doadores
cat("\nDoadores presentes:\n")
print(table(expr_metadata$Donor))

# Verificar combinações de Treatment × Time × Donor
cat("\nCombinações de Treatment × Time × Donor:\n")
combinacoes <- expr_metadata %>%
  select(MatrixSample, Treatment, Time, Donor) %>%
  distinct() %>%
  arrange(Treatment, Time, Donor)
print(combinacoes)

cat("\nTotal de combinações esperadas: 4 × 2 × 3 = 24\n")
cat("Total de combinações encontradas:", nrow(combinacoes), "\n")

# ==========================================================
# BLOCO 8: Salvar dados processados
# ==========================================================

cat("\n=== SALVAMENTO DOS DADOS ===\n")

# Criar diretórios se não existirem
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# Salvar a submatriz de genes de interesse
write.csv(expr_genes_interesse, 
          file = "data/processed/expr_genes_interesse_matrix.csv")
cat("✓ Submatriz salva em: data/processed/expr_genes_interesse_matrix.csv\n")

# Salvar a tabela integrada em formato longo
write.csv(expr_metadata, 
          file = "data/processed/expr_genes_interesse_long.csv",
          row.names = FALSE)
cat("✓ Tabela integrada salva em: data/processed/expr_genes_interesse_long.csv\n")

# Salvar relatório de validação
validation_report <- data.frame(
  Categoria = c(
    "Genes de interesse definidos",
    "Genes encontrados na matriz",
    "Genes ausentes",
    "Genes com múltiplas linhas",
    "Genes únicos utilizados",
    "Valores faltantes (NA)",
    "Total de amostras",
    "Total de combinações (Gene × Sample)"
  ),
  Valor = c(
    length(genes_interesse),
    length(genes_encontrados),
    length(genes_ausentes),
    length(genes_multiplas_linhas),
    length(genes_encontrados_unicos),
    n_na_join,
    length(amostras_unicas),
    nrow(expr_metadata)
  )
)

write.csv(validation_report,
          file = "results/lesson09_validation_report.csv",
          row.names = FALSE)
cat("✓ Relatório de validação salvo em: results/lesson09_validation_report.csv\n")

cat("\n=== RESUMO DA VALIDAÇÃO ===\n")
print(validation_report)

# ==========================================================
# BLOCO 9: Preparar objetos para próximas aulas
# ==========================================================

cat("\n=== OBJETOS DISPONÍVEIS PARA PRÓXIMAS AULAS ===\n")

cat("\n1. expr_genes_interesse\n")
cat("   - Matriz com 20 genes × 24 amostras\n")
cat("   - Valores em CPM\n")
cat("   - Dimensões:", dim(expr_genes_interesse), "\n")

cat("\n2. expr_metadata\n")
cat("   - Tabela em formato longo\n")
cat("   - Colunas: Gene, MatrixSample, CPM, Treatment, Time, Donor\n")
cat("   - Dimensões:", dim(expr_metadata), "\n")

cat("\n3. metadata\n")
cat("   - Informações de todas as 24 amostras\n")
cat("   - Colunas: MatrixSample, Treatment, Time, Donor\n")
cat("   - Dimensões:", dim(metadata), "\n")

cat("\nTodos os objetos estão prontos para a Lesson 10!\n")

# ==========================================================
# FIM DO SCRIPT
# ==========================================================
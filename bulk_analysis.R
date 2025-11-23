library("DESeq2")
library("apeglm")
library("dplyr")
library("tidyverse")
library("tibble")
library("devtools")
library("ggplot2")
library("ggrepel")

#calculate tpm:
tpm <- function(raw_counts, gene_lengths) {
  
  x <- raw_counts*1e3 / gene_lengths
  return(t(t(x)*1e6 / colSums(x)))
  
}

tpm

# import the samples to conditions correspondence.
ivt <- read.csv("GSE116987_ivt_ms_data_exon_count_matrix_mm10_noadjust.csv", 
                 header = T, 
                 stringsAsFactors = F,
)

#read in gene names
gene <- read.csv("genename.csv", 
                    header = T, 
                    stringsAsFactors = F,
)

#read in design
ivt_design <- read.csv("GSE116987_ivt_design.csv", 
                       header = T, 
                       stringsAsFactors = F,
)

#merge to one dataset
ivt_merge <- cbind(ivt, gene)

colnames(ivt_merge)[1] <- "Gene.ID"
colnames(ivt_merge)[2] <- "2hrs.A"
colnames(ivt_merge)[3] <- "2hrs.B"
colnames(ivt_merge)[4] <- "2hrs.C"
colnames(ivt_merge)[14] <- "12days.A"
colnames(ivt_merge)[15] <- "12days.B"
colnames(ivt_merge)[16] <- "12days.C"

#select relevant time points
new_ivt <- ivt_merge %>% 
  select(Gene_name, `2hrs.A`, `2hrs.B`, `2hrs.C`, `12days.A`, `12days.B`, `12days.C`)

colnames(ref)[3] <- "Gene_name"

new_ivt <- merge(new_ivt, ref, by = "Gene_name")
gene_lengths <- as.data.frame(new_ivt$zz)
test_ivt <- new_ivt %>% 
  select(`2hrs.A`, `2hrs.B`, `2hrs.C`, `12days.A`, `12days.B`, `12days.C`)
rownames(test_ivt) <- make.names(new_ivt$Gene_name, unique=TRUE)
rownames(gene_lengths) <- make.names(new_ivt$Gene_name, unique = TRUE)

#calculate tpm values
two_a <- as.data.frame(tpm(test_ivt$`2hrs.A`, gene_lengths))
two_b <- as.data.frame(tpm(test_ivt$`2hrs.B`, gene_lengths))
two_c <- as.data.frame(tpm(test_ivt$`2hrs.C`, gene_lengths))
twelve_a <- as.data.frame(tpm(test_ivt$`12days.A`, gene_lengths))
twelve_b <- as.data.frame(tpm(test_ivt$`12days.B`, gene_lengths))
twelve_c <- as.data.frame(tpm(test_ivt$`12days.C`, gene_lengths))

ivt_tpm <- merge(two_a,two_b, by = "row.names")
rownames(ivt_tpm) <- ivt_tpm$Row.names
ivt_tpm <- ivt_tpm[,-1]
colnames(ivt_tpm)[1] <- "2hrs.A"
colnames(ivt_tpm)[2] <- "2hrs.B"
ivt_tpm <- merge(ivt_tpm, two_c, by = "row.names")
rownames(ivt_tpm) <- ivt_tpm$Row.names
ivt_tpm <- ivt_tpm[,-1]
colnames(ivt_tpm)[3] <- "2hrs.C"
ivt_tpm <- merge(ivt_tpm, twelve_a, by = "row.names")
rownames(ivt_tpm) <- ivt_tpm$Row.names
ivt_tpm <- ivt_tpm[,-1]
colnames(ivt_tpm)[4] <- "12days.A"
ivt_tpm <- merge(ivt_tpm, twelve_b, by = "row.names")
rownames(ivt_tpm) <- ivt_tpm$Row.names
ivt_tpm <- ivt_tpm[,-1]
colnames(ivt_tpm)[5] <- "12days.B"
ivt_tpm <- merge(ivt_tpm, twelve_c, by = "row.names")
rownames(ivt_tpm) <- ivt_tpm$Row.names
ivt_tpm <- ivt_tpm[,-1]
colnames(ivt_tpm)[6] <- "12days.C"

geneID <- new_ivt$Gene_name
rownames(new_ivt) <- make.names(geneID,unique = TRUE)
new_ivt <-subset(new_ivt, select = c("2hrs.A", "2hrs.B", "2hrs.C", "12days.A", "12days.B", "12days.C"))
head(new_ivt)

head(ivt_design)
rownames(ivt_design) <- ivt_design$Run
ivt_design <- subset(ivt_design, select = c("Cell.Type"))
colnames(ivt_design) <- c("celltype")
head(ivt_design)
#ivt_design$batch <- c("Aq", "Bq", "Cq", "Aa", "Ba", "Ca")

# Creation of the DESeqDataSet
DEseq2_ivt <- DESeqDataSetFromMatrix(countData = new_ivt, 
                              colData = ivt_design, 
                              design = ~ celltype)

DEseq2_ivt

DEseq2_ivt <- DESeq(DEseq2_ivt)

res_ivt <- results(DEseq2_ivt)

res_ivt

all_genes_results_ivt <- results(DEseq2_ivt, contrast = c("celltype",                      # name of the factor
                                               "qHSCs",    # name of the numerator level for fold change
                                               "aHSCs"))                          # name of the denominator level    

all.equal(res_ivt, all_genes_results_ivt)

head(all_genes_results_ivt)

diff_genes_ivt = all_genes_results_ivt %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))

head(diff_genes_ivt)

p1_ivt <- ggplot(diff_genes_ivt, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"fold change")) + 
  ylab(expression("-log"[10]*"P"))
p1_ivt +ylim(0,40)

diff_genes_ivt <- diff_genes_ivt %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(diff_genes_ivt)

p2_ivt <- ggplot(diff_genes_ivt, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"fold change")) + 
  ylab(expression("-log"[10]*"P")) +
  scale_color_manual(values = c("pink", "gray50", "lightgreen")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2_ivt +xlim (-10, +10) + ylim (0,300)

top_ivt <- 10
top_genes_ivt <- bind_rows(
  diff_genes_ivt %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top_ivt),
  diff_genes_ivt %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top_ivt)
)

p3_ivt <-  p2_ivt +
  geom_label_repel(data = top_genes_ivt,
                   mapping = aes(log2FoldChange, -log(padj,10), label = genes),
                   size = 2)
p3_ivt

options(ggrepel.max.overlaps = Inf)

#datacheck:
diff_genes_ivt %>%
  filter(genes %in% c("F3", "Adam19", "Chst11", "Acta2", "Akap12", "Mmd", "Thsd7a", "Col5a2", "Rgs5", "Fcna", "Nrxn1", "Cygb", "Apoe", "Col14a1", "Steap4"))


# import the samples to conditions correspodence.
ccl4 <- read.csv("GSE116987_ccl4_ms_data_exon_count_matrix_mm10_noadjust.csv", 
                 header = T, 
                 stringsAsFactors = F,
)

ccl4_design <- read.csv("ccl4_design.csv", 
                       header = T, 
                       stringsAsFactors = F,
)
colnames(ccl4)[1] <- "Gene.ID"
colnames(ccl4)[2] <- "0weeks.A"
colnames(ccl4)[3] <- "0weeks.B"
colnames(ccl4)[4] <- "0weeks.C"
colnames(ccl4)[5] <- "0weeks.D"
colnames(ccl4)[13] <- "8weeks.A"
colnames(ccl4)[14] <- "8weeks.B"
colnames(ccl4)[15] <- "8weeks.C"

ccl4_merge <-cbind(ccl4,gene)

new_ccl4 <- subset(ccl4_merge, select = c("Gene_name", "0weeks.A", "0weeks.B", "0weeks.C", "0weeks.D","8weeks.A", "8weeks.B", "8weeks.C"))

geneID <- new_ccl4$Gene_name
rownames(new_ccl4) <- make.names(geneID,unique = TRUE)
new_ccl4 <-subset(new_ccl4, select = c("0weeks.A", "0weeks.B", "0weeks.C", "0weeks.D", "8weeks.A", "8weeks.B", "8weeks.C"))
head(new_ccl4)

head(ccl4_design)
rownames(ccl4_design) <- ccl4_design$Run
ccl4_design <- subset(ccl4_design, select = c("Cell.Type"))
colnames(ccl4_design) <- c("celltype")
head(ccl4_design)

# Creation of the DESeqDataSet
DEseq2_ccl4 <- DESeqDataSetFromMatrix(countData = new_ccl4, 
                                     colData = ccl4_design, 
                                     design = ~ celltype)

DEseq2_ccl4

DEseq2_ccl4 <- DESeq(DEseq2_ccl4)

res_ccl4 <- results(DEseq2_ccl4)

res_ccl4

all_genes_results_ccl4 <- results(DEseq2_ccl4, contrast = c("celltype",                      # name of the factor
                                                          "qHSCs",    # name of the numerator level for fold change
                                                          "aHSCs"))                          # name of the denominator level    

all.equal(res_ccl4, all_genes_results_ccl4)

head(all_genes_results_ccl4)

diff_genes_ccl4 = all_genes_results_ccl4 %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))

head(diff_genes_ccl4)

p1_ccl4 <- ggplot(diff_genes_ccl4, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"fold change")) + 
  ylab(expression("-log"[10]*"P"))
p1_ccl4

diff_genes_ccl4 <- diff_genes_ccl4 %>% 
  mutate(
    Expression = case_when(log2FoldChange >= log(2) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(2) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(diff_genes_ccl4)

p2_ccl4 <- ggplot(diff_genes_ccl4, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"fold change")) + 
  ylab(expression("-log"[10]*"P")) +
  scale_color_manual(values = c("pink", "gray50", "lightgreen")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2_ccl4<- p2_ccl4 +xlim (-10, +10) + ylim (0,300)


top_ccl4 <- 10
top_genes_ccl4 <- bind_rows(
  diff_genes_ccl4 %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top_ccl4),
  diff_genes_ccl4 %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top_ccl4)
)

p3_ccl4 <-  p2_ccl4 +
  geom_label_repel(data = top_genes_ccl4,
                   mapping = aes(log2FoldChange, -log(padj,10), label = genes),
                   size = 2)
p3_ccl4

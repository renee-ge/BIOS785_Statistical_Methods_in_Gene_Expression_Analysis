#read in data
liver <- readRDS("humanLiver_scRNAseq_seuratObject_Wang2021.rds")

#plot umap
DimPlot(liver, reduction = "umap", label = TRUE)

#find marker genes
aHSC_markers <- FindMarkers(liver, ident.1 = "aHSCs", ident.2 = "qHSCs")
head(aHSC_markers, 10)

aHSC_markers <- aHSC_markers %>%
  arrange(avg_log2FC)
tail(aHSC_markers, 10)

#Visualize expression of top marker genes
FeaturePlot(liver_hsc, features = c("SERPINF1", "C3", "TIMP1", "COL1A2"))
FeaturePlot(liver, features = c("RGS5", "MYH11", "ADIRF", "RERGL"))

#subset the original dataset by cell type
liver_hsc <- subset(liver, ident = c("aHSCs", "qHSCs"))

#Violin plots of expression of top markers genes
VlnPlot(liver_hsc, features = c("SERPINF1", "C3", "TIMP1", "COL1A2", "DCN", "IGFBP3"))
VlnPlot(liver_hsc, features = c("RGS5", "MYH11", "ADIRF", "RERGL", "NDUFA4L2", "CRIP1"))

#Count number of significantly upregulated/downregulated genes
sum(aHSC_markers$avg_log2FC>0 & aHSC_markers$p_val<0.05)
sum(aHSC_markers$avg_log2FC<0 & aHSC_markers$p_val<0.05)

#Compute cell type proportions by sample
cell.prop <- prop.table(table(liver$Cell.type, liver$orig.ident), margin = 2)
cell.prop <- as.data.frame(cell.prop)
cell.prop <- spread(cell.prop, Var2, Freq)
colnames(cell.prop) <- c("Var1", "cirrhotic_1","cirrhotic_2", "cirrhotic_3", "healthy_1", "healthy_2", "healthy_3", "healthy_4", "healthy_5", "healthy_6", "healthy_7", "healthy_8", "healthy_9", "healthy_10")
cell.prop <- gather(cell.prop, Var2, Freq, cirrhotic_1:healthy_10)
cell.prop <- cell.prop %>%
  filter(Var1 %in% c("qHSCs", "aHSCs")) %>%
  filter(Var2 %in% c("cirrhotic_1", "cirrhotic_2", "cirrhotic_3", "cirrhotic_4", "healthy_1", "healthy_2", "healthy_3", "healthy_4"))

library(rstatix)
cell.prop.ahsc <- filter(cell.prop, Var1 == "aHSCs")
cell.prop.qhsc <- filter(cell.prop, Var1 == "qHSCs")
cell.prop.ahsc$disease <- factor(c("Cirrhotic", "Cirrhotic", "Cirrhotic", "Healthy", "Healthy", "Healthy", "Healthy"))
cell.prop.qhsc$disease <- factor(c("Cirrhotic", "Cirrhotic", "Cirrhotic", "Healthy", "Healthy", "Healthy", "Healthy"))

#Perform wilcoxon rank sum test to compare proportions between the disease groups
cell.prop.qhsc %>% wilcox_test(Freq ~ disease)

#Make boxplots of cell proportions by disease status
sc_ahsc <- cell.prop.ahsc %>%
  ggplot() + 
  geom_boxplot(aes(disease, Freq*100, fill = disease)) +
  labs(x = "Disease", y = "Percentage of Cells")

sc_qhsc <- cell.prop.qhsc %>%
  ggplot() + 
  geom_boxplot(aes(disease, Freq*100, fill = disease)) +
  labs(x = "Disease", y = "Percentage of Cells")

#save plots
ggsave("sc_ahsc_plot.png", sc_ahsc)


ggplot(cell.prop, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=1))

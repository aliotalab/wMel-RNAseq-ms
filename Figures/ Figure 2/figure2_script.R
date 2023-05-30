##Load libraries----
library(tidyverse)
library(GenomicFeatures)
library(biomaRt)
library(tximport)
library(ensembldb)
library(cowplot)
library(DESeq2)
library(wesanderson)
library(plotly)


## Create annotation and read in kallisto data (from DIY Transcriptomics)----
ensembl_metazoa <- useEnsemblGenomes(biomart = "metazoa_mart")
ensembl_AAE <- useEnsemblGenomes(biomart = "metazoa_mart",
                                 dataset = "aalvpagwg_eg_gene")
AAE.filters <- listFilters(ensembl_AAE)
Tx.AAE <- getBM(attributes=c('ensembl_transcript_id_version',
                             'external_gene_name'),
                mart = ensembl_AAE)

Tx.AAE <- as_tibble(Tx.AAE)
#we need to rename the two columns we just retreived from biomart
Tx.AAE <- dplyr::rename(Tx.AAE, target_id = ensembl_transcript_id_version, 
                        gene_name = external_gene_name)


## Day 4 Samples
targets_carcass4 <- read_tsv("~/Desktop/mosquito_rnaseq/res_noinfect/wmel_carcass4/res_carcass4_studydesign.txt")
targets_midgut4 <- read_tsv("~/Desktop/mosquito_rnaseq/res_noinfect/wmel_midgut4/res_midgut4_studydesign.txt")

path_carcass4 <- file.path("wmel_carcass4", targets_carcass4$sample_pool, "abundance.tsv")
path_midgut4 <- file.path("wmel_midgut4", targets_midgut4$sample_pool, "abundance.tsv")

Txi_transcript_carcass4 <- tximport(path_carcass4, 
                           type = "kallisto", 
                           tx2gene = Tx.AAE, 
                           txOut = TRUE, #determines whether your data represented at transcript or gene level, false = gene level
                           countsFromAbundance = "lengthScaledTPM",
                           ignoreTxVersion = TRUE)
Txi_transcript_midgut4 <- tximport(path_midgut4, 
                     type = "kallisto", 
                     tx2gene = Tx.AAE, 
                     txOut = TRUE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreAfterBar = TRUE)

sampleLabels_carcass4 <- targets_carcass4$sample_pool
sampleLabels_midgut4 <- targets_midgut4$sample_pool

## Making a matrix of raw counts for DESeq2 (Zamanian lab code)
# start with myDGEList (from DIY transcriptomics code)
gene_count_carcass4 <- round(Txi_transcript_carcass4[["counts"]])
gene_count_midgut4 <- round(Txi_transcript_midgut4[["counts"]])

class(gene_count_carcass4) # check that gene_count is a matrix

colnames(gene_count_carcass4) <- sampleLabels_carcass4 # change the column names to the samples labels from metadata in targets (read from tsv)
colnames(gene_count_midgut4) <- sampleLabels_midgut4



## Day 7 Samples
targets_carcass7 <- read_tsv("~/Desktop/mosquito_rnaseq/res_noinfect/wmel_carcass7/res_carcass7_studydesign.txt")
targets_midgut7 <- read_tsv("~/Desktop/mosquito_rnaseq/res_noinfect/wmel_midgut7/res_midgut7_studydesign.txt")

path_carcass7 <- file.path("wmel_carcass7", targets_carcass7$sample_pool, "abundance.tsv")
path_midgut7 <- file.path("wmel_midgut7", targets_midgut7$sample_pool, "abundance.tsv")

Txi_transcript_carcass7 <- tximport(path_carcass7, 
                                   type = "kallisto", 
                                   tx2gene = Tx.AAE, 
                                   txOut = TRUE, #determines whether your data represented at transcript or gene level, false = gene level
                                   countsFromAbundance = "lengthScaledTPM",
                                   ignoreTxVersion = TRUE)
Txi_transcript_midgut7 <- tximport(path_midgut7, 
                                  type = "kallisto", 
                                  tx2gene = Tx.AAE, 
                                  txOut = TRUE,
                                  countsFromAbundance = "lengthScaledTPM",
                                  ignoreAfterBar = TRUE)

sampleLabels_carcass7 <- targets_carcass7$sample_pool
sampleLabels_midgut7 <- targets_midgut7$sample_pool

## Making a matrix of raw counts for DESeq2 (Zamanian lab code)
# start with myDGEList (from DIY transcriptomics code)
gene_count_carcass7 <- round(Txi_transcript_carcass7[["counts"]])
gene_count_midgut7 <- round(Txi_transcript_midgut7[["counts"]])

class(gene_count_carcass7) # check that gene_count is a matrix

colnames(gene_count_carcass7) <- sampleLabels_carcass7 # change the column names to the samples labels from metadata in targets (read from tsv)
colnames(gene_count_midgut7) <- sampleLabels_midgut7


## Running DESeq2 and generating DEG lists (DEGs due to the presence of wMel)----
#code adapted from Zamanian Lab

## DEG list 1 - carcass4
samples_res_carcass4 <- dplyr::filter(targets_carcass4) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_c4 <- DESeqDataSetFromMatrix(countData = gene_count_carcass4, colData = samples_res_carcass4, design = ~ wmel_status)
dds_c4 <- DESeq(dds_c4) 
resultsNames(dds_c4)
res_c4 <- results(dds_c4, contrast = c("wmel_status", "wmel", "tet"))
summary(res_c4)

# save csv of full deg result
res_c4.df <- as.data.frame(res_c4) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_c4.df, file ="./deg_lists/deg.res_carcass4_full.csv", row.names=FALSE)


## DEG list 2 - midgut4
samples_res_midgut4 <- dplyr::filter(targets_midgut4) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_mg4 <- DESeqDataSetFromMatrix(countData = gene_count_midgut4, colData = samples_res_midgut4, design = ~ wmel_status)
dds_mg4 <- DESeq(dds_mg4) 
resultsNames(dds_mg4)
res_mg4 <- results(dds_mg4, contrast = c("wmel_status", "wmel", "tet"))
summary(res_mg4)

# save csv of full deg result
res_mg4.df <- as.data.frame(res_mg4) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_mg4.df, file ="./deg_lists/deg.res_midgut4_full.csv", row.names=FALSE)


## DEG list 3 - carcass7
samples_res_carcass7 <- dplyr::filter(targets_carcass7) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_c7 <- DESeqDataSetFromMatrix(countData = gene_count_carcass7, colData = samples_res_carcass7, design = ~ wmel_status)
dds_c7 <- DESeq(dds_c7) 
resultsNames(dds_c7)
res_c7 <- results(dds_c7, contrast = c("wmel_status", "wmel", "tet"))
summary(res_c7)

# save csv of full deg result
res_c7.df <- as.data.frame(res_c7) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_c7.df, file ="./deg_lists/deg.res_carcass7_full.csv", row.names=FALSE)


## DEG list 4 - midgut7
samples_res_midgut7 <- dplyr::filter(targets_midgut7) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_mg7 <- DESeqDataSetFromMatrix(countData = gene_count_midgut7, colData = samples_res_midgut7, design = ~ wmel_status)
dds_mg7 <- DESeq(dds_mg7) 
resultsNames(dds_mg7)
res_mg7 <- results(dds_mg7, contrast = c("wmel_status", "wmel", "tet"))
summary(res_mg7)

# save csv of full deg result
res_mg7.df <- as.data.frame(res_mg7) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_mg7.df, file ="./deg_lists/deg.res_midgut7_full.csv", row.names=FALSE)


## Prepping the deg lists for combining ----

deg_midgut4_full <- read_csv2('deg_lists/deg.res_midgut4_full.csv')
deg_midgut7_full <- read_csv2('deg_lists/deg.res_midgut7_full.csv')
deg_carcass4_full <- read_csv2('deg_lists/deg.res_carcass4_full.csv')
deg_carcass7_full <- read_csv2('deg_lists/deg.res_carcass7_full.csv')

# mg4 prep
transcript_id <- gsub("-R.*", "", as.character(deg_midgut4_full$gene_id))
dpf <- 4
tissue_type <- "midgut"
deg_midgut4_full <- cbind(transcript_id, deg_midgut4_full, dpf, tissue_type)
deg_midgut4_full <- deg_midgut4_full %>%
  dplyr::select(transcript_id, dpf, tissue_type, log2FoldChange, padj, gene_id)

#mg7 prep
transcript_id <- gsub("-R.*", "", as.character(deg_midgut7_full$gene_id))
dpf <- 7
tissue_type <- "midgut"
deg_midgut7_full <- cbind(transcript_id, deg_midgut7_full, dpf, tissue_type)
deg_midgut7_full <- deg_midgut7_full %>%
  dplyr::select(transcript_id, dpf, tissue_type, log2FoldChange, padj, gene_id)

#c4 prep
transcript_id <- gsub("-R.*", "", as.character(deg_carcass4_full$gene_id))
dpf <- 4
tissue_type <- "carcass"
deg_carcass4_full <- cbind(transcript_id, deg_carcass4_full, dpf, tissue_type)
deg_carcass4_full <- deg_carcass4_full %>%
  dplyr::select(transcript_id, dpf, tissue_type, log2FoldChange, padj, gene_id)

#c7 prep
transcript_id <- gsub("-R.*", "", as.character(deg_carcass7_full$gene_id))
dpf <- 7
tissue_type <- "carcass"
deg_carcass7_full <- cbind(transcript_id, deg_carcass7_full, dpf, tissue_type)
deg_carcass7_full <- deg_carcass7_full %>%
  dplyr::select(transcript_id, dpf, tissue_type, log2FoldChange, padj, gene_id)


## Combining deg lists into one large list ----
immune_list <- data.frame(read_csv2('immune_lists/immune.list.txt'))

deg <- rbind(deg_midgut4_full, deg_midgut7_full, deg_carcass4_full, deg_carcass7_full) %>%
  mutate(sig = ifelse(is.na(padj), "no", 
                      ifelse(abs(log2FoldChange) > 1 & padj < 0.01, "yes","no"))) %>%
  mutate(is_immune = ifelse(transcript_id %in% immune_list$X, "immune_gene", "n/a"))


## Panel A = volcano plots ----

# New facet label names for dpf variable
date.labs <- c("4 dpf", "7 dpf")
names(date.labs) <- c("4","7")

# New facet label names for supp variable
tis.labs <- c("Carcass", "Midgut")
names(tis.labs) <- c("carcass", "midgut")

vplot <- ggplot(deg) +
  aes(y=-log10(padj), x=log2FoldChange) +
  geom_point(data = dplyr::filter(deg, sig == "no"), colour = "grey", size=0.2, alpha = 0.1) +
  #geom_point(data = dplyr::filter(deg, is_immune == "immune_gene"), color = "red", size = 1, alpha = 1) +
  geom_point(data = dplyr::filter(deg, sig == "yes"), colour = "black", size=1, alpha = 0.2) +
  geom_point(data = dplyr::filter(deg, is_immune == "immune_gene" & sig == "yes"), color = "magenta", size = 1, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="#e7cb4f", size=0.5) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#dd3b22", size=0.5) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#5498af", size=0.5) +
  theme_bw() +
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.75)) +
  facet_grid(tissue_type ~ dpf, labeller = labeller(tissue_type = tis.labs, dpf = date.labs))
vplot

## Panel B = immune gene heatmap----
#adapted from Zamanian Lab BrugiaSpatial-ms

# combining both timepoints (counts, targets, sampleLabels)
gene_count_carcass <- cbind(gene_count_carcass4, gene_count_carcass7)
gene_count_midgut <- cbind(gene_count_midgut4, gene_count_midgut7)
gene_count4 <- cbind(gene_count_carcass4, gene_count_midgut4)
gene_count_total <- cbind(gene_count_carcass, gene_count_midgut)

targets_carcass <- rbind(targets_carcass4, targets_carcass7)
targets_midgut <- rbind(targets_midgut4, targets_midgut7)
targets4 <- rbind(targets_carcass4, targets_midgut4)
targets_total <- rbind(targets_carcass, targets_midgut)

sampleLabels_carcass <- c(sampleLabels_carcass4, sampleLabels_carcass7)
sampleLabels_midgut <- c(sampleLabels_midgut4, sampleLabels_midgut7)
sampleLabels4 <- c(sampleLabels_carcass4, sampleLabels_midgut4)
sampleLabels <- c(sampleLabels_carcass, sampleLabels_midgut)


## get vsd data--
#gene_count <- as.matrix(counts.raw)
samples <- dplyr::filter(targets4) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds <- DESeqDataSetFromMatrix(countData = gene_count4, colData = samples, design = ~ wmel_status)

#filter 1 (require x reads total across all samples)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

#filter 2 (require 5 reads in each of at least 6 samples)
keep <- rowSums(counts(dds) >= 5) >= 6
dds <- dds[keep,]
nrow(dds)

#pick vst transform
vsd <- vst(dds, blind = FALSE)

##adjust immune list
colnames(immune_list) <- c("transcript_id")

#predicted innate immune genes differentially expressed in carcass4
deg_carcass4 <- deg %>%
  dplyr::filter(dpf == "4", tissue_type == "carcass", sig == "yes") %>%
  dplyr::select(transcript_id)

immune_crossover <- dplyr::intersect(immune_list, deg_carcass4)
immune_crossover

# read in counts (vsd normalized) 
vsd_degs <- as.data.frame(assay(vsd))
vsd_degs <- vsd_degs %>%
  rownames_to_column(var = "gene_id")
vsd_degs$gene_id<- gsub('-R.*','', as.character(vsd_degs$gene_id))
vsd_degs <- as.data.frame(vsd_degs) 
vsd_degs <- vsd_degs[!duplicated(vsd_degs$gene_id), ]
#write.csv2(vsd_degs, file = "./res_noinfect_counts_GSEA.csv", row.names = FALSE)
vsd_degs <- vsd_degs %>%
  #subset(nchar(as.character(gene_id)) == 13) %>%
  dplyr::filter(gene_id %in% unique(immune_crossover$transcript_id)) %>%  # option: filter for gene list
  column_to_rownames(var = "gene_id")
counts <- vsd_degs

#join counts and sample info, widen, declare and normalize matrix
df.m <- data.matrix(counts, rownames.force = TRUE)
ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
df.m <- df.m[!ind,]
df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
counts <- df.m

#calculate gene distances and use hclust to cluster samples based on dist
geneDists <- dist(counts, method = "euclidean")
gclust.dist <- hclust(geneDists, method="ward.D2")

#get list order of clustered genes from hclust output
ord <- gclust.dist$order

#convert genecount matrix to df, re-order, and tidy (long form + metadata) for plotting
counts<- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")

#add annotations (switch labels to gene_name)

counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]), labels = c(counts$gene_id[ord]))

#long form for plotting
counts <- counts %>%
  pivot_longer(2:37, names_to = "sample", values_to = "expression") %>%
  separate(sample, c("wmel_status", "dpf", "tissue", "rep")) %>%
  group_by(gene_id, wmel_status, dpf, tissue) %>%
  mutate(expression_mean = mean(expression)) %>%
  select(-c(rep, expression)) %>%
  distinct(gene_id, wmel_status, dpf, tissue, expression_mean)
counts <- counts %>%
  unite(sample, wmel_status:tissue, remove = FALSE)
counts$wmel_status <- factor(counts$wmel_status, levels = c("WB", "TB"))
counts$dpf <- factor(counts$dpf, levels = c("4D", "7D"))
counts$tissue <- factor(counts$tissue, levels = c("MG", "C"))

c4genenames <- read_csv("./immune_lists/c4immunegenes.csv")
colnames(c4genenames) <- c("gene_id", "gene_name")
counts <- merge(counts, c4genenames)

xlabs <- c("Mosquito innate immune genes")
grouplabs <- c("COL.tet carcass", "COL.tet midgut", "COL.wMel caracss", "COL.wMel midgut")

# plot heatmap
pal <- wes_palette("Zissou1", 100, type = "continuous")
heatmap <- ggplot(counts, aes(sample, gene_name)) +
  geom_tile(aes(fill = expression_mean)) +
  scale_fill_gradientn(colours = pal, "Z-Score",
                       guide = guide_colorbar(
                         direction = "vertical",
                         title.position = "top",
                         label.position = "right"
                       )) +
  scale_y_discrete(guide = guide_axis(n.dodge = 1)) +
  scale_x_discrete(labels = grouplabs, position = "bottom") +
  xlab("") + 
  ylab("") + 
  theme_bw() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8, angle =0, hjust = 1),
    strip.background = element_rect(fill="white", color = "white"),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
    panel.spacing=unit(-0.75, "lines"),
    legend.position = "right",
    legend.margin=margin(-5,0,0,0))
#ggsave("noinfect_heatmap.pdf", plot = heatmap, device = "pdf", path = "./figures/" )
heatmap

## Panel C = temporal and tissue lfc comparisons ----

deg_temporal <- deg %>%
  dplyr::filter(sig == "yes") %>%
  dplyr::select(-transcript_id,-padj, -sig) %>%
  tidyr::pivot_wider(names_from = "dpf", values_from = "log2FoldChange") %>%
  dplyr::rename(log2FC_4 = `4`, log2FC_7 = `7`) %>%
  dplyr::mutate(pattern = ifelse(log2FC_4 > 0 & log2FC_7 > 0, "constitutively increased",
                                 ifelse(log2FC_4 < 0 & log2FC_7 < 0, "constitutively decreased", 
                                        ifelse(log2FC_4 >0 & log2FC_7 < 0, "up down",
                                               ifelse(log2FC_4 < 0 & log2FC_7 >0, "down up", "other"))))) %>%
  dplyr::mutate(slope = log2FC_7-log2FC_4) %>%
  tidyr::pivot_longer(log2FC_4:log2FC_7, names_to = "log2FC_day", values_to = "log2FoldChange")

  
deg_temporal_list <- deg_temporal %>% dplyr::filter(pattern %in% c("up_up","down_down"))
deg_temporal_list <- deg_temporal_list$gene_id

deg_static <- deg %>%
  dplyr::filter(sig == "yes") %>%
  dplyr::select(-transcript_id,-padj, -sig) %>%
  tidyr::pivot_wider(names_from = "dpf", values_from = "log2FoldChange") %>%
  dplyr::rename(log2FC_4 = `4`, log2FC_7 = `7`) %>%
  dplyr::mutate(pattern = ifelse(log2FC_4 > 0 & log2FC_7 > 0, "up",
                                 ifelse(log2FC_4 > 0 & is.na(log2FC_7), "up",
                                        ifelse(is.na(log2FC_4) & log2FC_7 > 0, "up",
                                               ifelse(log2FC_4 < 0 & log2FC_7 < 0, "down", 
                                                      ifelse(log2FC_4 < 0 & is.na(log2FC_7), "down",
                                                             ifelse(is.na(log2FC_4) & log2FC_7 < 0, "down","other"))))))) %>%
  tidyr::pivot_longer(log2FC_4:log2FC_7, names_to = "log2FC_day", values_to = "log2FoldChange") %>%
  dplyr::filter(!gene_id %in% deg_temporal_list)

temp_text <- data.frame(
  label = c("Muscle lim protein", "GPRMTH6"),
  tissue_type = c("carcass", "carcass"),
  x = c("log2FC_7", "log2FC_7"),
  y = c(-7, -9),
  size = c(5, 5)
  
)

x_labs <- c("4", "7")
#temporal plot
temp_plot <- ggplot(data = dplyr::filter(deg_temporal, pattern %in% c("constitutively increased","constitutively decreased")), text = paste("Symbol:", gene_id)) +
  aes(y = log2FoldChange, x = factor(as.character(log2FC_day))) + # text = paste("Symbol:", gene_id)) +
  #geom_jitter(aes(color = pattern), alpha = 0.75, width = 0.05) +
  geom_point(aes(color = pattern), alpha = 0.5) +
  scale_color_manual(values = c("#5498af", "#dd3b22")) +
  geom_text(data = temp_text,
            mapping = aes(x = x, y = y, label = label), size = 2, nudge_x = 0.28) +
  geom_line(aes(x = log2FC_day, y = log2FoldChange, group = gene_id, color = pattern, alpha = abs(slope)), show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.01, 1)) +
  scale_x_discrete(labels = x_labs) +
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.75)) +
  facet_grid(. ~ tissue_type, labeller = labeller(tissue_type = tis.labs)) +
  labs(x = "days post-feeding (dpf)") +
  theme_bw() +
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.75), legend.position = "bottom")
temp_plot
#ggplotly(temp_plot)


## put all the panels together----
figure_2.6 <- plot_grid(vplot, heatmap, temp_plot, labels = c("A.","B.", "C."), nrow = 1, axis = "tblr")
figure_2.6

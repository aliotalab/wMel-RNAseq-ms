##Load libraries----
library(tidyverse)
library(GenomicFeatures)
library(biomaRt)
library(tximport)
library(ensembldb)
library(cowplot)
library(DESeq2)
library(pheatmap)
library(hexbin)
library(ggVennDiagram)
library(wesanderson)
library(here)
library(topGO)
library(conflicted)
library(Rgraphviz)

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


## Day 4
targets_carcass4 <- read_tsv("~/Desktop/mosquito_rnaseq/sus_zikvinfect/tet_carcass4/sus_carcass4_studydesign.txt")
targets_midgut4 <- read_tsv("~/Desktop/mosquito_rnaseq/sus_zikvinfect/tet_midgut4/sus_midgut4_studydesign.txt")

path_carcass4 <- file.path("tet_carcass4", targets_carcass4$sample_pool, "abundance.tsv")
path_midgut4 <- file.path("tet_midgut4", targets_midgut4$sample_pool, "abundance.tsv")

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

## Making a matrix of raw counts for Mostafa's code
# start with myDGEList (from DIY transcriptomics code)
gene_count_carcass4 <- round(Txi_transcript_carcass4[["counts"]])
gene_count_midgut4 <- round(Txi_transcript_midgut4[["counts"]])

class(gene_count_carcass4) # check that gene_count is a matrix

colnames(gene_count_carcass4) <- sampleLabels_carcass4 # change the column names to the samples labels from metadata in targets (read from tsv)
colnames(gene_count_midgut4) <- sampleLabels_midgut4



## Day 7
targets_carcass7 <- read_tsv("~/Desktop/mosquito_rnaseq/sus_zikvinfect/tet_carcass7/sus_carcass7_studydesign.txt")
targets_midgut7 <- read_tsv("~/Desktop/mosquito_rnaseq/sus_zikvinfect/tet_midgut7/sus_midgut7_studydesign.txt")

path_carcass7 <- file.path("tet_carcass7", targets_carcass7$sample_pool, "abundance.tsv")
path_midgut7 <- file.path("tet_midgut7", targets_midgut7$sample_pool, "abundance.tsv")

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

## Making a matrix of raw counts for Mostafa's code
# start with myDGEList (from DIY transcriptomics code)
gene_count_carcass7 <- round(Txi_transcript_carcass7[["counts"]])
gene_count_midgut7 <- round(Txi_transcript_midgut7[["counts"]])

class(gene_count_carcass7) # check that gene_count is a matrix

colnames(gene_count_carcass7) <- sampleLabels_carcass7 # change the column names to the samples labels from metadata in targets (read from tsv)
colnames(gene_count_midgut7) <- sampleLabels_midgut7


## Running DESeq2 and generating DEG lists (degs due to the presence of wMel)----

## DEG list 1 - carcass4
samples_res_carcass4 <- dplyr::filter(targets_carcass4) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_c4 <- DESeqDataSetFromMatrix(countData = gene_count_carcass4, colData = samples_res_carcass4, design = ~ infect_status)
dds_c4 <- DESeq(dds_c4) 
resultsNames(dds_c4)
res_c4 <- results(dds_c4, contrast = c("infect_status", "ZIKV", "blood"))
summary(res_c4)

# save csv of full deg result
res_c4.df <- as.data.frame(res_c4) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_c4.df, file ="./deg_lists/deg.carcass4_full.csv", row.names=FALSE)


## DEG list 2 - midgut4
samples_res_midgut4 <- dplyr::filter(targets_midgut4) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_mg4 <- DESeqDataSetFromMatrix(countData = gene_count_midgut4, colData = samples_res_midgut4, design = ~ infect_status)
dds_mg4 <- DESeq(dds_mg4) 
resultsNames(dds_mg4)
res_mg4 <- results(dds_mg4, contrast = c("infect_status", "ZIKV", "blood"))
summary(res_mg4)

# save csv of full deg result
res_mg4.df <- as.data.frame(res_mg4) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_mg4.df, file ="./deg_lists/deg.midgut4_full.csv", row.names=FALSE)


## DEG list 3 - carcass7
samples_res_carcass7 <- dplyr::filter(targets_carcass7) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_c7 <- DESeqDataSetFromMatrix(countData = gene_count_carcass7, colData = samples_res_carcass7, design = ~ infect_status)
dds_c7 <- DESeq(dds_c7) 
resultsNames(dds_c7)
res_c7 <- results(dds_c7, contrast = c("infect_status", "ZIKV", "blood"))
summary(res_c7)

# save csv of full deg result
res_c7.df <- as.data.frame(res_c7) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_c7.df, file ="./deg_lists/deg.carcass7_full.csv", row.names=FALSE)


## DEG list 4 - midgut7
samples_res_midgut7 <- dplyr::filter(targets_midgut7) %>%
  mutate(rownames = sample_pool) %>%
  column_to_rownames(var = "rownames")
dds_mg7 <- DESeqDataSetFromMatrix(countData = gene_count_midgut7, colData = samples_res_midgut7, design = ~ infect_status)
dds_mg7 <- DESeq(dds_mg7) 
resultsNames(dds_mg7)
res_mg7 <- results(dds_mg7, contrast = c("infect_status", "ZIKV", "blood"))
summary(res_mg7)

# save csv of full deg result
res_mg7.df <- as.data.frame(res_mg7) %>%
  rownames_to_column(var = "gene_id")
write.csv2(res_mg7.df, file ="./deg_lists/deg.midgut7_full.csv", row.names=FALSE)


## Prepping the deg lists for combining ----

deg_midgut4_full <- read_csv2('deg_lists/deg.midgut4_full.csv')
deg_midgut7_full <- read_csv2('deg_lists/deg.midgut7_full.csv')
deg_carcass4_full <- read_csv2('deg_lists/deg.carcass4_full.csv')
deg_carcass7_full <- read_csv2('deg_lists/deg.carcass7_full.csv')

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
#deg_full <- rbind(deg_midgut4_full, deg_midgut7_full, deg_carcass4_full, deg_carcass7_full)

deg <- rbind(deg_midgut4_full, deg_midgut7_full, deg_carcass4_full, deg_carcass7_full) %>%
  mutate(sig = ifelse(is.na(padj), "no", 
                      ifelse(abs(log2FoldChange) > 1 & padj < 0.01, "yes","no")))

## Panel A = volcano plots ----

#plotting MG4 and MG7, no carcass

deg_mg <- dplyr::filter(deg, tissue_type == "midgut")

# New facet label names for dpf variable
date.labs <- c("4 dpf", "7 dpf")
names(date.labs) <- c("4","7")

# New facet label names for supp variable
tis.labs <- c("Carcass", "Midgut")
names(tis.labs) <- c("carcass", "midgut")

vplot <- ggplot(deg_mg) +
  aes(y=-log10(padj), x = log2FoldChange, text = paste("Symbol:", gene_id)) + 
  geom_point(data = dplyr::filter(deg, sig == "no"), colour = "grey", size=0.5, alpha = 0.1) +
  geom_point(data = dplyr::filter(deg, sig == "yes"), colour = "black", size=2, alpha = 0.25) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="#e7cb4f", size=0.75) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#5498af", size=0.75) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#dd3b22", size=0.75) +
  theme_bw() +
  facet_grid(.~dpf, labeller = labeller(dpf = date.labs)) +
  theme(strip.background = element_rect(color = "black", fill = "white", linewidth = 0.75))
  #facet_grid(tissue_type ~ dpf, labeller = labeller(tissue_type = tis.labs, dpf = date.labs))
vplot
ggplotly(vplot)


##GO-Figure prep for MG7 degs ----

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")


# Create GO ID --> gene/transcript ID mappings ----------------------------
# Note: this only needs to be performed if you have reason to believe that
#       the mappings have been updated, otherwise you can use the previously
#       generated files

### VectorBase species

# aedes
path <- 'https://vectorbase.org/common/downloads/Current_Release/AaegyptiLVP_AGWG/gaf/VectorBase-61_AaegyptiLVP_AGWG_GO.gaf'
aedes_go <- read_tsv(path,
                     comment = '!',
                     trim_ws = TRUE,
                     col_names = c('db', 'gene_id', 'name', 'null', 'go_id', 'db:reference', 'evidence_code', 'null', 'aspect', 'db_object_name', 'null', 'db_object_type', 'taxon', 'date', 'assigned_by')) %>%
  janitor::remove_empty('cols') 

aedes_go_out <- dplyr::select(aedes_go, gene_id, go_id) %>%
  group_by(gene_id) %>%
  distinct() %>%
  summarise(go_ids = list(go_id)) %>%
  mutate(go_ids = str_c(go_ids)) %>%
  mutate(go_ids = str_remove_all(go_ids, "c\\(")) %>%
  mutate(go_ids = str_remove_all(go_ids, '\\"')) %>%
  mutate(go_ids = str_remove_all(go_ids, '\\)'))

write.table(aedes_go_out, here('aedes_gene_go.txt'), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


# GO enrichment -----------------------------------------------------------
# load in the text file that has GO IDs matched to gene or transcript IDs
# Note: this must be a text file (can't be a data frame already in memory)
#       with a very specific structure; you can regenerate these files with
#       the commands listed in the first section

# load in the mappings (this must be a file, can't be a df)
gene_go <- topGO::readMappings(here('aedes_gene_go.txt'))

# read the df for later perusal
genes_go <- read_tsv(here('aedes_gene_go.txt'), col_names = c('gene_id', 'go_ids')) %>%
  separate_rows(go_ids, sep = ', ') %>%
  rename(go_id = go_ids)

# get the list of possible gene_ids
gene_ids <- names(gene_go)

# read in your list of genes/transcripts of interest
interest_genes <- read_csv2("./deg_lists/deg.midgut7_up.csv")
interest_genes <- interest_genes[,"gene_id"]
interest_genes$gene_id <- gsub("-R.*","",as.character(interest_genes$gene_id))
interest_go <- left_join(interest_genes, genes_go, by = "gene_id")
go_summary <- group_by(interest_go, go_id) %>%
  summarize(n = n()) %>%
  filter(!is.na(go_id))

# the final data.frame needs to have one column with all the transcript_ids
# and a second column denoting whether or not it is a transcript of interest
final_genes <- distinct(select(genes_go, gene_id)) %>%
  mutate(interest = case_when(
    gene_id %in% interest_genes$gene_id ~ 1,
    TRUE ~ 0
  )) 

# get the interest column as a factor
final_genes_tg <- as.factor(final_genes$interest)

# convert to a factor with names
names(final_genes_tg) <- final_genes$gene_id

# create the topGOdata objects
# MF == molecular function
# BP == biological process
# CC == cellular component
go_data_mf <- new("topGOdata", ontology = "MF", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)
go_data_bp <- new("topGOdata", ontology = "BP", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)
go_data_cc <- new("topGOdata", ontology = "CC", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)

# create the statistical test
fisher_test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

all_go_data <- tibble(top_go_object = c(go_data_mf, go_data_bp, go_data_cc)) %>% # make a tibble of all 3 topGOobjects
  mutate(test = map(top_go_object, getSigGroups, fisher_test)) %>% # run the fisher test on each topGOobject
  mutate(result = map2(top_go_object, test, GenTable, ranksOf = "classic", topNodes = 10)) %>% # extract significant GO IDs to a df
  mutate(result = map2(top_go_object, result, ~mutate(.y, class = .x@ontology))) # add the GO class as a column for each nested df

# view the graph for a given class
showSigOfNodes(all_go_data[[1]][[2]], all_go_data[[2]][[2]]@score, firstSigNodes = 5, useInfo = 'all')

plot_data <- select(all_go_data, result) %>%
  unnest(cols = c(result)) %>%
  janitor::clean_names() %>%
  rename(result = result1) %>%
  mutate(class = case_when(
    class == 'BP' ~ 'Biological Process',
    class == 'MF' ~ 'Molecular Function',
    class == 'CC' ~ 'Cellular Component'
  ))
plot_data <- as.data.frame(plot_data)
write_tsv(plot_data, file = "~/Desktop/mosquito_rnaseq/sus_zikvinfect/GOFigure_inputs/topGO_MG7up.tsv")




library(tidyverse)
library(here)
library(topGO)
library(conflicted)
library(Rgraphviz)

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
interest_genes <- read_csv2(here("./deg_lists/deg.carcass7_down.csv")) ## change file for each tissue type, dpf, and expression direction"./deg.tissuetypedpf.csv"
colnames(interest_genes) <- c("gene_id")
interest_genes <- interest_genes[,"gene_id"]
#interest_genes$gene_id <- gsub("-R.*","",as.character(interest_genes$gene_id))
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
plot_data <- plot_data[,c("go_id","term","annotated","significant","expected","result")]

write_tsv(plot_data, file = "~/Desktop/mosquito_rnaseq/res_v_sus_zikvinfect/GOFigure_inputs/topGO_C7down.tsv") ## change file for each tissue type, dpf, and expression direction"./deg.tissuetypedpf.csv"




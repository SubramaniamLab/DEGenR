suppressMessages(library (shiny))
#library(markdown)
#library(ggplot2)
suppressMessages(library(DT))
suppressMessages(library(limma))
suppressMessages(library(shinydashboard))
suppressMessages(library(edgeR))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
suppressMessages(library(fgsea))
suppressMessages(library(data.table))
#library(Glimma)
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(msigdbr))
#library(shinyBS)
suppressMessages(library(topconfects))
suppressMessages(library(enrichR))
suppressMessages(library(viper))
suppressMessages(library(dorothea))
suppressMessages(library("tmod"))
suppressMessages(library(GEOquery))
suppressMessages(library(pheatmap))
suppressMessages(library(umap))
suppressMessages(library("maptools"))
suppressMessages(library(dplyr))
suppressMessages(library(recount))
suppressMessages(library(SummarizedExperiment))
#suppressMessages(library("tximport"))
#suppressMessages(library("scaterlegacy"))
#suppressMessages(library("fastqcr"))


Hs.GO <- msigdbr(species = "Homo sapiens", category = "C5")
Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
Hs.Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
#Hs.motif <- msigdbr(species = "Homo sapiens", category = "C3")
Hs.GOCC.full <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
# Gene Ontology: Molecular Function (Full)
Hs.GOMF.full <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
# Human Phenotype Ontology
Hs.HPO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO")
Hs.Biocarta <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
# KEGG
Hs.KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
# Pathway Interaction Database
Hs.PID <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
# WikiPathways
Hs.WikiPathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
# MSigDB Chemical and Genetic Perturbations 
Hs.CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
# MSigDB Computational Genesets
Hs.Comp <- msigdbr(species = "Homo sapiens", category = "C4")
# MSigDB Oncogenic Signature Genesets
Hs.Oncogenic <- msigdbr(species = "Homo sapiens", category = "C6")
# MSigDB Immunologic signature Genesets
Hs.Immune <- msigdbr(species = "Homo sapiens", category = "C7")
# MSigDB Cell Types
Hs.CellType <- msigdbr(species = "Homo sapiens", category = "C8")
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")


#Hs.GOBP.Entrez <-  split(x = as.numeric(Hs.GOBP$entrez_gene), f = Hs.GOBP$gs_name)
Hs.GOBP.Symbol <- split(x = Hs.GOBP$gene_symbol, f = Hs.GOBP$gs_name)
#Hs.Hallmark.Entrez <- Hs.Hallmark %>% split(x = .$entrez_gene, f = .$gs_name)
Hs.Hallmark.Symbol <- Hs.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
#Hs.Reactome.Entrez <- Hs.Reactome %>% split(x = .$entrez_gene, f = .$gs_name)
Hs.Reactome.Symbol <- Hs.Reactome %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.GO.Symbol <- Hs.GO %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.CellType.Symbol <- Hs.CellType %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Immune.Symbol <- Hs.Immune %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Oncogenic.Symbol <- Hs.Oncogenic %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Comp.Symbol <- Hs.Comp %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.CGP.Symbol <- Hs.CGP %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.WikiPathways.Symbol <- Hs.WikiPathways %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.PID.Symbol <- Hs.PID %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.KEGG.Symbol <- Hs.KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Biocarta.Symbol <- Hs.Biocarta %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.HPO.Symbol <- Hs.HPO %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.GOMF.Symbol <- Hs.GOMF.full %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.GOCC.full.Symbol <- Hs.GOCC.full %>% split(x = .$gene_symbol, f = .$gs_name)
							
# For Enrichr 

Enrichrdbs <- listEnrichrDbs()
dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
            "ENCODE_TF_ChIP-seq_2015", 
            "ChEA_2016",
            "TRANSFAC_and_JASPAR_PWMs", 
            "TargetScan_microRNA",
            "ARCHS4_TFs_Coexp",
            "TRRUST_Transcription_Factors_2019",
            "TargetScan_microRNA_2017",
            "miRTarBase_2017")
		 
dbs_ontology <- c("GO_Biological_Process_2018",
                  "GO_Biological_Process_2017b", 
                  "GO_Molecular_Function_2018",
                  "GO_Molecular_Function_2017b", 
                  "GO_Cellular_Component_2018",
                  "GO_Cellular_Component_2017b",
                  "MSigDB_Hallmark_2020",
                  "Reactome_2016",
                  "BioCarta_2016",
                  "KEGG_2019_Human",
                  "Panther_2016",
                  "WikiPathways_2019_Human",
                  "BioPlanet_2019")
		 

# For DoRothEA TF analysis
#
data(dorothea_hs, package = "dorothea")
#Subset the regulon for stringency 
regulon_a = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A"))
regulon_b = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))
regulon_c = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))
regulon_d = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C", "D"))
regulon_e = dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C", "D", "E"))


###fgsea TF
# ENCODE-ChEA Consensus
Hs.ECC <- qusage::read.gmt(file.path("./data", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt"))
# ChEA 2016 (Enrichr)
Hs.ChEA2016 <- qusage::read.gmt(file.path("./data", "ChEA_2016.gmt"))
# ENCODE 2015 (Enrichr)
Hs.ENCODE <- qusage::read.gmt(file.path("./data", "ENCODE_TF_ChIP-seq_2015.gmt"))
# ReMap ChIP-seq 2018 Human
Hs.ReMap <- qusage::read.gmt(file.path("./data", "ReMap_ChIP-seq.gmt"))
# TRRUST 2019 Human
Hs.TRRUST <- qusage::read.gmt(file.path("./data", "TRRUST_Transcription_Factors_2019.gmt"))
# ChEA3 Literature ChIP-Seq
Hs.Literature <- qusage::read.gmt(file.path("./data", "Literature_ChIP-seq.gmt"))
# TRANSFAC/JASPAR PWMs (Enrichr)
Hs.TRANSFACJASPAR <- qusage::read.gmt(file.path("./data", "TRANSFAC_and_JASPAR_PWMs.gmt"))
# Gene Transcription Regulation Database (GTRD v20.06) 
Hs.GTRD <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
Hs.GTRD.sel = split(x = Hs.GTRD$gene_symbol, f = Hs.GTRD$gs_name)
# MSigDB Legacy TF targets
Hs.TFLegacy <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy")
Hs.TFLegacy.sel = split(x = Hs.TFLegacy$gene_symbol, f = Hs.TFLegacy$gs_name)

# miRTarBase 2017 
Hs.miR <- qusage::read.gmt(file.path("./data", "miRTarBase_2017.gmt"))
# miRNA TargetScan 2017
Hs.miRTargetScan <- qusage::read.gmt(file.path("./data", "TargetScan_microRNA_2017.gmt"))
# miRDB v6.0
Hs.miRNA <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "MIR:MIRDB")
Hs.miRNA.sel = split(x = Hs.miRNA$gene_symbol, f = Hs.miRNA$gs_name)


##cerno

msig <- tmodImportMSigDB(file.path("./data", "msigdb_v7.2.xml"))
#Select the Hallmark gene sets (you can selectwhichever geneset db you would like)
Hallmark.sel <- msig$MODULES$Category == "H"
# Gene Ontology Biological Process (MSigDB Filtered)
GOBP.sel <- msig$MODULES$Subcategory == "C5_GO:BP"
# Gene Ontology Biological Process (Full)
GOBPfull.sel <- msig$MODULES$Subcategory == "GO:BP"
#reactome
Reactome.sel <- (msig$MODULES$Subcategory == "CP:REACTOME") | (msig$MODULES$Subcategory == "C2_CP:REACTOME")
# Biocarta
BioCARTA.sel <- msig$MODULES$Subcategory == "C2_CP:BIOCARTA"

# Gene Ontology Molecular Function (MSigDB Filtered)
GOMF.sel <- msig$MODULES$Subcategory == "C5_GO:MF"
# Gene Ontology Molecular Function (Full)
GOMFfull.sel <- msig$MODULES$Subcategory == "GO:MF"
# Gene Ontology Cellular Compartment (MSigDB Filtered)
GOCC.sel <- msig$MODULES$Subcategory == "C5_GO:CC"
# Gene Ontology Cellular Compartment (Full)
GOCCfull.sel <- msig$MODULES$Subcategory == "C5_GO:CC"
# Human Phenotype Ontology
HPO.sel <- msig$MODULES$Subcategory == "HPO"

# KEGG
Kegg.sel <- msig$MODULES$Subcategory == "CP:KEGG"
# Pathway Interaction Database
PID.sel <- msig$MODULES$Subcategory == "CP:PID"
#Wikipathways
Wiki.sel <- msig$MODULES$Subcategory == "CP:WIKIPATHWAYS"
# MSigDB Chemical and Genetic Perturbations 
CGP.sel <- msig$MODULES$Subcategory == "C2_CGP"
# MSigDB Computational Genesets
CM.sel <- msig$MODULES$Subcategory == "CM"
# MSigDB Oncogenic Signature Genesets
Onc.sel <- msig$MODULES$Category == "C6"
# MSigDB Immunologic signature Genesets
Imm.sel <- msig$MODULES$Category == "C7"
# MSigDB Cell Types
CellType.sel <- msig$MODULES$Category == "C8"


data(tmod)
ECC.gmt <-  qusage::read.gmt(file.path("./data", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt"))
ECC.MODULES <- as.data.frame(read.csv(file.path("./data", "ECCTFlist_MODULES.csv")))
ECC.MODULES2GENES <- as.list(ECC.gmt)
rownames(ECC.MODULES) <- ECC.MODULES$ID


# Make ReMap ChIP-Seq mset for use with tmod
ReMap.gmt <-  qusage::read.gmt(file.path("./data", "ReMap_ChIP-seq.gmt"))
ReMap.MODULES <- as.data.frame(read.csv(file.path("./data", "ReMapTFlist_MODULES.csv")))
ReMap.MODULES2GENES <- as.list(ReMap.gmt)
rownames(ReMap.MODULES) <- ReMap.MODULES$ID

# TRRUST 2019 Human
TRRUST.gmt <-  qusage::read.gmt(file.path("./data", "TRRUST_Transcription_Factors_2019.gmt"))
TRRUST.MODULES <- as.data.frame(read.csv(file.path("./data", "TRRUST_Transcription_Factors_2019_MODULES.csv")))
TRRUST.MODULES2GENES <- as.list(TRRUST.gmt)
rownames(TRRUST.MODULES) <- TRRUST.MODULES$ID


# TRANSFAC/JASPAR PWMs (Enrichr)
TFJ.gmt <-  qusage::read.gmt(file.path("./data", "TRANSFAC_and_JASPAR_PWMs.gmt"))
TFJ.MODULES <- as.data.frame(read.csv(file.path("./data", "TRANSFAC_and_JASPAR_PWMs_MODULES.csv")))
TFJ.MODULES2GENES <- as.list(TFJ.gmt)
rownames(TFJ.MODULES) <- TFJ.MODULES$ID


# Gene Transcription Regulation Database (GTRD v20.06) 
GTRD.sel <- msig$MODULES$Subcategory == "TFT:GTRD"
# miRTarBase 2017
miRTarBase.gmt <-  qusage::read.gmt(file.path("./data", "miRTarBase_2017.gmt"))
miRTarBase.MODULES <- as.data.frame(read.csv(file.path("./data", "miRTarBase_MODULES.csv")))
miRTarBase.MODULES2GENES <- as.list(miRTarBase.gmt)
rownames(miRTarBase.MODULES) <- miRTarBase.MODULES$ID

# use msetmiRTarBase
# miRDB v6.0
MIRDB.sel <- msig$MODULES$Subcategory == "MIR:MIRDB"


# create transcript to gene (t2g)
# martGRCh38.99 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                                  # dataset = "hsapiens_gene_ensembl",
                                  # host = 'jan2020.archive.ensembl.org',
                                  # path="/biomart/martservice")
# GRCh38.99t2g<- biomaRt::getBM(attributes = c("ensembl_transcript_id_version",
                                             # "ensembl_gene_id",
                                             # "external_gene_name"),
                              # mart = martGRCh38.99)
# GRCh38.99t2g <- dplyr::rename(GRCh38.99t2g, 
                              # TXNAME = ensembl_transcript_id_version,
                              # ENSEMBL = ensembl_gene_id,
                              # Symbol = external_gene_name)



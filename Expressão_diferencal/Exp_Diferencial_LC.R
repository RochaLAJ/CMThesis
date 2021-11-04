#-----------PACOTES NECESSÁRIOS-----------#
require(NOISeq)
require(dplyr)
require(biomaRt)

#-----------CONSTANTES--------------------#
design <- as.factor(rep(c("1", "2"), each = 4))
alpha <- 0.95
myfactors <- data.frame(
    Thyroid    = c("1", "1", "1","1", "2", "2", "2", '2'),
    ThyroidRun = c("1_1", "1_1", "1_1","1_1","2_2", "2_2", "2_2","2_2"),
    Run        = c(rep("R2", 4), rep("R2", 4))
                       )
interest <- c("ensembl_gene_id", "hgnc_symbol", "transcript_length", 
              "gene_biotype", "entrezgene_id","percentage_gene_gc_content",
              "version","description")

#----------SISTEMA DE ANOTAÇÃO------------#
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

#-----------IMPORTAÇÃO & FILTRO-----------#
data <- read.csv(file.choose(), header=T, row.names=1)
row.names(data) <- gsub("\\.[0-9]*$", "", row.names(data))
depth <- list(sum(data$CTL1), sum(data$CTL2), sum(data$CTL3), sum(data$CTL4),
              sum(data$LNC1), sum(data$LNC2), sum(data$LNC3), sum(data$LNC4))
bol <- filterByExpr(data, design)
data <- data[bol,]
rm(bol)

#----------------NOISeq-------------------#
data <- tmm(data)
mydata <- readData(data = data, factors = myfactors)
mynoiseq <- noiseqbio(
  mydata,
  norm   = "n",
  factor = "Thyroid",
  filter = 3,
  a0per  = 0.8,
  depth  = c(depth[[1]],depth[[2]],depth[[3]],depth[[4]],depth[[5]],depth[[6]],depth[[7]],depth[[8]])
)
GDE <- degenes(mynoiseq, q = alpha)
GDE <- GDE %>% filter(log2FC> 0.4| log2FC< -0.4)

#------------APLICANDO ANOTAÇÃO------------#
holder <- getBM(attributes = interest,
                values     = row.names(GDE), 
                uniqueRows = TRUE,
                verbose    = TRUE,
                mart       = mart)
holder <- holder[!duplicated(holder$ensembl_gene_id), ]
holder <- data.frame(holder, row.names=1)
GDE <- merge(holder, GDE, by="row.names")
GDE <- data.frame(GDE, row.names=1)
GDE_LNC <- merge(GDE, data, by="row.names")
rm(holder)
rm(mart)
rm(depth)
rm(GDE)
#Coruja

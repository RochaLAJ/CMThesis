#-------------PACOTES NECESSÁRIOS-----------#
require(TissueEnrich)
require(biomaRt)

#-------------CONSTANTES--------------------#
interest <- c("ensembl_gene_id", "hgnc_symbol", "transcript_length", 
"gene_biotype", "entrezgene_id","percentage_gene_gc_content",
"version","description")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

#-------------IMPORTAÇÃO & FILTRO------------#
data <- read.csv("ExpListLNC.csv", header=T, row.names=1)
data <- data[!duplicated(data$hgnc_symbol), ]
background <- read.csv("counts_local.csv", header=T, row.names=1)
row.names(background) <- gsub("\\.[0-9]*$", "", row.names(background))

#------------ENRIQUECIMENTO-----------------#
holder <- getBM(attributes = interest,
                values     = row.names(background),
                uniqueRows = TRUE,
                verbose    = TRUE,
                mart       = mart)
holder <- holder[!duplicated(holder$ensembl_gene_id), ]
background <- data.frame(holder, row.names=1)
background <- background[!duplicated(background$hgnc_symbol), ]

gs <- GeneSet(geneIds    = data$hgnc_symbol,
              organism   = "Homo Sapiens",
              geneIdType = SymbolIdentifier())

gs_background <- GeneSet(geneIds    = background$hgnc_symbol,
                         organism   = "Homo Sapiens",
                         geneIdType = SymbolIdentifier())

gse <- teEnrichment(inputGenes             = gs,
                    rnaSeqDataset          = 1,
                    tissueSpecificGeneType = 1,
                    multiHypoCorrection    = TRUE,
                    backgroundGenes        = gs_background)

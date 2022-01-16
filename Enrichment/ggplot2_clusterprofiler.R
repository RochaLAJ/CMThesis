#--------------------PACKAGES LOADING
require(dplyr)
require(ggplot2)
require(org.Hs.eg.db)
require(clusterProfiler)
require(RColorBrewer)



#--------------------DATA LOADING
all <- list.files(pattern='.csv')
dataLoader <- lapply(all, read.csv)

lc <- dataLoader[[1]]
mdt <- dataLoader[[2]]

all_tissues <- list.files(pattern='.csv')
barplot <- lapply(all_tissues, read.csv)


#---------------------TRANSLATE
lc <- lc %>% dplyr::select(log2FC, gene_biotype, entrezgene_id) %>% 
  mutate(Biotype = case_when(gene_biotype == 0 ~ 'RNA Codificante de Proteína',
                             gene_biotype == 'antisense' ~ 'Antisenso',
                             gene_biotype == 'lincRNA' ~ 'RNA Intergênico Não Codificante',
                             gene_biotype == 'Mt_rRNA' ~ 'RNA Ribossômico Mitocondrial',
                             gene_biotype == 'nonsense_mediated_decay' ~ 'RNA Degradado por Mutação Não Senso',
                             gene_biotype == 'processed_pseudogene' ~ 'Pseudogene Processado',
                             gene_biotype == 'processed_transcript' ~ 'Transcrito Processado',
                             gene_biotype == 'protein_coding' ~  'Codificante de Proteína',
                             gene_biotype == 'retained_intron' ~ 'Intron Retido',
                             gene_biotype == 'scaRNA' ~ 'RNA Pequeno de Corpos Cajais',
                             gene_biotype == 'snoRNA' ~ 'RNA Pequeno Nucleolar',
                             gene_biotype == 'TEC' ~ 'Tipo Ainda Será Confirmado [To Be Confirmed Experimentally]',
                             gene_biotype == 'transcribed_unprocessed_pseudogene' ~ 'Pseudogene Transcrito Não Processado',
                             gene_biotype == 'unprocessed_pseudogene' ~ 'Pseudogene Não Processado',
                             gene_biotype == 'sense_intronic' ~ 'Senso Intrônico',
                             gene_biotype == 'lncRNA' ~ 'RNA Longo Não Codificante',
                             gene_biotype == 'transcribed_processed_pseudogene' ~ 'Pseudogene Transcrito Processado',
                             gene_biotype == 'transcribed_unitary_pseudogene' ~ 'Pseudogene Transcrito Unitário'
  )
  )

mdt <- mdt %>% dplyr::select(log2FC, gene_biotype, entrezgene_id) %>% 
  mutate(Biotype = case_when(gene_biotype == 0 ~ 'RNA Codificante de Proteína',
                             gene_biotype == 'antisense' ~ 'Antisenso',
                             gene_biotype == 'lincRNA' ~ 'RNA Intergênico Não Codificante',
                             gene_biotype == 'Mt_rRNA' ~ 'RNA Ribossômico Mitocondrial',
                             gene_biotype == 'nonsense_mediated_decay' ~ 'RNA Degradado por Mutação Não Senso',
                             gene_biotype == 'processed_pseudogene' ~ 'Pseudogene Processado',
                             gene_biotype == 'processed_transcript' ~ 'Transcrito Processado',
                             gene_biotype == 'protein_coding' ~  'Codificante de Proteína',
                             gene_biotype == 'retained_intron' ~ 'Intron Retido',
                             gene_biotype == 'scaRNA' ~ 'RNA Pequeno de Corpos Cajais',
                             gene_biotype == 'snoRNA' ~ 'RNA Pequeno Nucleolar',
                             gene_biotype == 'TEC' ~ 'Tipo Ainda Será Confirmado [To Be Confirmed Experimentally]',
                             gene_biotype == 'transcribed_unprocessed_pseudogene' ~ 'Pseudogene Transcrito Não Processado',
                             gene_biotype == 'unprocessed_pseudogene' ~ 'Pseudogene Não Processado',
                             gene_biotype == 'sense_intronic' ~ 'Senso Intrônico',
                             gene_biotype == 'lncRNA' ~ 'RNA Longo Não Codificante',
                             gene_biotype == 'transcribed_processed_pseudogene' ~ 'Pseudogene Transcrito Processado',
                             gene_biotype == 'transcribed_unitary_pseudogene' ~ 'Pseudogene Transcrito Unitário',
                             gene_biotype == 'TR_C_gene' ~ 'TR_C_gene',
                             gene_biotype == 'polymorphic_pseudogene' ~ 'Pseudogene Polimórfico'
  )
  )

#--------------------DATA CLEANING
mdt$gene_biotype <- NULL
lc$gene_biotype <- NULL

#---------------------COLORS PALETTE
mypalette<-brewer.pal(10,"Paired")
#-----------------GENERATING ENRICHMENT DATA
ggo_lc <- enrichGO(gene     = as.character(lc$entrezgene_id),
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",
                  pvalueCutoff = 0.01)

ggo_mdt <- enrichGO(gene     = as.character(mdt$entrezgene_id),
                   OrgDb    = org.Hs.eg.db,
                   ont      = "BP",
                   pvalueCutoff = 0.01)

#-----------------GENERATING PLOTS
ggo_lc = dotplot(ggo_lc, showCategory = 7, title="Pricpais Procesos Biológicos pelo GO \nCTLxLC")
ggo_mdt = dotplot(ggo_mdt, showCategory = 7, title="CTLxMDT")


lc_polar <- ggplot(lc, aes(x = factor(1), fill = factor(Biotype))) +
  geom_bar(width = 1) +
 coord_polar(theta = "y") + theme_bw() + ggtitle("Proporção de Biotipos de RNAs \nCTLxLC")

mdt_polar <- ggplot(mdt, aes(x = factor(1), fill = factor(Biotype))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y") + theme_bw() + ggtitle("CTLxMDT") +   scale_fill_manual(values=mypalette)


lc_bar <- ggplot(barplot[[1]], aes(x=reorder(Tecidos_LC, Tecidos_LC, function(x)-length(x)),fill=Tecidos_LC)) +
  geom_bar()  + theme_bw() + labs(title="Número de Genes por Tecido \nCTLxLC", x ="Tecidos", y = "Quantidade Absoluta")


mdt_bar <- ggplot(barplot[[2]], aes(x=reorder(Tecido_MDT, Tecido_MDT, function(x)-length(x)),fill=Tecido_MDT)) +
  geom_bar()  + theme_bw() + labs(title="CTLxMDT", x ="Tecidos", y = "Quantidade Absoluta") + scale_fill_manual(values=mypalette)

  
#---------------PLOT PAGE
ggarrange(lc_polar, ggo_lc, 
          lc_bar, mdt_polar, 
          ggo_mdt, mdt_bar,
          labels = "AUTO", ncol=3, nrow=2)


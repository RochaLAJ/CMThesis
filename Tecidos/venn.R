require(VennDiagram)

data <- read.csv(file.choose(), header=T)

venn.diagram(
  x = list(data$mdt, data$lc),
  category.names = c("MDT Tecidos" , "LC Tecidos"),
  filename = '#14_venn_diagramm.png',
  height= 4500, width= 4500,
  print.mode = c('raw','percent'),
  hyper.test = TRUE,
  total.population = 86,
  sub = 'n = 86',
  sub.cex = 2.5,
  cat.cex = 2.5,
  cat.pos = 2,
  col = 4:5,
  fill = 3:4,
  alpha = 0.3,
  cex = 2,
  output = TRUE
)

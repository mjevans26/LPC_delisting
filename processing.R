nwe <- read.csv("NewWellEval.csv", header = TRUE, sep = ",")
nwe$X <- NULL
nwe$X.1 <- NULL
nwe <- nwe[,c(3,6,7,8,9,10)]
nwe$RATIO <- 1-(nwe$POS)/nwe$NEG

cde <- read.csv("ChDetEval.csv", header = TRUE, sep = ",")
cde$X <- NULL
cde$X.1 <- NULL
cde <- cde[,c(3,6,7,8,9,10)]
cde$RATIO <- 1-(cde$POS)/cde$NEG

plot_ly(data = cde, x = ~1-NEG, y = ~POS, color = ~ndvi, type = "scatter", mode = "markers")

pcc(cde[,3:6], cde[,1], nboot = 0)
src(cde[,3:6], cde[,2], nboot = 0)
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

plot_ly(data = cde, x = ~1-NEG, y = ~POS, color = ~ndvi,
        type = "scatter", mode = "markers",
        text = ~ paste('CV:', cv, "RCV:", rcv, "NDVI:", ndvi, "NSDI:", ndsi),
        hoverinfo = text)

pcc(cde[,3:6], cde[,1], nboot = 0)
src(cde[,3:6], cde[,2], nboot = 0)

rawcd <- read.csv("SRange_chngData.csv", header = TRUE, sep = ",")

plot_ly()%>%
  add_trace(data = rawcd[rawcd$Index == 1,], x = ~ mode, y = ~ndsi_z, type = "box")%>%
  add_trace(data = rawcd[rawcd$Index == 0,], x = ~ mode, y = ~ndsi_z, type = "box")
  layout(barmode = "group")
  
rawcd_stats <- read.table("SRange_ChangeStats.txt", header = TRUE, sep = "\t")

rawcd$cv_z <- (rawcd$cv - rawcd_stats$Min[rawcd_stats$Metric == "cv"])/rawcd_stats$SD[rawcd_stats$Metric == "cv"]
rawcd$rcv_z <- (rawcd$rcvmax - rawcd_stats$Min[rawcd_stats$Metric == "rcv"])/rawcd_stats$SD[rawcd_stats$Metric == "rcv"]
rawcd$ndvi_z <- (rawcd$ndvi - rawcd_stats$Mean[rawcd_stats$Metric == "ndvi"])/rawcd_stats$SD[rawcd_stats$Metric == "ndvi"]
rawcd$ndsi_z <- (rawcd$ndsi - rawcd_stats$Mean[rawcd_stats$Metric == "ndsi"])/rawcd_stats$SD[rawcd_stats$Metric == "ndsi"]

ct <- ctree(Index~cv_z+rcv_z+ndvi_z+ndsi_z+mode, data = rawcd[rawcd$mode == 82,], controls = ctree_control(mincriterion = 0.95))

auc <- roc(response = goddamn$group, predictor = unlist(predict(ct, newdata = goddamn, type = "prob"))[seq(1,nrow(goddamn)*2,2)])

scrub <- lda(Index ~ cv_z + rcv_z + ndvi_z + ndsi_z, data = rawcd)
scrub_pred <- predict(scrub, type = "prob")
ldahist(data = scrub_pred$x[,1], g = scrub_pred$class)
scrub_out <- data.frame("predict" = scrub_pred$class, "LD1" = scrub_pred$x[,1], "p0" = scrub_pred$posterior[,1], "p1" = scrub_pred$posterior[,2], "observed" = rawcd$Index)
scrub_out$predict <- as.numeric(as.character(scrub_out$predict))
scrub_roc <- roc(response = scrub_out$observed, predictor = scrub_out$p1)
plot_ly(data = as.data.frame(scrub_roc[2:4]), y = ~sensitivities, x = ~ specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)

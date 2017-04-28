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
rawcdN <- read.csv("NRange_chngData.csv", header = TRUE, sep = ",")

plot_ly()%>%
  add_trace(data = rawcd[rawcd$Index == 1,], x = ~ mode, y = ~ndsi_z, type = "box")%>%
  add_trace(data = rawcd[rawcd$Index == 0,], x = ~ mode, y = ~ndsi_z, type = "box")
  layout(barmode = "group")
  
rawcd_stats <- read.table("SRange_ChangeStats.txt", header = TRUE, sep = "\t")
rawcd_statsN <- read.table("NRange_ChangeStats.txt", header = TRUE, sep = "\t")

rawcdN$cv_z <- (rawcdN$cv - rawcd_statsN$Min[rawcd_statsN$Metric == "cv"])/rawcd_statsN$SD[rawcd_statsN$Metric == "cv"]
rawcdN$rcv_z <- (rawcdN$rcvmax - rawcd_statsN$Min[rawcd_statsN$Metric == "rcv"])/rawcd_statsN$SD[rawcd_statsN$Metric == "rcv"]
rawcdN$ndvi_z <- (rawcdN$ndvi - rawcd_statsN$Mean[rawcd_statsN$Metric == "ndvi"])/rawcd_statsN$SD[rawcd_statsN$Metric == "ndvi"]
rawcdN$ndsi_z <- (rawcdN$ndsi - rawcd_statsN$Mean[rawcd_statsN$Metric == "ndsi"])/rawcd_statsN$SD[rawcd_statsN$Metric == "ndsi"]

ct <- ctree(Index~cv_z+rcv_z+ndvi_z+ndsi_z+mode, data = rawcd[rawcd$mode == 82,], controls = ctree_control(mincriterion = 0.95))

auc <- roc(response = goddamn$group, predictor = unlist(predict(ct, newdata = goddamn, type = "prob"))[seq(1,nrow(goddamn)*2,2)])

scrub <- lda(Index ~ cv_z + rcv_z + ndvi_z + ndsi_z, data = rawcd)
scrub_pred <- predict(scrub, type = "prob")
ldahist(data = scrub_pred$x[,1], g = scrub_pred$class)
scrub_out <- data.frame("predict" = scrub_pred$class, "LD1" = scrub_pred$x[,1], "p0" = scrub_pred$posterior[,1], "p1" = scrub_pred$posterior[,2], "observed" = rawcd$Index)
scrub_out$predict <- as.numeric(as.character(scrub_out$predict))
scrub_roc <- roc(response = scrub_out$observed, predictor = scrub_out$LD1)
plot_ly(data = as.data.frame(scrub_roc[2:4]), y = ~sensitivities, x = ~ 1- specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)

N <- lda(Index ~ cv_z + rcv_z + ndvi_z + ndsi_z, data = rawcdN)
N_pred <- predict(N, type = "prob")
ldahist(data = N_pred$x[,1], g = N_pred$class)
N_out <- data.frame("predict" = N_pred$class, "LD1" = N_pred$x[,1], "p0" = N_pred$posterior[,1], "p1" = N_pred$posterior[,2], "observed" = rawcdN$Index)
N_out$predict <- as.numeric(as.character(N_out$predict))
N_roc <- roc(response = N_out$observed, predictor = N_out$LD1)
plot_ly(data = as.data.frame(N_roc[2:4]), y = ~sensitivities, x = ~ 1 - specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)

LSag <- read.csv(file = "KS_ag_raw.csv", header = TRUE)

library(plotly)
library(pROC)
library(MASS)
library(pscl)
library(dplyr)
library(EnvStats)
library(fitdistrplus)
#inverse gamma function to estimate variance given sample size
invgammapdf <- function(x, n, v){
  alpha <- n/2
  beta <- n*v/2
  return((beta^alpha/factorial(alpha-1)*(x^(-alpha-1))*exp(-1*beta/x)))
}

(beta^alpha)/factorial(alpha-1)*(0.5^(-alpha-1))*exp(-1*beta/0.5)

ndvi_raw <- read.csv("cdl_ndvivalues.csv", header = TRUE, sep = ",")
ndvi_raw$system.index <- as.character(ndvi_raw$system.index)
ndvi_raw$feat <- 0
for(i in 1:nrow(ndvi_raw)){
  if(ndvi_raw$Class[i] == 0){
    ndvi_raw$feat[i] <- paste(str_split(ndvi_raw$system.index[i], "_")[[1]][1:2], collapse = "")
  }else{
    ndvi_raw$feat[i] <- paste(str_split(ndvi_raw$system.index[i], "_")[[1]][1:3], collapse = "")
  }
}

new2 <- group_by(ndvi_raw, feat)%>%summarise(Class = first(Class),
                                             NDVI_max = median(NDVI_max),
                                             NDVI_var = median(NDVI_variance),
                                             NDVI_mean = median(NDVI_mean),
                                             NDVI_disp = median(disp),
                                             count = median(NDVI_count))
-3.33734 + ( "RATIO"*(1-"CONVEX")*-0.2876616) +(94.4642983*("Shape_Area"/(POWER("Shape_Leng",2)/4*3.14))/0.000247)< -0.4
##Shape and Landcover data from Arc, after filtering for well detection output polygons
## 1 > x < 15 acres, and falling in scrub or grassland.
NWells <- read.csv("NWell_ID_thresh.csv", header = TRUE, sep = ",")
SWells <- read.csv("SWell_ID_thresh.csv", header = TRUE, sep = ",")

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


##LDA DOESNT PRINT INTERCEPT!!
scrub <- lda(Index ~ cv_z + rcv_z + ndvi_z + ndsi_z, data = rawcd)
scrub_pred <- predict(scrub, type = "prob")
ldahist(data = scrub_pred$x[,1], g = scrub_pred$class)
scrub_out <- data.frame("predict" = scrub_pred$class, "LD1" = scrub_pred$x[,1], "p0" = scrub_pred$posterior[,1], "p1" = scrub_pred$posterior[,2], "observed" = rawcd$Index)
scrub_out$predict <- as.numeric(as.character(scrub_out$predict))
scrub_roc <- roc(response = scrub_out$observed, predictor = scrub_out$LD1)
plot_ly(data = as.data.frame(scrub_roc[2:4]), y = ~sensitivities, x = ~ specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)

N <- lda(Index ~ cv_z + rcv_z + ndvi_z + ndsi_z, data = rawcdN)
N_pred <- predict(N, type = "prob")
ldahist(data = N_pred$x[,1], g = N_pred$class)
N_out <- data.frame("predict" = N_pred$class, "LD1" = N_pred$x[,1], "p0" = N_pred$posterior[,1], "p1" = N_pred$posterior[,2], "observed" = rawcdN$Index)
N_out$predict <- as.numeric(as.character(N_out$predict))
N_roc <- roc(response = N_out$observed, predictor = N_out$LD1)
plot_ly(data = as.data.frame(shape_roc[2:4]), y = ~sensitivities, x = ~ 1 - specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)

LSag <- read.csv(file = "KS_ag_raw.csv", header = TRUE)
LSag$ID <- as.character(LSag$ID)
LSag$Feat <- substr(LSag$ID, 3, nchar(LSag$ID))
LSag$Feat <- as.factor(LSag$Feat)

LSchange <- LSag[LSag$Season == 2 & duplicated(LSag$Feat), 9:11] - LSag[LSag$Season == 1, 9:11]
LSchange$Feat <- LSag$Feat[LSag$Season == 1]
LSchange$Class <- LSag$Class[LSag$Season == 1]
LSchange$cv <- sum((LSag[LSag$Season == 2 & duplicated(LSag$Feat), 2:7] - LSag[LSag$Season == 1, 2:7])^2)
for (i in LSag$Feat[LSag$Season == 2 & duplicated(LSag$Feat)]){
  LSchange$cv[LSchange$Feat == i] <- sum((LSag[LSag$Season == 2 & LSag$Feat == i,2:7] - LSag[LSag$Season == 1 & LSag$Feat == i, 2:7])^2)
    b <- LSag[LSag$Season == 2 & LSag$Feat == i, 2:7]
    a <- LSag[LSag$Season == 1 & LSag$Feat == i, 2:7]
    m <- c(max(b[1],a[1]), max(b[2],a[2]), max(b[3],a[3]), max(b[4],a[4]), max(b[5],a[5]), max(b[6],a[6]))
  LSchange$rcvmax[LSchange$Feat == i] <- sum(((b-a)/m)^2)
}

LSchange$cv_z <- (LSchange$cv - LSchange_stats$Min[LSchange_stats$Metric == "cv"])/LSchange_stats$SD[LSchange_stats$Metric == "cv"]
LSchange$rcv_z <- (LSchange$rcvmax - LSchange_stats$Min[LSchange_stats$Metric == "rcv"])/LSchange_stats$SD[LSchange_stats$Metric == "rcv"]
LSchange$ndvi_z <- (LSchange$ndvi - LSchange_stats$Mean[LSchange_stats$Metric == "ndvi"])/LSchange_stats$SD[LSchange_stats$Metric == "ndvi"]
LSchange$ndsi_z <- (LSchange$ndsi - LSchange_stats$Mean[LSchange_stats$Metric == "ndsi"])/LSchange_stats$SD[LSchange_stats$Metric == "ndsi"]

##Intercept is -0.957484
change <- lda(Class ~ cv_z + rcv_z + ndvi_z + ndsi_z, data = LSchange)
change_pred <- predict(change, type = "prob")
ldahist(data = change_pred$x[,1], g = change_pred$class)
change_out <- data.frame("predict" = change_pred$class, "LD1" = change_pred$x[,1], "p0" = change_pred$posterior[,1], "p1" = change_pred$posterior[,2], "observed" = LSchange$Class)
change_out$predict <- as.numeric(as.character(change_out$predict))
change_roc <- roc(response = change_out$observed, predictor = change_out$LD1)
plot_ly(data = as.data.frame(change_roc[2:4]), y = ~sensitivities, x = ~ 1 - specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)


AG_on <- lda(Class ~ cv + + ndvi + ndsi, data = LSag[LSag$Season == 2,])
AGon_pred <- predict(AG_on, type = "prob")
ldahist(data = AGon_pred$x[,1], g = AGon_pred$class)
AGon_out <- data.frame("predict" = AGon_pred$class, "LD1" = AGon_pred$x[,1], "p0" = AGon_pred$posterior[,1], "p1" = AGon_pred$posterior[,2], "observed" = LSag[LSag$Season == 2]$Class)
AGon_out$predict <- as.numeric(as.character(AGon_out$predict))
AGon_roc <- roc(response = AGon_out$observed, predictor = AGon_out$LD1)
plot_ly(data = as.data.frame(AGon_roc[2:4]), y = ~sensitivities, x = ~ 1 - specificities, type = "scatter", text = ~paste("Thresh:", thresholds), hoverinfo = text)

AG_off <- lda(Class ~ cv + + ndvi + ndsi, data = LSag[LSag$Season == 1,])
AGon_pred <- predict(AG_off, type = "prob")
ldahist(data = AGoff_pred$x[,1], g = AGoff_pred$class)
AGoff_out <- data.frame("predict" = AGoff_pred$class, "LD1" = AGoff_pred$x[,1], "p0" = AGoff_pred$posterior[,1], "p1" = AGoff_pred$posterior[,2], "observed" = LSag[LSag$Season == 1]$Class)
AGoff_out$predict <- as.numeric(as.character(AGoff_out$predict))
AGoff_roc <- roc(response = AGoff_out$observed, predictor = AGoff_out$LD1)

intercept <- function(input, model, prediction, i){
  prediction$LD1[i] -
  sum(input$cv_z[i]*model$scaling["cv_z",1], input$rcv_z[i]*model$scaling["rcv_z",1],
      input$ndvi_z[i]*model$scaling["ndvi_z",1], input$ndsi_z[i]*model$scaling["ndsi_z",1])
}

distdiff <- function(d1, d2, iters){
  v <- rep(NA, iters)
  for(i in 1:iters){
    v[i] <- d1-d2
  }
  d <- density(v)
  pdf <- cbind(d$x, d$y/sum(d$y))
  colnames(pdf) <- c("x", "y")
  return(pdf)
}

#difference in best fit parametric distributions *done on variances
pdiff11 <- distdiff(rweibull(1, fit1$estimate[1], fit1$estimate[2]),
                     rweibull(1, fit1$estimate[1], fit1$estimate[2]),
                     5000)

pdiff121 <- distdiff(rweibull(1, fit1$estimate[1], fit1$estimate[2]),
                     rgevd(1, fit21$estimate[1], fit21$estimate[2], fit21$estimate[3]),
                     5000)

pdiff122 <- distdiff(rweibull(1, fit1$estimate[1], fit1$estimate[2]),
                     rgevd(1, fit22$estimate[1], fit22$estimate[2], fit22$estimate[3]),
                     5000)

pdiff123 <- distdiff(rweibull(1, fit1$estimate[1], fit1$estimate[2]),
                     rgevd(1, fit23$estimate[1], fit23$estimate[2], fit23$estimate[3]),
                     5000)

pdiff11var <- distdiff(rgamma(1, fit1var$estimate[1], fit1var$estimate[2]),
                        rgamma(1, fit1var$estimate[1], fit1var$estimate[2]),
                        5000)

pdiff31var <- distdiff(rgamma(1, fit23var$estimate[1], fit23var$estimate[2]),
                     rgamma(1, fit21var$estimate[1], fit21var$estimate[2]),
                     5000)

pdiff21var <- distdiff(rgamma(1, fit22var$estimate[1], fit22var$estimate[2]),
                     rgamma(1, fit21var$estimate[1], fit21var$estimate[2]),
                     5000)

pdiff23var <- distdiff(rgamma(1, fit22var$estimate[1], fit22var$estimate[2]),
                     rgamma(1, fit23var$estimate[1], fit23var$estimate[2]),
                     5000)

p1 <- plot_ly(ndvi_raw, type = "histogram", alpha = 0.6)%>%
  add_trace(x = ~NDVI[cropland == 1 & NDVI_count > 5], name = 'corn')%>%
  add_trace(x = ~NDVI[cropland == 4& NDVI_count > 5], name = 'sorghum')%>%
  add_trace(x = ~NDVI[cropland == 24 & NDVI_count > 5], name = 'winter wheat')%>%
  add_trace(x = ~NDVI[cropland == 61 & NDVI_count > 5], name = 'barren')%>%
  add_trace(x = ~NDVI[cropland == 36 & NDVI_count > 5], name = 'alfalfa')%>%
  add_trace(x = ~NDVI[cropland == 152|cropland == 176 & NDVI_count > 5], name = 'habitat')%>%
  layout(barmode = "overlay")
p2 <- plot_ly(ndvi_raw, type = "histogram", alpha = 0.6)%>%
  add_trace(x = ~NDVI_disp[cropland == 1 & NDVI_count > 5], name = 'corn')%>%
  add_trace(x = ~NDVI_disp[cropland == 4& NDVI_count > 5], name = 'sorghum')%>%
  add_trace(x = ~NDVI_disp[cropland == 24 & NDVI_count > 5], name = 'winter wheat')%>%
  add_trace(x = ~NDVI_disp[cropland == 61 & NDVI_count > 5], name = 'barren')%>%
  add_trace(x = ~NDVI_disp[cropland == 36 & NDVI_count > 5], name = 'alfalfa')%>%
  add_trace(x = ~NDVI_disp[cropland == 152|cropland == 176 & NDVI_count > 5], name = 'habitat')%>%
  layout(barmode = "overlay")
subplot(p1, p2, nrows = 2)

#difference in observed distributions of dispersion values
d_hab_hab <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                          |ndvi_raw$cropland == 176
                                          & ndvi_raw$NDVI_count > 5],
                              5000, replace = TRUE),
                         sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                                   |ndvi_raw$cropland == 176
                                                   & ndvi_raw$NDVI_count > 5],
                                5000, replace = TRUE)
                       )

dvar_hab_corn <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                    |ndvi_raw$cropland == 176
                                                    & ndvi_raw$NDVI_count > 5],
                                 5000, replace = TRUE),
                          sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                                    & ndvi_raw$NDVI_count > 5],
                                 5000, replace = TRUE)
                        )

dvar_hab_wheat <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 24
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
                           )

dvar_hab_gum <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 4
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
                           )

d_hab_bare <- distdiff2(sample(ndvi_raw$NDVI[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI[ndvi_raw$cropland == 61
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
                           )

dvar_hab_alph <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                     |ndvi_raw$cropland == 176
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE),
                           sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 36
                                                     & ndvi_raw$NDVI_count > 5],
                                  5000, replace = TRUE)
                           )

dvar_crp_crp <- distdiff2(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                            |ndvi_raw$cropland == 24
                                           |ndvi_raw$cropland == 36
                                           |ndvi_raw$cropland == 4
                                            & ndvi_raw$NDVI_count > 5],
                              20000, replace = TRUE),
                       sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                            |ndvi_raw$cropland == 24
                                            |ndvi_raw$cropland == 36
                                            |ndvi_raw$cropland == 4
                                            & ndvi_raw$NDVI_count > 5],
                              20000, replace = TRUE)
                       )

dvar_ag_crp <- distdiff(sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 152
                                                   |ndvi_raw$cropland == 176
                                                   & ndvi_raw$NDVI_count > 5],
                                1, replace = TRUE),
                         sample(ndvi_raw$NDVI_disp[ndvi_raw$cropland == 1
                                                   |ndvi_raw$cropland == 24
                                                   |ndvi_raw$cropland == 36
                                                   |ndvi_raw$cropland == 4
                                                   & ndvi_raw$NDVI_count > 5],
                                1, replace = TRUE),
                         20000)


p3<-plot_ly(type = "scatter", mode = "lines")%>%
  add_trace(x = ~d_crp_crp[,1], y = ~cumsum(d_crp_crp[,2]), line = list(color = "black"),
            name = "crop")%>%
  add_trace(x = ~d_hab_hab[,1], y = ~d_hab_hab[,2]/max(d_hab_hab[,2])*0.5, line = list(color = "black"),
            name = "hab")%>%
  add_trace(x = ~d_hab_corn[,1], y = ~1-cumsum(d_hab_corn[,2]), line = list(color = "yellow"),
            name = 'corn')%>%
  add_trace(x = ~d_hab_wheat[,1], y = ~1-cumsum(d_hab_wheat[,2]), line = list(color = "brown"),
            name = 'wheat')%>%
  add_trace(x = ~d_hab_bare[,1], y = ~cumsum(d_hab_bare[,2]), line = list(color = "grey"),
            name = 'bare')%>%
  add_trace(x = ~d_hab_alph[,1], y = ~1-cumsum(d_hab_alph[,2]), line = list(color = "green"),
            name = 'alphalpha')%>%
  add_trace(x = ~d_hab_gum[,1], y = ~1-cumsum(d_hab_gum[,2]), line = list(color = "purple"),
            name = 'sorghum')

p4<-plot_ly(type = "scatter", mode = "lines")%>%
  add_trace(x = ~dvar_crp_crp[,1], y = ~cumsum(dvar_crp_crp[,2]), line = list(color = "black"),
            name = "crop")%>%
  add_trace(x = ~dvar_hab_hab[,1], y = ~dvar_hab_hab[,2]/max(dvar_hab_hab[,2])*0.5, line = list(color = "black"),
            name = "hab")%>%
  add_trace(x = ~dvar_hab_corn[,1], y = ~1-cumsum(dvar_hab_corn[,2]), line = list(color = "yellow"),
            name = 'corn')%>%
  add_trace(x = ~dvar_hab_wheat[,1], y = ~1-cumsum(dvar_hab_wheat[,2]), line = list(color = "brown"),
            name = 'wheat')%>%
  add_trace(x = ~dvar_hab_bare[,1], y = ~cumsum(dvar_hab_bare[,2]), line = list(color = "grey"),
            name = 'bare')%>%
  add_trace(x = ~dvar_hab_alph[,1], y = ~1-cumsum(dvar_hab_alph[,2]), line = list(color = "green"),
            name = 'alphalpha')%>%
  add_trace(x = ~dvar_hab_gum[,1], y = ~1-cumsum(dvar_hab_gum[,2]),line = list(color = "purple"),
            name = 'sorghum')

p5<-plot_ly(type = "scatter", mode = "lines")%>%
  add_trace(x = ~d_crp_crp[,1], y = ~d_crp_crp[,2], line = list(color = "black"),
            name = "crop")%>%
  add_trace(x = ~d_hab_hab[,1], y = ~d_hab_hab[,2], line = list(color = "black"),
            name = "hab")%>%
  add_trace(x = ~d_hab_corn[,1], y = ~d_hab_corn[,2], line = list(color = "yellow"),
            name = 'corn')%>%
  add_trace(x = ~d_hab_wheat[,1], y = ~d_hab_wheat[,2], line = list(color = "brown"),
            name = 'wheat')%>%
  add_trace(x = ~d_hab_bare[,1], y = ~d_hab_bare[,2], line = list(color = "grey"),
            name = 'bare')%>%
  add_trace(x = ~d_hab_alph[,1], y = ~d_hab_alph[,2], line = list(color = "green"),
            name = 'alphalpha')%>%
  add_trace(x = ~d_hab_gum[,1], y = ~d_hab_gum[,2], line = list(color = "purple"),
            name = 'sorghum')

p6<-plot_ly(type = "scatter", mode = "lines")%>%
  add_trace(x = ~dvar_crp_crp[,1], y = ~dvar_crp_crp[,2], line = list(color = "black"),
            name = "crop")%>%
  add_trace(x = ~dvar_hab_hab[,1], y = ~dvar_hab_hab[,2], line = list(color = "black"),
            name = "hab")%>%
  add_trace(x = ~dvar_hab_corn[,1], y = ~dvar_hab_corn[,2], line = list(color = "yellow"),
            name = 'corn')%>%
  add_trace(x = ~dvar_hab_wheat[,1], y = ~dvar_hab_wheat[,2], line = list(color = "brown"),
            name = 'wheat')%>%
  add_trace(x = ~dvar_hab_bare[,1], y = ~dvar_hab_bare[,2], line = list(color = "grey"),
            name = 'bare')%>%
  add_trace(x = ~dvar_hab_alph[,1], y = ~dvar_hab_alph[,2], line = list(color = "green"),
            name = 'alphalpha')%>%
  add_trace(x = ~dvar_hab_gum[,1], y = ~dvar_hab_gum[,2],line = list(color = "purple"),
            name = 'sorghum')

subplot(p3, p4, nrows = 2)

distdiff2 <-function(d1, d2){
  v <- d1-d2
  d <- density(v)
  pdf <- cbind(d$x, d$y/sum(d$y))
  colnames(pdf) <- c("x", "y")
  return(pdf)
}

auc_diff <- function(d1, d2, d3, d4, d5, val){
max((sum(d2[d2[,1] > val, 2]))*(1-sum(d1[d1[,1] > val, 2])),
     (sum(d3[d3[,1] > val, 2]))*(1-sum(d1[d1[,1] > val, 2])),
     (sum(d4[d4[,1] > val, 2]))*(1-sum(d1[d1[,1] > val, 2])),
     (sum(d5[d5[,1] > val, 2]))*(1-sum(d1[d1[,1] > val, 2])))
}
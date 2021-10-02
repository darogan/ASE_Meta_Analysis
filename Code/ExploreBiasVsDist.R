

library("readxl")
library("reshape2")
library("ggplot2")
library("ggridges")
library("tidyverse")
library("ggdist")
library("gghalves")
library("ggExtra")
library("ggpubr")
library("lattice")
require("tigerstats")



baseDir <- "/Users/rhamilto/Documents/CTR-Manuscripts/2021-Carol_Paper/"

setwd(baseDir)


message("+-------------------------------------------------------------------------------")
message("+ Read in TADS Excel Table ")
message("+-------------------------------------------------------------------------------")

xlsx <- read_excel("TADs_for_Russell.xlsx", col_names=T, skip=1)
str(xlsx)

xlsx.df <- as.data.frame(xlsx)
xlsx.df.m <- melt(xlsx.df)
str(xlsx.df.m)


message("+-------------------------------------------------------------------------------")
message("+ Read in TADS Excel Table ")
message("+-------------------------------------------------------------------------------")




ggplot(xlsx.df.m, aes(y = variable, x = value)) +
  geom_density_ridges( 
    #stat = "binline",
    fill="blue", alpha=.25, size=0.5, rel_min_height = 0.025, na.rm=T,
    panel_scaling=TRUE,
    jittered_points = T, quantile_lines=T, quantiles=2, scale = 0.75, 
    vline_size = 1, vline_color = "red",
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE, height=0.15)
  #  position = position_points_jitter(width = 0.2, height = 0)
    ) +
  #coord_flip() +
  theme_bw() 
  



ggplot(xlsx.df.m, aes(y = variable, x = value)) +
  geom_boxplot(fill = "grey92") +
  geom_point(
    size = 2,
    alpha = .3,
    position = position_jitter(  seed = 1, width = .1  )
  )


xlsx.df.m.na <- na.omit(xlsx.df.m) 


# Basic histogram plot with mean line and marginal rug
gghistogram(xlsx.df.m.na, x = "value", bins = 50, group="variable",  #facet.by="variable",
            fill = "variable", color = "variable", add_density=TRUE, 
            add = "mean", rug = TRUE) +
  coord_flip() +
  facet_wrap(~variable, nrow=2) +
  theme_bw() +
  theme(legend.position="none")
















set.seed(1)
d <- data.frame(a = rnorm(100), b = rnorm(100, 1), c = rnorm(100, 2),
                d = rnorm(100, 3), e = rnorm(100, 4))



m111survey

d <- xlsx.df.m
colnames(d) <- c("A", "B", "C", "D", "E", 
                 "F", "G", "H", "I", "J")

densityplot(~ value|variable, data = d, groups = variable,
  layout=c(5,2), pch = 25)


library("utils")
#install.packages("remotes")
library("remotes")
#remotes::install_github("chihlinwei/ddecay")
library("ddecay")

## https://rdrr.io/github/chihlinwei/ddecay/man/dist.decay.html
data(os)
dd <- dist.decay(gradient=os$dist, counts=os[, -1:-7], coords=os[, c("longitude", "latitude")], nboots=1000, dis.fun = "vegdist", method = "bray", like.pairs=T)
x <- vegdist(os$dist, method = "euclidean")
y <- 1-vegdist(os[, -1:-7], method = "bray")
plot(x, y)
lines(dd$Predictions[, "x"], dd$Predictions[,"mean"], col="red", lwd=2)



#install.packages("betapart")
library("betapart")

test.data <- xlsx.df.m
test.data$variable <- gsub("50.*", "55", test.data$variable)
test.data$variable <- gsub("60.*", "65", test.data$variable)
test.data$variable <- gsub("70.*", "75", test.data$variable)
test.data$variable <- gsub("80.*", "85", test.data$variable)
test.data$variable <- gsub("90.*", "95", test.data$variable)
test.data$variable <- as.numeric(test.data$variable)
head(test.data)
str(test.data)


plot(test.data$value, test.data$variable )


decay.model(test.data$value, test.data$variable, model.type="exponential", y.type="similarities", perm=100)



require(vegan)

data(BCI)
## UTM Coordinates (in metres)
UTM.EW <- rep(seq(625754, 626654, by=100), each=5)
UTM.NS <- rep(seq(1011569,  1011969, by=100), len=50)

spat.dist<-dist(data.frame(UTM.EW, UTM.NS))

dissim.BCI<-beta.pair.abund(BCI)$beta.bray.bal

plot(spat.dist, dissim.BCI, ylim=c(0,1), xlim=c(0, max(spat.dist)))

BCI.decay.exp<-decay.model(dissim.BCI, spat.dist, y.type="dissim", model.type="exp", perm=100)

BCI.decay.pow<-decay.model(dissim.BCI, spat.dist, y.type="dissim", model.type="pow", perm=100)

plot.decay(BCI.decay.exp, col=rgb(0,0,0,0.5))
plot.decay(BCI.decay.exp, col="red", remove.dots=TRUE, add=TRUE)
plot.decay(BCI.decay.pow, col="blue", remove.dots=TRUE, add=TRUE)





library(MASS)
set.seed(1) 
testData <- rnorm(1000) 


fitData <- function(data, fit="gamma", sample=0.5){
  distrib = list()
  numfit <- length(fit)
  results = matrix(0, ncol=5, nrow=numfit)
  
  for(i in 1:numfit){
    if((fit[i] == "gamma") | 
       (fit[i] == "poisson") | 
       (fit[i] == "weibull") | 
       (fit[i] == "exponential") |
       (fit[i] == "logistic") |
       (fit[i] == "normal") | 
       (fit[i] == "geometric")
    ) 
      distrib[[i]] = fit[i]
    else stop("Provide a valid distribution to fit data" )
  }
  
  # take a sample of dataset
  n = round(length(data)*sample)
  data = sample(data, size=n, replace=F)
  
  for(i in 1:numfit) {
    if(distrib[[i]] == "gamma") {
      gf_shape = "gamma"
      fd_g <- fitdistr(data, "gamma")
      est_shape = fd_g$estimate[[1]]
      est_rate = fd_g$estimate[[2]]
      
      ks = ks.test(data, "pgamma", shape=est_shape, rate=est_rate)
      
      # add to results
      results[i,] = c(gf_shape, est_shape, est_rate, ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "poisson"){
      gf_shape = "poisson"
      fd_p <- fitdistr(data, "poisson")
      est_lambda = fd_p$estimate[[1]]
      
      ks = ks.test(data, "ppois", lambda=est_lambda)
      # add to results
      results[i,] = c(gf_shape, est_lambda, "NA", ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "weibull"){
      gf_shape = "weibull"
      fd_w <- fitdistr(data,densfun=dweibull,start=list(scale=1,shape=2))
      est_shape = fd_w$estimate[[1]]
      est_scale = fd_w$estimate[[2]]
      
      ks = ks.test(data, "pweibull", shape=est_shape, scale=est_scale)
      # add to results
      results[i,] = c(gf_shape, est_shape, est_scale, ks$statistic, ks$p.value) 
    }
    
    else if(distrib[[i]] == "normal"){
      gf_shape = "normal"
      fd_n <- fitdistr(data, "normal")
      est_mean = fd_n$estimate[[1]]
      est_sd = fd_n$estimate[[2]]
      
      ks = ks.test(data, "pnorm", mean=est_mean, sd=est_sd)
      # add to results
      results[i,] = c(gf_shape, est_mean, est_sd, ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "exponential"){
      gf_shape = "exponential"
      fd_e <- fitdistr(data, "exponential")
      est_rate = fd_e$estimate[[1]]
      ks = ks.test(data, "pexp", rate=est_rate)
      # add to results
      results[i,] = c(gf_shape, est_rate, "NA", ks$statistic, ks$p.value)
    }
    
    else if(distrib[[i]] == "logistic"){
      gf_shape = "logistic"
      fd_l <- fitdistr(data, "logistic")
      est_location = fd_l$estimate[[1]]
      est_scale = fd_l$estimate[[2]]
      ks = ks.test(data, "plogis", location=est_location, scale=est_scale)
      # add to results
      results[i,] = c(gf_shape, est_location, est_scale, ks$statistic,    ks$p.value) 
    }
  }
  results = rbind(c("distribution", "param1", "param2", "ks stat", "ks    pvalue"),   results)
  #print(results)
  return(results)
}


res = fitData(testData, fit=c("logistic","normal","poisson"),
              sample=1)
res



library(fitdistrplus) 
require(MASS) 
set.seed(1) 
testData <- rnorm(1000) 
fitdist(testData, "gamma", method = "mme", start = list(shape = 0.1, rate = 0.1))



x <- c(15.771062,14.741310,9.081269,11.276436,11.534672,17.980860,13.550017,13.853336,11.262280,11.049087,14.752701,4.481159,11.680758,11.451909,10.001488,11.106817,7.999088,10.591574,8.141551,12.401899,11.215275,13.358770,8.388508,11.875838,3.137448,8.675275,17.381322,12.362328,10.987731,7.600881,14.360674,5.443649,16.024247,11.247233,9.549301,9.709091,13.642511,10.892652,11.760685,11.717966,11.373979,10.543105,10.230631,9.918293,10.565087,8.891209,10.021141,9.152660,10.384917,8.739189,5.554605,8.575793,12.016232,10.862214,4.938752,14.046626,5.279255,11.907347,8.621476,7.933702,10.799049,8.567466,9.914821,7.483575,11.098477,8.033768,10.954300,8.031797,14.288100,9.813787,5.883826,7.829455,9.462013,9.176897,10.153627,4.922607,6.818439,9.480758,8.166601,12.017158,13.279630,14.464876,13.319124,12.331335,3.194438,9.866487,11.337083,8.958164,8.241395,4.289313,5.508243,4.737891,7.577698,9.626720,16.558392,10.309173,11.740863,8.761573,7.099866,10.032640)
qqnorm(x); qqline(x)

fitdistr(x, "lognormal")
fitdistr(x, "beta", start=) 
fitdistr(x, "cauchy")
fitdistr(x, "chi-squared", start=)
fitdistr(x, "exponential")
fitdistr(x, "gamma")
fitdistr(x, "geometric")
fitdistr(x, "log-normal")
fitdistr(x, "lognormal")
fitdistr(x, "logistic")
fitdistr(x, "negative binomial")
fitdistr(x, "normal" )
fitdistr(x, "Poisson")
fitdistr(x, "t")
fitdistr(x, "weibull")



#https://www.statmethods.net/advgraphs/probability.html
qqplot(qchisq(ppoints(500), df = 3), x); qqline(x)
qqplot(qexp(ppoints(500)), x); qqline(x)
qqplot(qbinom(ppoints(500), size=10, prob=.5), x); qqline(x)
qqplot(qpois(ppoints(500), lambda=1), x); qqline(x)
qqplot(qnbinom(ppoints(500), size=10, prob=.5), x); qqline(x)
qqplot(qlnorm(ppoints(500)), x); qqline(x)
qqplot(qgamma(ppoints(500), shape=1), x); qqline(x)
qqplot(qgeom(ppoints(500), prob=0.5), x); qqline(x)
qqplot(qbeta(ppoints(500), shape1=1, shape2=1), x); qqline(x)
qqplot(qhyper(ppoints(500), m=2, n=2, k=2), x); qqline(x)
qqplot(qcauchy(ppoints(500)), x); qqline(x)
qqplot(qf(ppoints(500), df1=3, df2=3), x); qqline(x)
qqplot(qlogis(ppoints(500)), x); qqline(x)





library("fitdistrplus")
data("groundbeef")
str(groundbeef)
head(groundbeef)


plotdist(groundbeef$serving, histo = TRUE, demp = TRUE)
descdist(groundbeef$serving, boot=1000)


fw <- fitdist(groundbeef$serving, "weibull")
summary(fw)

fg <- fitdist(groundbeef$serving,"gamma")
summary(fg)
fln <- fitdist(groundbeef$serving,"lnorm")
summary(fln)

par(mfrow=c(2, 2))
denscomp(list(fw,fln,fg), legendtext=c("Weibull", "lognormal", "gamma"))
qqcomp(list(fw,fln,fg), legendtext=c("Weibull", "lognormal", "gamma"))
cdfcomp(list(fw,fln,fg), legendtext=c("Weibull", "lognormal", "gamma"))
ppcomp(list(fw,fln,fg), legendtext=c("Weibull", "lognormal", "gamma"))


message("+-------------------------------------------------------------------------------")
message("+ Read in Tables S2 ")
message("+-------------------------------------------------------------------------------")

pairs = read.table(paste0(baseDir, "Tables_S2.pairs.txt"), sep="\t", header=F)
head(pairs)
str(pairs)

print(pairs)


ggplot(data=subset(pairs, abs(V7) < 2500000 & V2 > 0), aes(x=(abs(V7)), y=abs(V2-V5), color=paste0(V3, "_", V6))) +
  geom_smooth(method = "lm", alpha = .1, aes(fill = paste0(V3, "_", V6))) +
  geom_point(alpha=.5) +
  
  xlab("Distance") + 
  ylab("Bias Difference") +
  # xlim(-5,20) +
  facet_wrap( ~ paste0(V3, "_", V6), nrow=3, scales="free_x") +
  theme_bw() +
  theme(legend.position="none")



df <- subset(pairs, abs(V7) < 2500000 & V2 > 0)
ggplot(data=df, aes(x=abs(V7), y=abs(V2-V5), color=paste0(V3, "_", V6))) +
  geom_jitter(alpha=.5, size=4) +
  
  ylab("Between gene pair Bias Difference") +
  xlab("Between gene pair Distance") + 
  
  ggtitle("Do nearby Genes have same bias?") +
  theme_bw() +
  theme(legend.position="right")



df <- subset(pairs, abs(V7) < 250000 & V2 > 0)
ggplot(data=df, aes(x=abs(V7), y=abs(V2-V5), fill=as.factor(abs(V2-V5)))) +
  geom_density_ridges(scale = 0.9, jittered_points = TRUE) + 
  ylab("Between gene pair Bias Difference") +
  xlab("Between gene pair Distance") + 
  ggtitle("Do nearby genes have same bias?") +
  theme_bw() +
  theme(legend.position="none")



dlk1.dists = read.table(paste0(baseDir, "test_ICR/Distances_Dlk1.bed"), sep="\t", header=F)
head(dlk1.dists)


ggplot(data=dlk1.dists, aes(x=V12, y=V8)) +
  geom_vline(xintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=subset(dlk1.dists, V9=="+"), aes(x=(V12), y=V8, shape=as.factor(V11)), alpha=.75, size=4, colour="blue") +
  geom_point(data=subset(dlk1.dists, V9=="-"), aes(x=(V12), y=V8, shape=as.factor(V11)), alpha=.75, size=4, colour="pink") +
  geom_smooth(se=F, colour="black") +
  xlab("Distance from ICR") +
  ylab("Bias (mean)") +
  ylim(50,100) +
  ggtitle("Dlk1") +
  scale_color_manual(name="Cell Type", 
                     values=c("+"="blue", "-"="pink") ) +
  theme_bw() 





ggplot(data=dlk1.dists, aes(x=abs(V12), y=V8)) +
  geom_vline(xintercept = 0, colour="red", linetype="dashed") +
  geom_point(data=subset(dlk1.dists, V9=="+"), aes(x=abs(V12), y=V8, shape=as.factor(V11)), alpha=.75, size=4, colour="blue") +
  geom_point(data=subset(dlk1.dists, V9=="-"), aes(x=abs(V12), y=V8, shape=as.factor(V11)), alpha=.75, size=4, colour="pink") +
  geom_smooth(se=F, colour="black") +
  xlab("Distance from ICR") +
  ylab("Bias (mean)") +
  ggtitle("Dlk1") +
  scale_color_manual(name="Cell Type", 
                     values=c("+"="blue", "-"="pink") ) +
  theme_bw() 



head(dlk1.dists)

ggplot(dlk1.dists, aes(x = abs(V12), y = as.factor(V8), fill=as.factor(V8))) +
  geom_density_ridges(scale = 0.9, jittered_points = TRUE) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + 
  xlab("Distance From ICR") +
  ylab("Bias") +
  ggtitle("Dlk1") +
  theme_ridges() +
  theme(legend.position="none")


ggplot(dlk1.dists, aes(x = abs(V12), y = as.factor(V8), fill=as.factor(V8)) ) +
  geom_density_ridges( jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7,
                       vline_size = 1, vline_color = "red",
                       point_size = 0.99, point_alpha = 1, position = position_raincloud(adjust_vlines=T) ) +
  scale_y_discrete(expand = c(0, 0)) +     
  scale_x_continuous(expand = c(0, 0), breaks=c(0, 500000, 1000000, 1500000), labels=c("0", "0.5Mb", "1Mb", "1.5Mb")) +  
  coord_cartesian(clip = "off") + 
  xlab("Distance From ICR") +
  ylab("Bias") +
  ggtitle("Dlk1") +
  theme_ridges() +
  theme(legend.position="none")


ggplot(dlk1.dists, aes(x = abs(V12), y = as.character(V8), alpha=V8), colour='black', fill='black') +
  geom_density_ridges( jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, #alpha = 0.7,
                       vline_size = 1, vline_color = "red",
                       point_size = 0.99, point_alpha = 1, position = position_raincloud(adjust_vlines=T) ) +
  scale_alpha(range = c(0.25, 1.0)) +
  scale_y_discrete(expand = c(0, 0)) +     
  scale_x_continuous(expand = c(0, 0), breaks=c(0, 500000, 1000000, 1500000), labels=c("0", "0.5Mb", "1Mb", "1.5Mb")) +  
  coord_cartesian(clip = "off") + 
  xlab("Distance From ICR") +
  ylab("Bias") +
  ggtitle("Dlk1") +
  theme_ridges() +
  theme(legend.position="none")


message("+-------------------------------------------------------------------------------")
message("+ Read in ISOLDE Results ")
message("+-------------------------------------------------------------------------------")

dlk1.dists.test           <- dlk1.dists

isolde.ARN.ASE        <- read.table( paste0(baseDir, "Isolde/ARN/", "ISoLDE_result_ASE_04-07-2020_18-27-13.tsv"), header=T, sep="\t")#, row.names=1)
isolde.ARN.BA         <- read.table( paste0(baseDir, "Isolde/ARN/", "ISoLDE_result_BA_04-07-2020_18-27-13.tsv"),  header=T, sep="\t")#, row.names=1)
isolde.ARN.UN         <- read.table( paste0(baseDir, "Isolde/ARN/", "ISoLDE_result_UN_04-07-2020_18-27-13.tsv"),  header=T, sep="\t")#, row.names=1)
isolde.ARN.ASE$flag   <- "none"
isolde.ARN.BA$origin  <- "none"
isolde.ARN.BA$flag    <- "none"
isolde.ARN.all        <- rbind(isolde.ARN.ASE, isolde.ARN.BA, isolde.ARN.UN)
isolde.ARN.all$tissue <- "ARN"

isolde.DRN.ASE        <- read.table( paste0(baseDir, "Isolde/DRN/", "ISoLDE_result_ASE_04-20-2020_16-15-29.tsv"), header=T, sep="\t")#, row.names=1)
isolde.DRN.BA         <- read.table( paste0(baseDir, "Isolde/DRN/", "ISoLDE_result_BA_04-20-2020_16-15-29.tsv"), header=T, sep="\t")#, row.names=1)
isolde.DRN.UN         <- read.table( paste0(baseDir, "Isolde/DRN/", "ISoLDE_result_UN_04-20-2020_16-15-29.tsv"), header=T, sep="\t")#, row.names=1)
isolde.DRN.ASE$flag   <- "none"
isolde.DRN.BA$origin  <- "none"
isolde.DRN.BA$flag    <- "none"
isolde.DRN.all        <- rbind(isolde.DRN.ASE, isolde.DRN.BA, isolde.DRN.UN)
isolde.DRN.all$tissue <- "DRN"

isolde.LVR.ASE        <- read.table( paste0(baseDir, "Isolde/Liver/", "ISoLDE_result_ASE_04-08-2020_11-40-18.tsv"), header=T, sep="\t")#, row.names=1)
isolde.LVR.BA         <- read.table( paste0(baseDir, "Isolde/Liver/", "ISoLDE_result_BA_04-08-2020_11-40-18.tsv"), header=T, sep="\t")#, row.names=1)
isolde.LVR.UN         <- read.table( paste0(baseDir, "Isolde/Liver/", "ISoLDE_result_UN_04-08-2020_11-40-18.tsv"), header=T, sep="\t")#, row.names=1)
isolde.LVR.ASE$flag   <- "none"
isolde.LVR.BA$origin  <- "none"
isolde.LVR.BA$flag    <- "none"
isolde.LVR.all        <- rbind(isolde.LVR.ASE, isolde.LVR.BA, isolde.LVR.UN)
isolde.LVR.all$tissue <- "Liver"

subset(isolde.LVR.all, name=="Meg3")


isolde.MSC.ASE        <- read.table( paste0(baseDir, "Isolde/Muscle/", "ISoLDE_result_ASE_04-08-2020_15-28-14.tsv"), header=T, sep="\t")#, row.names=1)
isolde.MSC.BA         <- read.table( paste0(baseDir, "Isolde/Muscle/", "ISoLDE_result_BA_04-08-2020_15-28-14.tsv"), header=T, sep="\t")#, row.names=1)
isolde.MSC.UN         <- read.table( paste0(baseDir, "Isolde/Muscle/", "ISoLDE_result_UN_04-08-2020_15-28-14.tsv"), header=T, sep="\t")#, row.names=1)
isolde.MSC.ASE$flag   <- "none"
isolde.MSC.BA$origin  <- "none"
isolde.MSC.BA$flag    <- "none"
isolde.MSC.all        <- rbind(isolde.MSC.ASE, isolde.MSC.BA, isolde.MSC.UN)
isolde.MSC.all$tissue <- "Muscle"

isolde.P8C.ASE        <- read.table( paste0(baseDir, "Isolde/P8_Cerebellum_merge/", "ISoLDE_result_ASE_06-23-2020_16-21-15.tsv"), header=T, sep="\t")#, row.names=1)
isolde.P8C.BA         <- read.table( paste0(baseDir, "Isolde/P8_Cerebellum_merge/", "ISoLDE_result_BA_06-23-2020_16-21-15.tsv"), header=T, sep="\t")#, row.names=1)
isolde.P8C.UN         <- read.table( paste0(baseDir, "Isolde/P8_Cerebellum_merge/", "ISoLDE_result_UN_06-23-2020_16-21-15.tsv"), header=T, sep="\t")#, row.names=1)
isolde.P8C.ASE$flag   <- "none"
isolde.P8C.BA$origin  <- "none"
isolde.P8C.BA$flag    <- "none"
isolde.P8C.all        <- rbind(isolde.P8C.ASE, isolde.P8C.BA, isolde.P8C.UN)
isolde.P8C.all$tissue <- "P8_Cerebellum"


isolde.P60C.ASE        <- read.table( paste0(baseDir, "Isolde/P60_Cerebellum_merge/", "ISoLDE_result_ASE_06-23-2020_16-27-10.tsv"), header=T, sep="\t")#, row.names=1)
isolde.P60C.BA         <- read.table( paste0(baseDir, "Isolde/P60_Cerebellum_merge/", "ISoLDE_result_BA_06-23-2020_16-27-10.tsv"), header=T, sep="\t")#, row.names=1)
isolde.P60C.UN         <- read.table( paste0(baseDir, "Isolde/P60_Cerebellum_merge/", "ISoLDE_result_UN_06-23-2020_16-27-10.tsv"), header=T, sep="\t")#, row.names=1)
isolde.P60C.ASE$flag   <- "none"
isolde.P60C.BA$origin  <- "none"
isolde.P60C.BA$flag    <- "none"
isolde.P60C.all        <- rbind(isolde.P60C.ASE, isolde.P60C.BA, isolde.P60C.UN)
isolde.P60C.all$tissue <- "P60_Cerebellum"


isolde.all     <- rbind(isolde.ARN.all, isolde.DRN.all, isolde.LVR.all, isolde.MSC.all, isolde.P8C.all, isolde.P60C.all)
isolde.all.ann <- merge(dlk1.dists.test, isolde.all, by.x="V7", by.y="name")

isolde.all.ann$tissue <- factor(isolde.all.ann$tissue,      # Reordering group factor levels
                         levels = c("Liver", "Muscle", "ARN", "DRN", "P8_Cerebellum", "P60_Cerebellum"))

pdf(paste0(baseDir, "ICR_Dist_vs_Bias_tissue.pdf"), width=10,height=12.5, onefile=FALSE)
par(bg=NA)
ggplot(isolde.all.ann, aes(x = abs(V12), y = diff_prop, alpha=variability, size=criterion, label=V7, colour=origin)) +
  geom_smooth(inherit.aes=F, data=subset(isolde.all.ann, origin=="P" | origin=="none")[, c("V12", "diff_prop")], aes(x = abs(V12), y = diff_prop), alpha=0.2, colour="blue", se=T ) +
  geom_smooth(inherit.aes=F, data=subset(isolde.all.ann, origin=="M" | origin=="none")[, c("V12", "diff_prop")], aes(x = abs(V12), y = diff_prop), alpha=0.2, colour="pink", se=T ) +
  
  geom_label_repel(size=3, alpha=1,  max.overlaps=25) +
  geom_point() +
  scale_color_manual(name="Bias", breaks=c("P", "M", "none"), values= c("P"="blue", "M"="pink", "none"="black")) +
  xlab("Distance From ICR") +
  ylab("Parental Bias") +
  ggtitle("ICR") +
  facet_wrap(~ tissue, ncol=2) +
  theme_bw() +
  theme(legend.position="bottom")
dev.off()






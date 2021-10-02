
library("ggplot2")


baseDir <- "/Users/rhamilto/Desktop/2021-Carol_Paper/"


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

  
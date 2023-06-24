################################################################################
#
#   Article : Multiclass classification of growth curves using random change 
#             points and heterogeneous random effects
#
#   Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
#             Scott A. Sisson
#
#   Required R packages : readr, mcclust, ggplot2, dplyr, gridExtra
#
################################################################################

# Required packages
library(readr)
library(mcclust)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Find optimal clustering
est_grp_label = read_csv("../results/application/est_grp_label.csv", 
                          col_names = FALSE)
psm     = comp.psm(t(est_grp_label))
result  = maxpear(psm)
cluster = result$cl
cluster = cluster + max(cluster)

# Ordering cluster and saving output
cluster_order = order(table(cluster), decreasing=T)
for(i in 1:length(cluster_order)){
  cluster[cluster==cluster_order[i]+length(cluster_order)] = i
}
write.csv(cluster, file="../results/application/cluster.csv", row.names=F)

# Loading attributes of each individuals
load("../data/attributes.Rdata")
dataset = cbind(attributes, cluster)

# Plotting figure 6
p1=ggplot(rbind(dataset, cbind(dataset[,1:6], data.frame(cluster="Pop"))), 
          aes(x=factor(cluster), group=SEX, fill=SEX)) + 
  geom_bar(position = "fill") +
  xlab('Subgroup') +
  ylab('Proportion') +
  scale_fill_manual("Gender", values=c("Female"="#D6D6D6", "Male"="#000000")) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18))

p2=ggplot(rbind(dataset, cbind(dataset[,1:6], data.frame(cluster="Pop"))), 
          aes(x=factor(cluster), group=MEDUCYRS, fill=MEDUCYRS)) + 
  geom_bar(position = "fill") +
  xlab('Subgroup') +
  ylab('Proportion') +
  scale_fill_manual("Mother's\nEducation\n(Years)", values=c("0"="#D6D6D6", 
                    "5"="#AAAAAA", "8 or more"="#000000")) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18))
 
dataset = dataset %>% mutate(color=if_else(cluster<7, "transparent", "black"))

p4=ggplot(dataset, aes(x=as.factor(cluster), y=PIQ)) + 
  geom_boxplot(varwidth=TRUE, outlier.shape=NA, col=c(rep("black",6), 
               rep("transparent",2))) +
  geom_point(col=dataset$color) +
  xlab('Subgroup') +
  ylab('Performance IQ') +  
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))

p5=ggplot(dataset, aes(x=as.factor(cluster), y=VIQ)) + 
  geom_boxplot(varwidth=TRUE, outlier.shape=NA, col=c(rep("black",6), 
               rep("transparent",2))) +
  geom_point(col=dataset$color) +
  xlab('Subgroup') +
  ylab('Verbal IQ') +  
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18))

p3=ggplot(dataset, aes(x=as.factor(cluster), y=GENIQ)) + 
  geom_boxplot(varwidth=TRUE, outlier.shape=NA, col=c(rep("black",6), 
               rep("transparent",2))) +
  geom_point(col=dataset$color) +
  xlab('Subgroup') +
  ylab('General Intelligence IQ') +  
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=18)) 

ggsave(file="../figures/figure6.eps", plot=grid.arrange(p1,p3,p5,p2,p4, ncol=3),
       width=18, height=9)
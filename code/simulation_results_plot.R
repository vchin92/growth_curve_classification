################################################################################
#
#   Article : Multiclass classification of growth curves using random change 
#             points and heterogeneous random effects
#
#   Authors : Vincent Chin, Jarod Lee, Louise M. Ryan, Robert Kohn, and
#             Scott A. Sisson
#
#   Required R packages : mcclust, ggplot2, ggh4x, dplyr, data.table, gridExtra
#
################################################################################

# Required packages
library(mcclust)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(data.table)
library(gridExtra)

# Function for generating all permutations of elements in a vector
perm = function(v){
  n = length(v)
  if(n==1){
    X = v
  }else{
    X = NULL
    for(i in 1:n) X = rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

# Function for performance evaluation in terms of the proportion of correct
# classification for each truth label and the maximum relative size of 
# additional subgroups
perf_measure = function(x1, x2, permute){
  # x1 is the truth label
  # x2 is the classification label
  
  group_size = data.frame(x2=x2) %>% 
    group_by(x2) %>%
    tally()
  
  if(nrow(group_size)>=4){
    permutation = permute[[length(unique(x2))]]
    permutation = as.vector(t(permutation))
    
    accuracy = sapply(1:(length(permutation)/4), function(j){
      cont_tab = table(x1, x2)[,permutation[((j-1)*4+1):(j*4)]]
      c(diag(cont_tab)/tabulate(x1), sum(diag(cont_tab))/length(x1))
    })
    
    perm_max_acc     = permutation[((which.max(accuracy[5,])-1)*4+1):
                                     (which.max(accuracy[5,])*4)]
    class_acc        = accuracy[1:4,which.max(accuracy[5,])]
    names(class_acc) = as.character(1:4)
    
    add_group = group_size %>%
      filter(!x2 %in% perm_max_acc)
    
    if(nrow(add_group)==0){
      add_group = 0
    }else{
      add_group = max(add_group$n)/length(x2)
    }
  }else{
    cont_tab       = table(x1, x2)
    cont_tab       = cont_tab[apply(cont_tab, 2, which.max),]
    max_acc        = diag(cont_tab)
    names(max_acc) = row.names(cont_tab)
    class_acc      = c(max_acc, rep(0, 4-length(max_acc)))
    
    names(class_acc)[names(class_acc)==""] = as.character(setdiff(1:4, 
                                             as.numeric(names(max_acc))))
    
    class_acc = class_acc[order(as.numeric(names(class_acc)))]/tabulate(x1)
    add_group = 0
  }
  return(list(cbind(as.data.frame(t(class_acc)), N=length(x2), 
                    rel_size=add_group)))
}

# A list to store permutations
permute = vector("list", 12)
for(i in 1:12){
  if(i < 4){
    permute[[i]] = NULL
  }else{
    combination  = t(combn(1:i, 4))
    permute[[i]] = do.call("rbind", lapply(1:nrow(combination), 
                                           function(j) perm(combination[j,])))
  }
}

################################################################################
# Data for truth label
load("../data/label.Rdata")

################################################################################
# Results for D_fixed and M_fixed
lapply(1:8, function(x) load(paste0("../results/simulation/df_mf_",x,".Rdata"),
                             .GlobalEnv))
df_mf = c(df_mf_1, df_mf_2, df_mf_3, df_mf_4, df_mf_5, df_mf_6, df_mf_7, df_mf_8)
rm(df_mf_1, df_mf_2, df_mf_3, df_mf_4, df_mf_5, df_mf_6, df_mf_7, df_mf_8)
df_mf_result = vector("list", length=8)

for(i in 1:8){
  N = (i-1)*50+50
  df_mf_result[[i]] = sapply(1:100, function(j){
    psm = comp.psm(t(df_mf[[i]][,,j]))
    pear = maxpear(psm)
    list(list(as.data.frame(pear$cl)), 
         list(data.frame(Data="D[fixed]", Model="M[fixed]", N=N, 
                         ARI=arandi(label[[i]][,j], pear$cl), G=max(pear$cl))))
  })
}

df_mf_cl = lapply(df_mf_result, 
                  function(j) as.data.frame(bind_cols(j[seq(1,200,2)])))
df_mf_result = do.call("rbind", lapply(df_mf_result, 
                                       function(j) bind_rows(j[seq(2,200,2)])))

rm(df_mf)

################################################################################
# Results for D_fixed and M_random
lapply(1:8, function(x) load(paste0("../results/simulation/df_mr_",x,".Rdata"),
                             .GlobalEnv))
df_mr = c(df_mr_1, df_mr_2, df_mr_3, df_mr_4, df_mr_5, df_mr_6, df_mr_7, df_mr_8)
rm(df_mr_1, df_mr_2, df_mr_3, df_mr_4, df_mr_5, df_mr_6, df_mr_7, df_mr_8)
df_mr_result = vector("list", length=8)

for(i in 1:8){
  N = (i-1)*50+50
  df_mr_result[[i]] = sapply(1:100, function(j){
    psm = comp.psm(t(df_mr[[i]][,,j]))
    pear = maxpear(psm)
    list(list(as.data.frame(pear$cl)), 
         list(data.frame(Data="D[fixed]", Model="M[random]", N=N, 
                         ARI=arandi(label[[i]][,j], pear$cl), G=max(pear$cl))))
  })
}

df_mr_cl = lapply(df_mr_result, 
                  function(j) as.data.frame(bind_cols(j[seq(1,200,2)])))
df_mr_result = do.call("rbind", lapply(df_mr_result, 
                                       function(j) bind_rows(j[seq(2,200,2)])))

rm(df_mr)

################################################################################
# Results for D_random and M_fixed
lapply(1:8, function(x) load(paste0("../results/simulation/dr_mf_",x,".Rdata"),
                             .GlobalEnv))
dr_mf = c(dr_mf_1, dr_mf_2, dr_mf_3, dr_mf_4, dr_mf_5, dr_mf_6, dr_mf_7, dr_mf_8)
rm(dr_mf_1, dr_mf_2, dr_mf_3, dr_mf_4, dr_mf_5, dr_mf_6, dr_mf_7, dr_mf_8)
dr_mf_result = vector("list", length=8)

for(i in 1:8){
  N = (i-1)*50+50
  dr_mf_result[[i]] = sapply(1:100, function(j){
    psm = comp.psm(t(dr_mf[[i]][,,j]))
    pear = maxpear(psm)
    list(list(as.data.frame(pear$cl)),
         list(data.frame(Data="D[random]", Model="M[fixed]", N=N, 
                         ARI=arandi(label[[i]][,j], pear$cl), G=max(pear$cl))))
  })
}

dr_mf_cl = lapply(dr_mf_result, 
                  function(j) as.data.frame(bind_cols(j[seq(1,200,2)])))
dr_mf_result = do.call("rbind", lapply(dr_mf_result, 
                                       function(j) bind_rows(j[seq(2,200,2)])))

rm(dr_mf)

################################################################################
# Results for D_random and M_random
lapply(1:8, function(x) load(paste0("../results/simulation/dr_mr_",x,".Rdata"),
                             .GlobalEnv))
dr_mr = c(dr_mr_1, dr_mr_2, dr_mr_3, dr_mr_4, dr_mr_5, dr_mr_6, dr_mr_7, dr_mr_8)
rm(dr_mr_1, dr_mr_2, dr_mr_3, dr_mr_4, dr_mr_5, dr_mr_6, dr_mr_7, dr_mr_8)
dr_mr_result = vector("list", length=8)

for(i in 1:8){
  N = (i-1)*50+50
  dr_mr_result[[i]] = sapply(1:100, function(j){
    psm = comp.psm(t(dr_mr[[i]][,,j]))
    pear = maxpear(psm)
    list(list(as.data.frame(pear$cl)),
         list(data.frame(Data="D[random]", Model="M[random]", N=N, 
                         ARI=arandi(label[[i]][,j], pear$cl), G=max(pear$cl))))
  })
}

dr_mr_cl = lapply(dr_mr_result, 
                  function(j) as.data.frame(bind_cols(j[seq(1,200,2)])))
dr_mr_result = do.call("rbind", lapply(dr_mr_result, 
                                       function(j) bind_rows(j[seq(2,200,2)])))

rm(dr_mr)

################################################################################
# ARI and G
result = rbind(df_mf_result, df_mr_result, dr_mf_result, dr_mr_result)

# Figure 2a
p1 = ggplot(result, aes(x=N, y=ARI)) + 
  geom_boxplot(aes(group=N), outlier.shape=NA, lwd=0.2)+
  facet_nested_wrap(~Data+Model, nrow=1, labeller=label_parsed) +
  scale_x_continuous(breaks=seq(50,400,by=50)) +
  ylim(c(0.25,1))+
  theme_bw() +
  theme(strip.text=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.x=element_text(size=8, angle=90, vjust=0.5, hjust=1))

# Figure 2b
p2 = ggplot(result, aes(x=N, y=G)) + 
  geom_count() +
  scale_size_continuous(range=c(0.25,2.5)) + 
  facet_nested_wrap(~Data+Model, nrow=1, labeller=label_parsed) +
  scale_x_continuous(breaks=seq(50,400,by=50)) +
  scale_y_continuous(breaks=seq(2,12,by=2)) +
  theme_bw() +
  theme(strip.text=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.x=element_text(size=8, angle=90, vjust=0.5, hjust=1),
        legend.position="none")

# Accuracy for each label and relative size of extra G
df_mf_cl_perf = lapply(1:8, function(i){
  sapply(1:100, function(j) 
    perf_measure(label[[i]][,j], df_mf_cl[[i]][,j], permute))})
df_mr_cl_perf = lapply(1:8, function(i){
  sapply(1:100, function(j) 
    perf_measure(label[[i]][,j], df_mr_cl[[i]][,j], permute))})
dr_mf_cl_perf = lapply(1:8, function(i){
  sapply(1:100, function(j) 
    perf_measure(label[[i]][,j], dr_mf_cl[[i]][,j], permute))})
dr_mr_cl_perf = lapply(1:8, function(i){
  sapply(1:100, function(j) 
    perf_measure(label[[i]][,j], dr_mr_cl[[i]][,j], permute))})

df_mf_cl_perf = cbind(Data="D[fixed]", Model="M[fixed]", 
                      bind_rows(df_mf_cl_perf))
df_mr_cl_perf = cbind(Data="D[fixed]", Model="M[random]", 
                      bind_rows(df_mr_cl_perf))
dr_mf_cl_perf = cbind(Data="D[random]", Model="M[fixed]", 
                      bind_rows(dr_mf_cl_perf))
dr_mr_cl_perf = cbind(Data="D[random]", Model="M[random]", 
                      bind_rows(dr_mr_cl_perf))

result_perf = rbind(
  melt(setDT(df_mf_cl_perf), id.vars=c(1,2,7,8), variable.name="Truth"),
  melt(setDT(df_mr_cl_perf), id.vars=c(1,2,7,8), variable.name="Truth"),
  melt(setDT(dr_mf_cl_perf), id.vars=c(1,2,7,8), variable.name="Truth"),
  melt(setDT(dr_mr_cl_perf), id.vars=c(1,2,7,8), variable.name="Truth")
)
result_perf$N = as.factor(result_perf$N)
result_perf = result_perf %>% filter(value!=0)

# Figure 2c
p3 = ggplot(result_perf, aes(x=Truth, y=value, fill=N)) + 
  geom_boxplot(outlier.shape=NA, lwd=0.1) +
  ylab("Proportion of Correct Classification") +
  facet_nested_wrap(~Data+Model, nrow=1, labeller=label_parsed) +
  ylim(c(0.25,1))+
  theme_bw() +
  theme(strip.text=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.x=element_text(size=8, angle=90, vjust=0.5, hjust=1),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8))

# Figure 3
p4 = ggplot(dr_mr_cl_perf %>% filter(rel_size>0), aes(x=N, y=rel_size)) + 
  geom_boxplot(aes(group=N), outlier.shape=NA) +
  stat_summary(fun=mean, colour="red", geom="point", 
               shape=19, size=3, show.legend=FALSE) +
  ylab('Maximum Size of Additional Subgroups Relative to N') +
  scale_x_continuous(breaks=seq(50,400,by=50)) +
  theme_bw() +
  theme(strip.text=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x=element_text(size=12, angle=90, vjust=0.5, hjust=1),
        legend.position="none")

ggsave(file="../figures/figure2.eps", plot=grid.arrange(p1,p2,p3, nrow=3), width=18, 
       height=18, units="cm")
ggsave(file="../figures/figure3.eps", plot=p4, width=24, 
       height=12, units="cm")

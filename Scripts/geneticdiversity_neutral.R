##Plotting thetas calculated in ANGSD with neutral SNPs and transitions EXCLUDED

library(tidyverse)
library(cowplot)
library(boot)

##Load in theta outputs by population (abol, amta, cbol, cmta)

abol_thetas_neutral <- read_table("abol_neutral_notrans.thetas.idx.pestPG")
amta_thetas_neutral <- read_table("amta_neutral_notrans.thetas.idx.pestPG")
cbol_thetas_neutral <- read_table("cbol_neutral_notrans.thetas.idx.pestPG")
cmta_thetas_neutral <- read_table("cmta_neutral_notrans.thetas.idx.pestPG")

#Add a column for site to each dataset

abol_thetas_neutral$Site <- "ABol"
amta_thetas_neutral$Site <- "AMta"
cbol_thetas_neutral$Site <- "CBol"
cmta_thetas_neutral$Site <- "CMta"

#Add a column for era to each dataset
abol_thetas_neutral$Era <- "Historical"
amta_thetas_neutral$Era <- "Historical"
cbol_thetas_neutral$Era <- "Contemporary"
cmta_thetas_neutral$Era <- "Contemporary"

#Merge datasets for plotting

angsd_thetas_neutral <- rbind(abol_thetas_neutral, cbol_thetas_neutral, amta_thetas_neutral, cmta_thetas_neutral)

##Import depth information

angsd_depth_neutral <- read_table("all_neutral.pos.gz")
#This is by site, we need it by chromosome to match with the theta calculations

#Take the average of total Depth for all the sites within a chromosome
angsd_avgdepth_neutral <- angsd_depth_neutral %>% group_by(chr) %>% summarise(avg_Depth = mean(totDepth))
#This gets us average total depth (summed across individuals) for each chromosome

#Add a column with average depth by individual
angsd_avgdepth_neutral$ind_depth <- angsd_avgdepth_neutral$avg_Depth/222

#Rename Chr column to match angsd_thetas dataframe
names(angsd_avgdepth_neutral)[1] <- "Chr"

#Merge with angsd_thetas dataframe

angsd_thetas_depth_neutral <- merge(angsd_thetas_neutral, angsd_avgdepth_neutral, by="Chr")

#Add a column with theta by site (theta for the contig divided by number of sites on the contig)
angsd_thetas_depth_neutral$tW_bysite <- angsd_thetas_depth_neutral$tW/angsd_thetas_depth_neutral$nSites

#Reorder for plotting
angsd_thetas_depth_neutral$Site <- factor(angsd_thetas_depth_neutral$Site, levels=c("ABol", "CBol", "AMta", "CMta"))

#Theta by Depth plot

angsd_thetas_depth_neutral %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

angsd_thetas_depth_neutral_bol <- subset(angsd_thetas_depth_neutral, Site=="ABol" | Site=="CBol")
angsd_thetas_depth_neutral_mta <- subset(angsd_thetas_depth_neutral, Site=="AMta" | Site=="CMta")

bol_plot <- angsd_thetas_depth_neutral_bol %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era), show.legend=FALSE) +
  labs(x="Mean Depth Per Individual",y="Theta") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(breaks=c(0.04, 0.06, 0.08, 0.10, 0.12)) +
  ylim(0.025, 0.125)

mta_plot <- angsd_thetas_depth_neutral_mta %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Mean Depth Per Individual",y="Theta") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(breaks=c(0.04, 0.06, 0.08, 0.10, 0.12)) +
  ylim(0.025, 0.125)

plot_grid(bol_plot, mta_plot,
          labels=c('A', 'B'))

plot_grid(bol_plot_all, mta_plot_all, bol_plot, mta_plot,
          labels=c('A', 'B', 'C', 'D'))
#bol_plot_all and mta_plot_all from the geneticdiversity.R script

##Pairwise theta (tP- also referred to as thetaD)

#Add a column with tP by site (theta for the contig divided by number of sites on the contig)
angsd_thetas_depth_neutral$tP_bysite <- angsd_thetas_depth_neutral$tP/angsd_thetas_depth_neutral$nSites

#Theta by Depth plot

angsd_thetas_depth_neutral %>%
  ggplot(aes(x=ind_depth, y=tP_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

angsd_thetas_depth_neutral_bol <- subset(angsd_thetas_depth_neutral, Site=="ABol" | Site=="CBol")
angsd_thetas_depth_neutral_mta <- subset(angsd_thetas_depth_neutral, Site=="AMta" | Site=="CMta")

bol_tP <- angsd_thetas_depth_neutral_bol %>%
  ggplot(aes(x=ind_depth, y=tP_bysite, color=Era), show.legend=FALSE) +
  labs(x="Mean Depth Per Individual",y="Nucleotide Diversity") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  ylim(0.02, 0.18)

mta_tP <- angsd_thetas_depth_neutral_mta %>%
  ggplot(aes(x=ind_depth, y=tP_bysite, color=Era)) +
  labs(x="Mean Depth Per Individual",y="Nucleotide Diversity") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  ylim(0.02, 0.18)


plot_grid(bol_tP, mta_tP,
          labels=c('A', 'B'))

plot_grid(bol_tp_all, mta_tp_all, bol_tP, mta_tP,
          labels=c('A', 'B', 'C', 'D'))
#bol_tp_all and mta_tp_all from the geneticdiversity.R script

#Subset to depth 3-6X

angsd_thetas_depth_neutral_36 <- subset(angsd_thetas_depth_neutral, ind_depth >= 3 & ind_depth<=6)

angsd_thetas_depth_neutral_36_bol <- subset(angsd_thetas_depth_neutral_36, Site=="ABol" | Site=="CBol")
angsd_thetas_depth_neutral_36_mta <- subset(angsd_thetas_depth_neutral_36, Site=="AMta" | Site=="CMta")

#Tajima's D

bol_tajima_neutral <- angsd_thetas_depth_neutral_36_bol %>%
  ggplot(aes(x=Tajima, fill=Era)) +
  labs(x="Tajima's D",y="Density") +
  theme_classic() +
  theme(legend.position="none") +
  geom_density(alpha=0.65) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  xlim(-2.5, 2.5) +
  ylim(0, 0.9)

mta_tajima_neutral <- angsd_thetas_depth_neutral_36_mta %>%
  ggplot(aes(x=Tajima, fill=Era)) +
  labs(x="Tajima's D",y="Density") +
  theme_classic() +
  theme(legend.position="none") +
  geom_density(alpha=0.65) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  xlim(-2.5, 2.5) +
  ylim(0, 0.9)


plot_grid(bol_tajima_neutral, mta_tajima_neutral,
          labels=c('A', 'B'))

#Subset by population
abol_notrans_neutral_36 <- subset(angsd_thetas_depth_neutral_36, Site=="ABol")
amta_notrans_neutral_36 <- subset(angsd_thetas_depth_neutral_36, Site=="AMta")
cbol_notrans_neutral_36 <- subset(angsd_thetas_depth_neutral_36, Site=="CBol")
cmta_notrans_neutral_36 <- subset(angsd_thetas_depth_neutral_36, Site=="CMta")

#Test to see if Watterson's theta is normally distributed

hist(abol_notrans_neutral_36$tW_bysite)
hist(amta_notrans_neutral_36$tW_bysite)
hist(cbol_notrans_neutral_36$tW_bysite)
hist(cmta_notrans_neutral_36$tW_bysite)

#Look's normally distributed

##Paired t-test

t.test(abol_notrans_neutral_36$tW_bysite, cbol_notrans_neutral_36$tW_bysite, paired=TRUE)
#Paired t-test
#data:  abol_notrans_neutral_36$tW_bysite and cbol_notrans_neutral_36$tW_bysite
#t = 36.172, df = 1914, p-value < 2.2e-16
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #0.005047505 0.005626223
#sample estimates:
  #mean difference
#0.005336864


t.test(amta_notrans_neutral_36$tW_bysite, cmta_notrans_neutral_36$tW_bysite, paired=TRUE)
#Paired t-test
#data:  amta_notrans_neutral_36$tW_bysite and cmta_notrans_neutral_36$tW_bysite
#t = 6.3556, df = 1914, p-value = 2.589e-10
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #0.000765754 0.001449265
#sample estimates:
  #mean difference
#0.00110751

#Bootstrapping of theta

library(boot)

x = as.vector(abol_notrans_neutral_36$tW_bysite)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

abol_boot = boot(x, samplemean, R=1000)

abol_boot
plot(abol_boot)
#original        bias     std. error
#t1* 0.09337913 -7.969634e-06 0.0002213565

boot.ci(boot.out=abol_boot, type="norm") #95% CI 0.0930, 0.0938
boot.ci(boot.out=abol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cbol_notrans_neutral_36$tW_bysite)
cbol_boot = boot(x, samplemean, R=1000)

cbol_boot
plot(cbol_boot)
#original       bias     std. error
#t1* 0.08804226 1.103813e-05 0.0002386666

boot.ci(boot.out=cbol_boot, type="norm") #95% CI 0.0876, 0.0885
boot.ci(boot.out=cbol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(amta_notrans_neutral_36$tW_bysite)
amta_boot = boot(x, samplemean, R=1000)

amta_boot
plot(amta_boot)
#original        bias     std. error
#t1* 0.09057536 -8.547077e-06 0.0002676779

boot.ci(boot.out=amta_boot, type="norm") #95% CI 0.0901, 0.0911
boot.ci(boot.out=amta_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cmta_notrans_neutral_36$tW_bysite)
cmta_boot = boot(x, samplemean, R=1000)

cmta_boot
#original      bias     std. error
#t1* 0.08946785 4.43883e-06 0.0002196294
boot.ci(boot.out=cmta_boot, type="norm") #95% CI 0.0890, 0.0899
boot.ci(boot.out=cmta_boot, type="bca") #Error estimated adjustment 'a' is NA

##Create data frame with means and 95% CI for plotting

population <- c("ABol", "CBol", "AMta", "CMta")
location <- c("Malampaya", "Malampaya", "Mantatao Island", "Mantatao Island")
mean_theta_neutral <- c(0.09337913, 0.08804226, 0.09057536, 0.08946785)
se_theta_neutral <- c(0.0002213565, 0.0002386666, 0.0002676779, 0.0002196294)
era <- c("Historical", "Contemporary", "Historical", "Contemporary")

theta_plotting_neutral <- data.frame(population, era, location, mean_theta_neutral, se_theta_neutral)
theta_plotting_neutral$population <- factor(theta_plotting_neutral$population, levels=c("ABol", "CBol", "AMta", "CMta"))
theta_plotting_neutral$location <- factor(theta_plotting_neutral$location, levels=c("Malampaya", "Mantatao Island"))
theta_plotting_neutral$era <- factor(theta_plotting_neutral$era, levels=c("Historical", "Contemporary"))

theta_neutral <- theta_plotting_neutral %>%
  ggplot(aes(x=era, y=mean_theta_neutral, color=location, shape=location)) +
  labs(x="Era",y="Mean Theta", color="Era") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_theta_neutral-1.96*se_theta_neutral, ymax=mean_theta_neutral+1.96*se_theta_neutral), width=0.0) +
  scale_color_manual(values = c("dodgerblue4", "darkseagreen2")) +
  ylim(0.06, 0.095)

##Nucleotide Diversity (pi)

#Test to see if theta is normally distributed

hist(abol_notrans_neutral_36$tP_bysite)
hist(amta_notrans_neutral_36$tP_bysite)
hist(cbol_notrans_neutral_36$tP_bysite)
hist(cmta_notrans_neutral_36$tP_bysite)

##Paired t-test

t.test(abol_notrans_neutral_36$tP_bysite, cbol_notrans_neutral_36$tP_bysite, paired=TRUE)
#Paired t-test

#t = -0.073435, df = 1914, p-value = 0.9415
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #-0.0003679528  0.0003413924
#sample estimates:
  #mean difference
#-1.328022e-05

t.test(amta_notrans_neutral_36$tP_bysite, cmta_notrans_neutral_36$tP_bysite, paired=TRUE)

#Paired t-test
#data:  amta_notrans_neutral_36$tP_bysite and cmta_notrans_neutral_36$tP_bysite
#t = 17.001, df = 1914, p-value < 2.2e-16
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #0.002596818 0.003274059
#sample estimates:
  #mean difference
#0.002935438

#Bootstrapping of pi

library(boot)

x = as.vector(abol_notrans_neutral_36$tP_bysite)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

abol_boot = boot(x, samplemean, R=1000)

abol_boot
plot(abol_boot)
#original        bias     std. error
#t1* 0.08579219 -1.908612e-06 0.0004680985

boot.ci(boot.out=abol_boot, type="norm") #95% CI 0.0849, 0.0867
boot.ci(boot.out=abol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cbol_notrans_neutral_36$tP_bysite)
cbol_boot = boot(x, samplemean, R=1000)

cbol_boot
plot(cbol_boot)
#original        bias     std. error
#t1* 0.08580547 -1.834982e-05 0.0005081295

boot.ci(boot.out=cbol_boot, type="norm") #95% CI 0.0848, 0.0868
boot.ci(boot.out=cbol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(amta_notrans_neutral_36$tP_bysite)
amta_boot = boot(x, samplemean, R=1000)

amta_boot
plot(amta_boot)
#original        bias    std. error
#t1* 0.08488255 -4.422332e-05 0.000481554

boot.ci(boot.out=amta_boot, type="norm") #95% CI 0.0840, 0.0859
boot.ci(boot.out=amta_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cmta_notrans_neutral_36$tP_bysite)
cmta_boot = boot(x, samplemean, R=1000)

cmta_boot
#original       bias     std. error
#t1* 0.08194711 6.487454e-06 0.0004913784
boot.ci(boot.out=cmta_boot, type="norm") #95% CI 0.0890, 0.0899
boot.ci(boot.out=cmta_boot, type="bca") #Error estimated adjustment 'a' is NA

##Create data frame with means and 95% CI for plotting

population <- c("ABol", "CBol", "AMta", "CMta")
location <- c("Malampaya", "Malampaya", "Mantatao Island", "Mantatao Island")
mean_pi_neutral <- c(0.08579219, 0.08580547, 0.08488255, 0.08194711)
se_pi_neutral <- c(0.0004680985, 0.0005081295, 0.000481554, 0.0004913784)
era <- c("Historical", "Contemporary", "Historical", "Contemporary")

pi_plotting_neutral <- data.frame(population, era, location, mean_pi_neutral, se_pi_neutral)
pi_plotting_neutral$population <- factor(pi_plotting_neutral$population, levels=c("ABol", "CBol", "AMta", "CMta"))
pi_plotting_neutral$era <- factor(pi_plotting_neutral$era, levels=c("Historical", "Contemporary"))
pi_plotting_neutral$location <- factor(pi_plotting_neutral$location, levels=c("Malampaya", "Mantatao Island"))

pi_neutral <- pi_plotting_neutral %>%
  ggplot(aes(x=era, y=mean_pi_neutral, color=location, shape=location)) +
  labs(x="Era",y="Mean Pi", color="Era") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_pi_neutral-1.96*se_pi_neutral, ymax=mean_pi_neutral+1.96*se_pi_neutral), width=0.0) +
  scale_color_manual(values = c("dodgerblue4", "darkseagreen2")) +
  ylim(0.06, 0.095)

plot_grid(theta_neutral, pi_neutral,
          labels=c('A', 'B'))

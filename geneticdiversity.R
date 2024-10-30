##Analyzing and Plotting genetic diversity metrics

library(tidyverse)
library(cowplot)
library(boot)

##Load in theta outputs by population (abol, amta, cbol, cmta)

abol_thetas_notrans <- read_table("abol_notrans.thetas.idx.pestPG")
amta_thetas_notrans <- read_table("amta_notrans.thetas.idx.pestPG")
cbol_thetas_notrans <- read_table("cbol_notrans.thetas.idx.pestPG")
cmta_thetas_notrans <- read_table("cmta_notrans.thetas.idx.pestPG")

#Add a column for site to each dataset

abol_thetas_notrans$Site <- "ABol"
amta_thetas_notrans$Site <- "AMta"
cbol_thetas_notrans$Site <- "CBol"
cmta_thetas_notrans$Site <- "CMta"

#Add a column for era to each dataset
abol_thetas_notrans$Era <- "Historical"
amta_thetas_notrans$Era <- "Historical"
cbol_thetas_notrans$Era <- "Contemporary"
cmta_thetas_notrans$Era <- "Contemporary"

#Merge datasets for plotting

angsd_thetas_notrans <- rbind(abol_thetas_notrans, cbol_thetas_notrans, amta_thetas_notrans, cmta_thetas_notrans)

##Import depth information

angsd_depth_notrans <- read_table("all_notrans.pos.gz")
#This is by site, we need it by chromosome to match with the theta calculations

#Take the average of total Depth for all the sites within a chromosome
angsd_avgdepth_notrans <- angsd_depth_notrans %>% group_by(chr) %>% summarise(avg_Depth = mean(totDepth))
#This gets us average total depth (summed across individuals) for each chromosome

#Add a column with average depth by individual
angsd_avgdepth_notrans$ind_depth <- angsd_avgdepth_notrans$avg_Depth/222

#Rename Chr column to match angsd_thetas dataframe
names(angsd_avgdepth_notrans)[1] <- "Chr"

#Merge with angsd_thetas dataframe

angsd_thetas_depth_notrans <- merge(angsd_thetas_notrans, angsd_avgdepth_notrans, by="Chr")

#Add a column with theta by site (theta for the contig divided by number of sites on the contig)
angsd_thetas_depth_notrans$tW_bysite <- angsd_thetas_depth_notrans$tW/angsd_thetas_depth_notrans$nSites

#Add a column with pi (nucleotide diversity) by site (theta for the contig divided by number of sites on the contig)
angsd_thetas_depth_notrans$tP_bysite <- angsd_thetas_depth_notrans$tP/angsd_thetas_depth_notrans$nSites

#Reorder for plotting
angsd_thetas_depth_notrans$Site <- factor(angsd_thetas_depth_notrans$Site, levels=c("ABol", "CBol", "AMta", "CMta"))

#Theta by Depth plot

angsd_thetas_depth_notrans %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

angsd_thetas_depth_notrans %>%
  ggplot(aes(x=ind_depth, y=tP_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

angsd_thetas_depth_notrans_bol <- subset(angsd_thetas_depth_notrans, Site=="ABol" | Site=="CBol")
angsd_thetas_depth_notrans_mta <- subset(angsd_thetas_depth_notrans, Site=="AMta" | Site=="CMta")

bol_plot_all <- angsd_thetas_depth_notrans_bol %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era), show.legend=FALSE) +
  labs(x="Mean Depth Per Individual",y="Theta") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(breaks=c(0.05, 0.075, 0.1, 0.125)) +
  ylim(0.025, 0.125)

mta_plot_all <- angsd_thetas_depth_notrans_mta %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Mean Depth Per Individual",y="Theta") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(breaks=c(0.05, 0.075, 0.1, 0.125)) +
  ylim(0.025, 0.125)

 bol_tp_all <- angsd_thetas_depth_notrans_bol %>%
  ggplot(aes(x=ind_depth, y=tP_bysite, color=Era), show.legend=FALSE) +
  labs(x="Mean Depth Per Individual",y="Nucleotide Diversity") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  ylim(0.02, 0.18)

mta_tp_all <- angsd_thetas_depth_notrans_mta %>%
  ggplot(aes(x=ind_depth, y=tP_bysite, color=Era)) +
  labs(x="Mean Depth Per Individual",y="Nucleotide Diversity") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  ylim(0.02, 0.18)

#Subset to depth 3-6X

angsd_thetas_depth_notrans_36 <- subset(angsd_thetas_depth_notrans, ind_depth >= 3 & ind_depth<=6)

angsd_thetas_depth_notrans_36 %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

angsd_thetas_depth_notrans_36_bol <- subset(angsd_thetas_depth_notrans_36, Site=="ABol" | Site=="CBol")
angsd_thetas_depth_notrans_36_mta <- subset(angsd_thetas_depth_notrans_36, Site=="AMta" | Site=="CMta")

bol_36_plot <- angsd_thetas_depth_notrans_36_bol %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(breaks=c(0.05, 0.075, 0.1, 0.125))

mta_36_plot <- angsd_thetas_depth_notrans_36_mta %>%
  ggplot(aes(x=ind_depth, y=tW_bysite, color=Era)) +
  labs(x="Average Individual Depth",y="Theta") +
  theme_bw() +
  theme(legend.position="none") +
  geom_point(size=1.5, alpha=0.5) +
  geom_smooth() +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(breaks=c(0.05, 0.075, 0.1, 0.125))

plot_grid(bol_plot, mta_plot, bol_36_plot, mta_36_plot, ncol=2,
          labels=c('A', 'B', 'C', 'D',
                    vjust=-1.5))

plot_grid(bol_plot, mta_plot,
          labels=c('A', 'B'))

plot_grid(bol_36_plot, mta_36_plot,
          labels=c('A', 'B'))


#Subset by population
abol_notrans_36 <- subset(angsd_thetas_depth_notrans_36, Site=="ABol")
amta_notrans_36 <- subset(angsd_thetas_depth_notrans_36, Site=="AMta")
cbol_notrans_36 <- subset(angsd_thetas_depth_notrans_36, Site=="CBol")
cmta_notrans_36 <- subset(angsd_thetas_depth_notrans_36, Site=="CMta")


#Test to see if theta is normally distributed

hist(abol_notrans_36$tP_bysite)
hist(amta_notrans_36$tP_bysite)
hist(cbol_notrans_36$tP_bysite)
hist(cmta_notrans_36$tP_bysite)

#Look's normally distributed

##Paired t-test

t.test(abol_notrans_36$tP_bysite, cbol_notrans_36$tP_bysite, paired=TRUE)
#Paired t-test

#data:  abol_notrans_36$tP_bysite and cbol_notrans_36$tP_bysite
#t = -42.344, df = 2290, p-value < 2.2e-16
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #-0.007941316 -0.007238323
#sample estimates:
 # mean difference
#-0.00758982

t.test(amta_notrans_36$tP_bysite, cmta_notrans_36$tP_bysite, paired=TRUE)
#Paired t-test

#data:  amta_notrans_36$tP_bysite and cmta_notrans_36$tP_bysite
#t = -47.459, df = 2290, p-value < 2.2e-16
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #-0.008333559 -0.007672206
#sample estimates:
  #mean difference
#-0.008002882

#Bootstrapping of pi

library(boot)

x = as.vector(abol_notrans_36$tP_bysite)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

abol_boot = boot(x, samplemean, R=1000)

abol_boot
plot(abol_boot)
#original       bias     std. error
#t1* 0.06590544 7.788715e-06 0.0003453374

boot.ci(boot.out=abol_boot, type="norm") #95% CI 0.0652, 0.0666
boot.ci(boot.out=abol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cbol_notrans_36$tP_bysite)
cbol_boot = boot(x, samplemean, R=1000)

cbol_boot
plot(cbol_boot)
#original       bias     std. error
#t1* 0.07349526 -6.64726e-06 0.0004091303

boot.ci(boot.out=cbol_boot, type="norm") #95% CI 0.0727, 0.0743
boot.ci(boot.out=cbol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(amta_notrans_36$tP_bysite)
amta_boot = boot(x, samplemean, R=1000)

amta_boot
plot(amta_boot)
#original        bias     std. error
#t1* 0.06437578 -8.452513e-06 0.0003441579

boot.ci(boot.out=amta_boot, type="norm") #95% CI 0.0637, 0.0651
boot.ci(boot.out=amta_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cmta_notrans_36$tP_bysite)
cmta_boot = boot(x, samplemean, R=1000)

cmta_boot
#original      bias     std. error
#t1* 0.07237866 2.45688e-05 0.0003814624
boot.ci(boot.out=cmta_boot, type="norm") #95% CI 0.0716, 0.0731
boot.ci(boot.out=cmta_boot, type="bca") #Error estimated adjustment 'a' is NA

##Create data frame with means and 95% CI for plotting

population <- c("ABol", "CBol", "AMta", "CMta")
location <- c("Malampaya", "Malampaya", "Mantatao Island", "Mantatao Island")
mean_pi <- c(0.06590544, 0.07349526, 0.06437578, 0.07237866)
se_pi <- c(0.0003453374, 0.0004091303, 0.0003441579, 0.0003814624)
era <- c("Historical", "Contemporary", "Historical", "Contemporary")

theta_plotting_pi <- data.frame(era, location, mean_pi, se_pi)
theta_plotting_pi$location <- factor(theta_plotting_pi$location, levels=c("Malampaya", "Mantatao Island"))
theta_plotting_pi$era <- factor(theta_plotting_pi$era, levels=c("Historical", "Contemporary"))

pi_all <- theta_plotting_pi %>%
  ggplot(aes(x=era, y=mean_pi, shape=location, color=location)) +
  labs(x="Era",y="Mean Pi", color="Era") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_pi-1.96*se_pi, ymax=mean_pi+1.96*se_pi), width=0.0) +
  ylim(0.06, 0.095) +
scale_color_manual(values = c("dodgerblue4", "darkseagreen2"))

#Test to see if theta is normally distributed

hist(abol_notrans_36$tW_bysite)
hist(amta_notrans_36$tW_bysite)
hist(cbol_notrans_36$tW_bysite)
hist(cmta_notrans_36$tW_bysite)

#Look's normally distributed

##Paired t-test

t.test(abol_notrans_36$tW_bysite, cbol_notrans_36$tW_bysite, paired=TRUE)
#Paired t-test

#data:  abol_notrans_36$tW_bysite and cbol_notrans_36$tW_bysite
#t = -9.6916, df = 2290, p-value < 2.2e-16
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #-0.001663957 -0.001103905
#sample estimates:
  #mean difference
#-0.001383931

t.test(amta_notrans_36$tW_bysite, cmta_notrans_36$tW_bysite, paired=TRUE)
#Paired t-test

#data:  amta_notrans_36$tW_bysite and cmta_notrans_36$tW_bysite
#t = -77.89, df = 2290, p-value < 2.2e-16
#alternative hypothesis: true mean difference is not equal to 0
#95 percent confidence interval:
  #-0.01446713 -0.01375655
#sample estimates:
  #mean difference
#-0.01411184

#Bootstrapping of theta

library(boot)

x = as.vector(abol_notrans_36$tW_bysite)

samplemean <- function(x, d) {
  return(mean(x[d]))
}

abol_boot = boot(x, samplemean, R=1000)

abol_boot
plot(abol_boot)
#original        bias     std. error
#t1* 0.07905219 -7.798933e-06 0.0002024426

boot.ci(boot.out=abol_boot, type="norm") #95% CI 0.0787,  0.0795
boot.ci(boot.out=abol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cbol_notrans_36$tW_bysite)
cbol_boot = boot(x, samplemean, R=1000)

cbol_boot
plot(cbol_boot)
#original       bias     std. error
#t1* 0.08043612 6.965127e-06 0.0001969585

boot.ci(boot.out=cbol_boot, type="norm") #95% CI 0.0800,  0.0808
boot.ci(boot.out=cbol_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(amta_notrans_36$tW_bysite)
amta_boot = boot(x, samplemean, R=1000)

amta_boot
plot(amta_boot)
#original        bias     std. error
#t1* 0.07443034 -1.082285e-05 0.0002415968

boot.ci(boot.out=amta_boot, type="norm") #95% CI 0.0740,  0.0749
boot.ci(boot.out=amta_boot, type="bca") #Error estimated adjustment 'a' is NA

x = as.vector(cmta_notrans_36$tW_bysite)
cmta_boot = boot(x, samplemean, R=1000)

cmta_boot
#original        bias     std. error
#t1* 0.08854218 -2.035286e-06 0.0001745086
boot.ci(boot.out=cmta_boot, type="norm") #95% CI 0.0882,  0.0889
boot.ci(boot.out=cmta_boot, type="bca") #Error estimated adjustment 'a' is NA

##Create data frame with means and 95% CI for plotting

population <- c("ABol", "CBol", "AMta", "CMta")
location <- c("Malampaya", "Malampaya", "Mantatao Island", "Mantatao Island")
mean_36 <- c(0.07905219, 0.08043612, 0.07443034, 0.08854218)
se_36 <- c(0.0002024426, 0.0001969585, 0.0002415968, 0.0001745086)
era <- c("Historical", "Contemporary", "Historical", "Contemporary")


theta_plotting_36 <- data.frame(population, era, location, mean_36, se_36)

theta_plotting_36$population <- factor(theta_plotting_36$population, levels=c("ABol", "CBol", "AMta", "CMta"))
theta_plotting_36$era <- factor(theta_plotting_36$era, levels=c("Historical", "Contemporary"))
theta_plotting_36$location <- factor(theta_plotting_36$location, levels=c("Malampaya", "Mantatao Island"))

theta_all <- theta_plotting_36 %>%
  ggplot(aes(x=era, y=mean_36, shape=location, color=location)) +
  labs(x="Era",y="Mean Theta", color="Era") +
  theme_classic() +
  theme(legend.position="none") +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_36-1.96*se_36, ymax=mean_36+1.96*se_36), width=0.0) +
  scale_color_manual(values = c("dodgerblue4", "darkseagreen2")) +
  ylim(0.06, 0.095)

plot_grid(theta_all, pi_all,
          labels=c('A', 'B'))

plot_grid(theta_all, pi_all, theta_neutral, pi_neutral,
          labels=c('A', 'B', 'C', 'D'))
#theta_neutral and pi_neutral are plotted in the geneticdiversity_neutral.R script

#Tajimas D Distribution Plot

angsd_thetas_depth_notrans_36 %>%
  ggplot(aes(x=Tajima, fill=Era)) +
  labs(x="Tajima's D",y="Density") +
  theme_bw() +
  geom_density(alpha=0.65) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))

bol_tajima_all <- angsd_thetas_depth_notrans_36_bol %>%
  ggplot(aes(x=Tajima, fill=Era)) +
  labs(x="Tajima's D",y="Density") +
  theme_classic() +
  theme(legend.position="none") +
  geom_density(alpha=0.65) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  xlim(-2.5, 2) +
  ylim(0, 1)

mta_tajima_all <- angsd_thetas_depth_notrans_36_mta %>%
  ggplot(aes(x=Tajima, fill=Era)) +
  labs(x="Tajima's D",y="Density") +
  theme_classic() +
  theme(legend.position="none") +
  geom_density(alpha=0.65) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  xlim(-2.5, 2) +
  ylim(0, 1)


plot_grid(bol_tajima_all, mta_tajima_all,
          labels=c('A', 'B'))

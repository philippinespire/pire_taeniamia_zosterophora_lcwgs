##Testing for Selection

#Adjusted chi-squared and CMH tests using the ACER package
library(ACER)

#Creating the frequency matrix
# MAFS
abol_mafs <- fread("abol_sites_notrans.mafs.gz", header=TRUE)
cbol_mafs <- fread("cbol_sites_notrans.mafs.gz", header=TRUE)
amta_mafs <- fread("amta_sites_notrans.mafs.gz", header=TRUE)
cmta_mafs <- fread("cmta_sites_notrans.mafs.gz", header=TRUE)

abol_cbol_sel <- abol_mafs$knownEM
abol_cbol_sel <- as.data.frame(abol_cbol_sel)
setnames(abol_cbol_sel, c("abol_cbol_sel"), c("abol"))
abol_cbol_sel$cbol <- cbol_mafs$knownEM

amta_cmta_sel <- as.data.frame(amta_mafs$knownEM)
setnames(amta_cmta_sel, c("amta_mafs$knownEM"), c("amta"))
amta_cmta_sel$cmta <- cmta_mafs$knownEM

#Coverage matrix
abol_pos <- fread("abol_sites_notrans.pos.gz", header=TRUE)
cbol_pos <- fread("cbol_sites_notrans.pos.gz", header=TRUE)
amta_pos <- fread("amta_sites_notrans.pos.gz", header=TRUE)
cmta_pos <- fread("cmta_sites_notrans.pos.gz", header=TRUE)

abol_cbol_cov_mat <- as.data.frame(abol_pos$totDepth)
setnames(abol_cbol_cov_mat, c("abol_pos$totDepth"), c("abol_totDepth"))
abol_cbol_cov_mat$cbol_totDepth <- cbol_pos$totDepth

amta_cmta_cov_mat <- as.data.frame(amta_pos$totDepth)
setnames(amta_cmta_cov_mat, c("amta_pos$totDepth"), c("amta_totDepth"))
amta_cmta_cov_mat$cmta_totDepth <- cmta_pos$totDepth

Ne <- c(1022, 1022)
gen <- c(0, 113)
rep <- c(1,2)


bol_pval <- adapted.chisq.test(freq=abol_cbol_sel, coverage=abol_cbol_cov_mat, Ne, gen=gen)
bol_pval <- as.data.frame(bol_pval)
sum(bol_pval$bol_pval<0.05) #419,676 SNPs with p-value <0.05
sum(bol_pval$bol_pval<0.01) #248,902 SNPs with p-value <0.01

bol_pval_vector <- as.numeric(bol_pval$bol_pval)
bol_pval$fdr <- p.adjust(bol_pval_vector, method='fdr', n=length(bol_pval_vector))
sum(bol_pval$fdr<0.05) #231,450 SNPs with q-value <0.05, 12.97% of SNPs

bol_freq <- cbind(bol_pval, all_sel)
bol_freq$chromo <- cmta_mafs$chromo
bol_freq$position <- cmta_mafs$position

bol_freq_fdr <- subset(bol_freq, fdr<0.05)
bol_freq_fdr$bol_change <- bol_freq_fdr$cbol - bol_freq_fdr$abol
sum(bol_freq_fdr$bol_change>0) #230,270
sum(bol_freq_fdr$cbol>0.25 & bol_freq_fdr$cbol<0.75) #33,755

mta_pval <- adapted.chisq.test(freq=amta_cmta_sel, coverage=amta_cmta_cov_mat, Ne, gen=gen)
mta_pval <- as.data.frame(mta_pval)
sum(mta_pval$mta_pval<0.05) #450,942 SNPs with p-value <0.05
sum(mta_pval$mta_pval<0.01) #306,265 SNPs with p-value <0.01

mta_pval_vector <- as.numeric(mta_pval$mta_pval)
mta_pval$fdr <- p.adjust(mta_pval_vector, method='fdr', n=length(mta_pval_vector))
sum(mta_pval$fdr<0.05) #275,831 SNPs with q-value <0.05, 15.45% of SNPs

mta_freq <- cbind(mta_pval, all_sel)
mta_freq$chromo <- cmta_mafs$chromo
mta_freq$position <- cmta_mafs$position

mta_freq_fdr <- subset(mta_freq, fdr<0.05)
mta_freq_fdr$mta_change <- mta_freq_fdr$cmta - mta_freq_fdr$amta
sum(mta_freq_fdr$mta_change>0) #274,219
sum(mta_freq_fdr$cmta>0.25 & mta_freq_fdr$cmta<0.75) #30,949

#Create matrices with both sampling sites (replicates for cmh test)

all_sel <- abol_cbol_sel
all_sel$amta <- amta_cmta_sel$amta
all_sel$cmta <- amta_cmta_sel$cmta
all_sel <- all_sel %>% relocate(amta, .after=abol)
all_sel <- as.matrix(all_sel)

all_cov_mat <- abol_cbol_cov_mat
all_cov_mat$amta_totDepth <- amta_cmta_cov_mat$amta_totDepth
all_cov_mat$cmta_totDepth <- amta_cmta_cov_mat$cmta_totDepth
all_cov_mat <- all_cov_mat %>% relocate(amta_totDepth, .after=abol_totDepth)
all_cov_mat <- as.matrix(all_cov_mat)

Ne <- 1022

cmh_pval <- adapted.cmh.test(freq=all_sel, coverage=all_cov_mat, Ne=rep(Ne, length(rep)), gen=gen, repl=rep)
cmh_pval <- as.data.frame(cmh_pval)
sum(cmh_pval$cmh_pval<0.05) #167,035 SNPs with p-value less than 0.05
sum(cmh_pval$cmh_pval<0.01) #108,839 SNPs with p-value less than 0.01

cmh_pval_vector <- as.numeric(cmh_pval$cmh_pval)
cmh_pval$fdr <- p.adjust(cmh_pval_vector, method='fdr', n=length(cmh_pval_vector))
sum(cmh_pval$fdr<0.05) #69,940 SNPs with q-value <0.05, 3.92% of SNPs

cmh_freq <- cbind(cmh_pval, all_sel)
cmh_freq$chromo <- cmta_mafs$chromo
cmh_freq$position <- cmta_mafs$position

cmh_freq_fdr <- subset(cmh_freq, fdr<0.05)
cmh_freq_fdr$bol_change <- cmh_freq_fdr$cbol - cmh_freq_fdr$abol
cmh_freq_fdr$mta_change <- cmh_freq_fdr$cmta - cmh_freq_fdr$amta

cmh_freq_fdr$results = ifelse(cmh_freq_fdr$bol_change > cmh_freq_fdr$mta_change, 'Bol',
                      ifelse(cmh_freq_fdr$bol_change < cmh_freq_fdr$mta_change, 'Mta', 'None'))
sum(cmh_freq_fdr$results=='Mta') #52,594 SNPs
sum(cmh_freq_fdr$results=='Bol') #17,337 SNPs

cmh_freq_fdr$results_abs = ifelse(abs(cmh_freq_fdr$bol_change) > abs(cmh_freq_fdr$mta_change), 'Bol',
                              ifelse(abs(cmh_freq_fdr$bol_change) < abs(cmh_freq_fdr$mta_change), 'Mta', 'None'))

sum(cmh_freq_fdr$results_abs =='Mta') #63,694 SNPs
sum(cmh_freq_fdr$results_abs =='Bol') #6237 SNPs

cmh_freq_fdr_pos <- subset(cmh_freq_fdr, bol_change>0 & mta_change>0) #4,733 SNPs
sum(cmh_freq_fdr_pos$results=='Mta') #4389 SNPs (92.73%)
sum(cmh_freq_fdr_pos$results=='Bol') #344 SNPs (7.27%)
sum(cmh_freq_fdr_pos$cbol>0.25 & bol_freq_fdr$cbol<0.75)

cmh_freq_fdr_neg <- subset(cmh_freq_fdr, bol_change<0 & mta_change<0) #17,334 SNPs
cmh_freq_fdr_bolpos <- subset(cmh_freq_fdr, bol_change>=0 & mta_change<=0) #2285 SNPs
cmh_freq_fdr_mtapos <- subset(cmh_freq_fdr, bol_change<=0 & mta_change>=0) #45597 SNPs

#Change the allele to the one that increases over time for cmh_freq_fdr_neg
setnames(cmh_freq_fdr_neg, c("abol"), c("abol_orig"))
setnames(cmh_freq_fdr_neg, c("amta"), c("amta_orig"))
setnames(cmh_freq_fdr_neg, c("cbol"), c("cbol_orig"))
setnames(cmh_freq_fdr_neg, c("cmta"), c("cmta_orig"))

cmh_freq_fdr_neg$abol <- 1-cmh_freq_fdr_neg$abol_orig
cmh_freq_fdr_neg$cbol <- 1-cmh_freq_fdr_neg$cbol_orig
cmh_freq_fdr_neg$amta <- 1-cmh_freq_fdr_neg$amta_orig
cmh_freq_fdr_neg$cmta <- 1-cmh_freq_fdr_neg$cmta_orig

cmh_freq_fdr_neg <- select(cmh_freq_fdr_neg, -3:-6)
cmh_freq_fdr_neg <- cmh_freq_fdr_neg %>% relocate(abol, .after=fdr)
cmh_freq_fdr_neg <- cmh_freq_fdr_neg %>% relocate(amta, .after=abol)
cmh_freq_fdr_neg <- cmh_freq_fdr_neg %>% relocate(cbol, .after=amta)
cmh_freq_fdr_neg <- cmh_freq_fdr_neg %>% relocate(cmta, .after=cbol)

cmh_freq_fdr_bind <- rbind(cmh_freq_fdr_neg, cmh_freq_fdr_pos)
cmh_freq_fdr_bind <- rbind(cmh_freq_fdr_bind, cmh_freq_fdr_bolpos)
cmh_freq_fdr_bind <- rbind(cmh_freq_fdr_bind, cmh_freq_fdr_mtapos)

cmh_freq_fdr_bind$bol_change <- cmh_freq_fdr_bind$cbol - cmh_freq_fdr_bind$abol
cmh_freq_fdr_bind$mta_change <- cmh_freq_fdr_bind$cmta - cmh_freq_fdr_bind$amta


cmh_freq_fdr_bind$results_abs = ifelse(abs(cmh_freq_fdr_bind$bol_change) > abs(cmh_freq_fdr_bind$mta_change), 'Bol',
                                  ifelse(abs(cmh_freq_fdr_bind$bol_change) < abs(cmh_freq_fdr_bind$mta_change), 'Mta', 'None'))

sum(cmh_freq_fdr_bind$results_abs =='Mta') #63,694 SNPs
sum(cmh_freq_fdr_bind$results_abs =='Bol') #6237 SNPs

cmh_freq_fdr_bind_pos <- subset(cmh_freq_fdr_bind, bol_change>=0 & mta_change>=0) #30,318 SNPs
sum(cmh_freq_fdr_bind_pos$results_abs=='Mta') #27336 SNPs (90.16%)
sum(cmh_freq_fdr_bind_pos$results_abs=='Bol') #2964 SNPs


##Plotting

cmh_bol_bind <- as.data.frame(cmh_freq_fdr_bind$abol)
cmh_bol_bind$abol <- cmh_freq_fdr_bind$abol
cmh_bol_bind$cbol <- cmh_freq_fdr_bind$cbol
cmh_bol_bind <- pivot_longer(cmh_bol_bind, abol:cbol)

cmh_bol <- cmh_bol_bind %>%
  ggplot(aes(x=value, color=name)) +
  labs(x="Allele Frequency",y="Density") +
  theme_classic() +
  theme(legend.position="none") +
  geom_density(size=1) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  ylim(0, 1)

cmh_mta_bind <- as.data.frame(cmh_freq_fdr_bind$amta)
cmh_mta_bind$amta <- cmh_freq_fdr_bind$amta
cmh_mta_bind$cmta <- cmh_freq_fdr_bind$cmta
cmh_mta_bind <- pivot_longer(cmh_mta_bind, amta:cmta)

cmh_mta <- cmh_mta_bind %>%
  ggplot(aes(x=value, color=name)) +
  labs(x="Allele Frequency",y="Density") +
  theme_classic() +
  theme(legend.position="none") +
  geom_density(size=1) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) +
  ylim(0, 1)

plot_grid(cmh_bol, cmh_mta,
          labels=c('A', 'B'))


##Creating the neutral SNP list

##Read in original SNP list

snp_list <- fread("global_snp_list_depth1_15.txt", header=FALSE)
snp_list$chrpos <- paste(snp_list$V1, snp_list$V2, sep="_")
cmh_freq_fdr$chrpos <- paste(cmh_freq_fdr$chromo, cmh_freq_fdr$position, sep="_")
mta_freq_fdr$chrpos <- paste(mta_freq_fdr$chromo, mta_freq_fdr$position, sep="_")
bol_freq_fdr$chrpos <- paste(bol_freq_fdr$chromo, bol_freq_fdr$position, sep="_")

all_freq_fdr <- as.vector(c(cmh_freq_fdr$chrpos, mta_freq_fdr$chrpos, bol_freq_fdr$chrpos))
length(all_freq_fdr[duplicated(all_freq_fdr)]) #96722 duplicates
#577221 - 96722 duplicates = 480449 unique SNPs to be removed

snp_list_nosel <- snp_list[!snp_list$chrpos %in% all_freq_fdr,]
#Removed 480,449 SNPs- perfect! This is the neutral SNP list to re-do theta and pop structure with.

#Check that number of chromosomes is the same
length(unique(snp_list_nosel$V1))
#Yes number of chromosomes (contigs) is still 11150 so no need to update chromosome list

#Export neutral SNP list as a txt file
snp_list_nosel_txt <- snp_list_nosel[,1:4]


write.table(snp_list_nosel_txt, "global_snp_list_neutral.txt", sep="\t", row.names=FALSE,
            col.names=FALSE)

##PCAngsd Admixture Results##

library(pophelper) #v.2.3.1

#K=2

k2_angsd_not <- read.table("angsd_admix_notrans.admix.2.Q")
k2_angsd_not <- as.data.frame(k2_angsd_not)

##Add in pop labels

meta.data <- data.frame(loc=pop_label_angsd)
meta.data$loc <- as.character(meta.data$loc)

q2_not <- list(k2_angsd_not)
plotQ(as.qlist(q2_not), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(), dpi=1000,
      clustercol = c("#F8766D", "#00BFC4"),
      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1)

#K=3

k3_angsd_not <- read.table("angsd_admix_notrans_k3.admix.3.Q")

q3_not <- list(k3_angsd_not)
plotQ(as.qlist(q3_not), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(),
      clustercol = c("#9999FF", "#2121D9", "#FF9329"),
      showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1)

#K=4

k4_angsd_not <- read.table("angsd_admix_notrans_k4.admix.4.Q")

q4_not <- list(k4_angsd_not)
plotQ(as.qlist(q4_not), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(),
      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
      showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1)

#K=5

k5_angsd_not <- read.table("angsd_admix_notrans_k5.admix.5.Q")

q5_not <- list(k5_angsd_not)
plotQ(as.qlist(q5_not), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(),
      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
      showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1)


##PCAngsd Admixture Results##

#Plot K=2

library(pophelper) #v.2.3.1

k2_neutral <- read.table("angsd_allpop_neutral_default.admix.2.Q")

##Add in pop labels

meta.data <- data.frame(loc=pop_label_angsd)
meta.data$loc <- as.character(meta.data$loc)

q2_neutral <- list(k2_neutral)
plotQ(as.qlist(q2_neutral), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(), dpi=1000,
      clustercol = c("#F8766D", "#00BFC4"),
      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = FALSE, grplabspacer = 0.1)

#K=3

k3_neutral <- read.table("angsd_allpop_neutral_k3.admix.3.Q")

q3_neutral <- list(k3_neutral)
plotQ(as.qlist(q3_neutral), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(),
      clustercol = c("#9999FF", "#2121D9", "#FF9329"),
      showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1)

#K=4

k4_neutral <- read.table("angsd_allpop_neutral_k4.admix.4.Q")

q4_neutral <- list(k4_neutral)
plotQ(as.qlist(q4_neutral), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(),
      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
      showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1)

#K=5

k5_neutral <- read.table("angsd_allpop_neutral_k5.admix.5.Q")

q5_neutral <- list(k5_neutral)
plotQ(as.qlist(q5_neutral), imgoutput = "sep", returnplot = TRUE, exportpath=getwd(),
      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
      showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1)


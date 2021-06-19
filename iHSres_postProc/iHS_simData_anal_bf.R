## Packages
library(dplyr)
library(ggplot2)

source("~/Documents/kaaj/selectSign/selSweeps/scripts/iHSres_postProc/propHiIHS_wind.R")


# Set parent result folder as the working directory 
setwd("~/Documents/kaaj/selectSign/selSweeps/selScan/simulations/results/")

iHS_sim_sweep_chr1_500N = data.frame(read.table("../sweep_drift_40Mb/sweep_chr1_500Ne_100ind_40Mb_biAll.ihs.out.100bins.norm.final", header = FALSE, skip = 0))
colnames(iHS_sim_sweep_chr1_500N) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")

iHS_sim_sweep_chr2_500N = data.frame(read.table("../sweep_drift_40Mb/sweep_chr2_500Ne_100ind_40Mb_biAll.ihs.out.100bins.norm.final", header = FALSE, skip = 0))
colnames(iHS_sim_sweep_chr2_500N) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")

iHS_simul_driftOnly_chr1_5kN = data.frame(read.table("driftOnly_chr1_5kNe_100ind_20Mb_biAll.ihs.out.100bins.norm.final", header = FALSE, skip = 0))
colnames(iHS_simul_driftOnly_chr1_5kN) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")

iHS_simul_driftOnly_chr2_5kN = data.frame(read.table("driftOnly_chr2_5kNe_100ind_20Mb_biAll.ihs.out.100bins.norm.final", header = FALSE, skip = 0))
colnames(iHS_simul_driftOnly_chr2_5kN) <- c("locID", "chrPos", "freq_d","iHH_d","iHH_a","iHS_unst","iHS","qScore","cmPos")



# Create an extra-column containing |iHS| scores
iHS_sim_sweep_chr1_500N = iHS_sim_sweep_chr1_500N %>% mutate(iHS_mod = Mod(iHS))
iHS_sim_sweep_chr2_500N = iHS_sim_sweep_chr2_500N %>% mutate(iHS_mod = Mod(iHS))

iHS_simul_driftOnly_chr1_5kN = iHS_simul_driftOnly_chr1_5kN %>% mutate(iHS_mod = Mod(iHS))
iHS_simul_driftOnly_chr2_5kN = iHS_simul_driftOnly_chr2_5kN %>% mutate(iHS_mod = Mod(iHS))
# iHS_simul_driftOnly_10kN_chr1 = iHS_simul_driftOnly_10kN %>% filter(chrPos < 1000000) %>% mutate(iHS_mod = Mod(iHS))
# iHS_simul_driftOnly_10kN_chr2 = iHS_simul_driftOnly_10kN %>% filter(chrPos > 1000000) %>% mutate(iHS_mod = Mod(iHS))

# Calculate the 99% quantile value of the |iHS|
val_top1_sweep_chr1_500N<-unname(quantile(iHS_sim_sweep_chr1_500N$iHS_mod,.99))
val_top1_sweep_chr2_500N<-unname(quantile(iHS_sim_sweep_chr2_500N$iHS_mod,.99))

val_top1_chr1_5kN<-unname(quantile(iHS_simul_driftOnly_chr1_5kN$iHS_mod,.99))
val_top1_chr2_5kN<-unname(quantile(iHS_simul_driftOnly_chr2_5kN$iHS_mod,.99))
# Select all the |iHS| scores that are higher than the 99% quantile value (top 1%)

iHS_top1_sweep_chr1_500N = iHS_sim_sweep_chr1_500N %>% filter(iHS_mod > val_top1_sweep_chr1_500N)
iHS_top1_sweep_chr2_500N = iHS_sim_sweep_chr2_500N %>% filter(iHS_mod > val_top1_sweep_chr2_500N)

iHS_top1_simul_chr1_5kN = iHS_simul_driftOnly_chr1_5kN %>% filter(iHS_mod > val_top1_chr1_5kN)
iHS_top1_simul_chr2_5kN = iHS_simul_driftOnly_chr2_5kN %>% filter(iHS_mod > val_top1_chr2_5kN)

propHiiHSwind_sweep_chr1_500N = find_HiiHS_SNPclusters(iHS_sim_sweep_chr1_500N,.5)
propHiiHSwind_sweep_chr2_500N = find_HiiHS_SNPclusters(iHS_sim_sweep_chr2_500N,.5)

propHiiHSwind_simul_chr1_5kN = find_HiiHS_SNPclusters(iHS_simul_driftOnly_chr1_5kN,.5)
propHiiHSwind_simul_chr2_5kN = find_HiiHS_SNPclusters(iHS_simul_driftOnly_chr2_5kN,.5)

# Filters the windows that contains at least 10% hiIHS SNPs and 20 SNPs
hiIHSwind_simul_chr1_500N_10p = propHiiHSwind_simul_chr1_500N %>% filter (propTop1 > 0.1 & nSNPs > 10)
hiIHSwind_simul_chr2_500N_10p = propHiiHSwind_simul_chr2_500N %>% filter (propTop1 > 0.1 & nSNPs > 10)

hiIHSwind_simul_chr1_5kN_10p = propHiiHSwind_simul_chr1_5kN %>% filter (propTop1 > 0.1 & nSNPs > 10)
hiIHSwind_simul_chr2_5kN_10p = propHiiHSwind_simul_chr2_5kN %>% filter (propTop1 > 0.1 & nSNPs > 10)

# Filters top 1% of the windows based on the proportion of hiIHS SNPs
xIHSwind_sweep_chr1_500N_top1 = propHiiHSwind_sweep_chr1_500N %>% filter(propTop1 > quantile(propTop1,.99) & nSNPs > 10)
xIHSwind_sweep_chr2_500N_top1 = propHiiHSwind_sweep_chr2_500N %>% filter(propTop1 > quantile(propTop1,.99) & nSNPs > 10)

hiIHSwind_simul_chr1_5kN_top1 = propHiiHSwind_simul_chr1_5kN %>% filter(propTop1 > quantile(propTop1,.99) & nSNPs > 10)
hiIHSwind_simul_chr2_5kN_top1 = propHiiHSwind_simul_chr2_5kN %>% filter(propTop1 > quantile(propTop1,.99) & nSNPs > 10)
# create the ggplot Object with chromosome positions as X and top 1% |iHS| as Y
# plotObj <- ggplot(data=tempDat, aes(x = chrPos/1000000, y = iHS_mod, colour=pOp)) 
ggplot() +  geom_rect(data = xIHSwind_sweep_chr1_500N_top1, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_sweep_chr1_500N, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = xIHSwind_sweep_chr2_500N_top1, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_sweep_chr2_500N, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr1_5kN_top1, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr1_5kN, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr2_5kN_top1, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr2_5kN, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr1_500N_10p, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr1_500N, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr2_500N_10p, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr2_500N, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr1_5kN_10p, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr1_5kN, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr2_5kN_10p, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr2_5kN, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")



hist(propHiiHSwind_simul_chr1_500N$propTop1, xlab = "Proportion of extreme IHS SNPs", main = "Distribution of extreme IHS SNPs proportions")
hist(propHiiHSwind_simul_chr2_500N$propTop1, xlab = "Proportion of extreme IHS SNPs", main = "Distribution of extreme IHS SNPs proportions")

hist(propHiiHSwind_simul_chr1_5kN$propTop1, xlab = "Proportion of extreme IHS SNPs", main = "Distribution of extreme IHS SNPs proportions")
hist(propHiiHSwind_simul_chr2_5kN$propTop1, xlab = "Proportion of extreme IHS SNPs", main = "Distribution of extreme IHS SNPs proportions")

ggplot() +  geom_rect(data = hiIHSwind_simul_chr2_top1, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr2, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

ggplot() +  geom_point(data=iHS_top1_simul_chr1, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  geom_rect(data = hiIHSwind_simul_chr1_top1, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")


ggplot() +  geom_rect(data = hiIHSwind_simul_chr1_10p, aes(xmin = fSNP/1000000, xmax = lSNP/1000000, ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_point(data=iHS_top1_simul_chr1, aes(x = chrPos/1000000, y = iHS_mod), size=.1) +
  theme_bw() + xlab("Position (Mb)") + ylab("|iHS|")

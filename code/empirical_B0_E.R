setwd('~/Documents/div_cue_paper/code')

graphics.off()
library(minpack.lm)
library(ggplot2)
library(reshape2)
library(grid)

###### Experimental growth and res data, Tr = 0 ########
###### but with only 33 strains in total #########
###### smith2020systematic ######
data = read.csv("../data/aerobic_tpc_data.csv")
data = data[-which(data$temps_before_peak_resp<=3 | data$temps_before_peak_growth<=3),]
res_data = data.frame(data$E_resp, data$B0_resp, data$E_D_resp, data$Tpk_resp, data$Ppk_resp, data$r_sq_resp)
names(res_data) = c('Ea', 'B0', 'E_D', 'T_pk', 'P_pk','r_sq')
grow_data = data.frame(data$E_growth, data$B0_growth, data$E_D_growth, data$Tpk_growth, data$Ppk_growth, data$r_sq_growth)
names(grow_data) = c('Ea', 'B0', 'E_D', 'T_pk', 'P_pk','r_sq')

mean_Br = mean(grow_data$B0)
mean_Er = mean(grow_data$Ea)
var_Br = var(log(grow_data$B0))/abs(mean(log(grow_data$B0))) # standardizing with mean value
var_Er = var(grow_data$Ea)/abs(mean(grow_data$Ea)) # standardizing with mean value
# correlation coefficient
rho_r = cov(log(grow_data$B0), grow_data$Ea)/(sqrt(var(log(grow_data$B0))) * sqrt(var(grow_data$Ea)))

mean_Bm = mean(res_data$B0)
mean_Em = mean(res_data$Ea)
var_Bm = var(log(res_data$B0))/abs(mean(log(res_data$B0))) # standardizing with mean value
var_Em = var(res_data$Ea)/abs(mean(res_data$Ea)) # standardizing with mean value
# correlation coefficient
rho_m = cov(log(res_data$B0), res_data$Ea)/(sqrt(var(log(res_data$B0))) * sqrt(var(res_data$Ea)))

# Carbon use efficiency
CUE_0 = grow_data$B0/(grow_data$B0+res_data$B0)
mean(CUE_0)

###### Metadata of growth ########
###### smith2019community ######
###### reanalysed in clegg2022variation, changed Tr = 10 #######
summ = read.csv("../data/summary.csv")
s_data = summ[summ$source == "meta", ]
smean_Br = mean(s_data$B0_ss)
smean_Er = mean(s_data$E)
svar_Br = var(log(s_data$B0_ss))/mean(log(s_data$B0_ss))
svar_Er = var(s_data$E)/mean(s_data$E)
# covariance coefficient
srho_r = cov(log(s_data$B0_ss), s_data$E)/(sqrt(var(log(s_data$B0_ss)))*sqrt(var(s_data$E)))

png('../result/B_Ea_empirical.png', width = 800, height = 600)
ggplot(s_data, aes(x=B0, y=E)) + labs(x = expression("log(B"[0]*")"), y = "E")+
  geom_point(size=4, shape=16,color = "black", alpha = 0.7) + theme_bw() + theme(panel.background = element_blank(), text = element_text(size = 35))
dev.off()

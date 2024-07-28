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
rho_r = cor(log(grow_data$B0), grow_data$Ea)

mean_Bm = mean(res_data$B0)
mean_Em = mean(res_data$Ea) 
median(res_data$Ea)
var_Bm = var(log(res_data$B0))/abs(mean(log(res_data$B0))) # standardizing with mean value
var_Em = var(res_data$Ea)/abs(mean(res_data$Ea)) # standardizing with mean value
# correlation coefficient
rho_m = cor(log(res_data$B0), res_data$Ea)

### Calc B0
Schoolfield = function(Temp, B0, T_pk, Ea, E_D){
  return(B0*exp(-Ea*(1/Temp-1/273.15)/k)/(1+(Ea/(E_D-Ea))*exp(E_D*(1/T_pk-1/Temp)/k)))
}

Bm = Schoolfield(283.15, res_data$B0, res_data$T_pk, res_data$Ea, res_data$E_D)
hist(Bm) 
mean(log (Bm))

CUE = function(G, R){
  return(G/(G+R))
}
CUE10 = c()
for(i in 1:nrow(res_data)){
  G = Schoolfield(Temp = 283.15, B0 = grow_data$B0[i], T_pk = grow_data$T_pk[i], Ea = grow_data$Ea[i], E_D = grow_data$E_D[i])
  R = Schoolfield(Temp = 283.15, B0 = res_data$B0[i], T_pk = res_data$T_pk[i], Ea = res_data$Ea[i], E_D = res_data$E_D[i])
  CUE10 = c(CUE10, G/(G+R))
}
# Carbon use efficiency
BCUE = median(CUE10)

L = 0.3
Bu = Bm / (1 - L - BCUE)
mean(log(Bu))

###### Metadata of growth ########
###### smith2019community ######

temp = 0:30+273.15
summ = read.csv("../data/summary.csv")
  
s_data = summ[summ$source == "growth", ]
m_data = summ[summ$source == "res", ]

smean_Br = mean(log(s_data$B0_ss))
median(log(s_data$B0_ss))
smean_Er = mean(s_data$E)
mean(s_data$E)
svar_Br = var(log(s_data$B0_ss))/mean(log(s_data$B0_ss))
svar_Er = var(s_data$E)/mean(s_data$E)
# covariance coefficient
srho_r = cor(log(s_data$B0_ss), s_data$E)

mean(m_data$E)

mean(m_data$B0_ss)
svar_Bm = var(log(m_data$B0_ss))/mean(log(m_data$B0_ss))
svar_Em = var(m_data$E)/mean(m_data$E)
# covariance coefficient
srho_m = cor(log(m_data$B0_ss), m_data$E)
Bu = mean(m_data$B0_ss) / (1 - L - BCUE)

median(s_data$Tpk)
median(r_data$Tpk)

Schoolfield = function(Temp, B0, T_pk, Ea, E_D){
  return(B0*exp(-Ea*(1/Temp-1/283.15)/k)/(1+(Ea/(E_D-Ea))*exp(E_D*(1/T_pk-1/Temp)/k)))
}

temps <- 273.15 + seq(0, 30, by = 1)
plot_data <- data.frame()
for(i in 1:nrow(s_data)){
  temp_p = log(Schoolfield(Temp = 273.15:303.15, B0 = s_data$B0_ss[i], T_pk = 273.15+s_data$Tpk[i], Ea = s_data$E[i], E_D = s_data$Ed[i]))
  temp_df <- data.frame(Temperature = temps, GrowthRate = temp_p, Species = as.factor(i))
  plot_data <- bind_rows(plot_data, temp_df)}

png('../result/TPC_empirical.png', width = 1200, height = 450)
ggplot(plot_data, aes(x = Temperature, y = GrowthRate, color = Species)) +
  geom_line(alpha = 0.4, size = 1.2) +
  labs(x = "Temperature (°C)",
       y = "Growth (log)",
       color = "Species") +
  scale_x_continuous(breaks = seq(273.15, 303.15, by = 5), labels = function(x) sprintf("%.0f", x - 273.15)) +  # Convert Kelvin to °C
  scale_y_continuous(limits = c(-20, -5)) +
  theme_bw()+
  theme(panel.background = element_blank(),  # Clean panel background
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        legend.position = "none",
        text = element_text(size = 50))
dev.off()

png('../result/B_Ea_empirical.png', width = 1200, height = 450)
ggplot(s_data, aes(x=log(B0_ss), y=E)) + labs(x = expression("log(B"[r]*")"), y = expression("E"[r]*""))+
  geom_point(size=4, shape=16,color = "black", alpha = 0.7) + theme_bw() + theme(panel.background = element_blank(), text = element_text(size = 50))
dev.off()
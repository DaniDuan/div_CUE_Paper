
EMP_data = read.csv("../data/EMP.csv", header = T)
names(EMP_data) = c("ID", "Sample", "Temp", "pH", "Richness")
EMP_data = subset(EMP_data,Temp != "" & Temp <= 30 & Temp > 0)
scatter.smooth(EMP_data$Temp, EMP_data$Richness)

m = max(EMP_data$Richness)
EMP_data$relative_rich = EMP_data$Richness/m
scatter.smooth(EMP_data$Temp, EMP_data$relative_rich)

write.csv(EMP_data, "../data/EMP_filtered.csv", row.names=T)


quantile(EMP_data$relative_rich, probs = 0.99)

quan = c()
for(i in unique(EMP_data$Temp)){
  sub = subset(EMP_data, Temp == i)
  quan = c(quan, quantile(sub$relative_rich, probs = 0.99))
}

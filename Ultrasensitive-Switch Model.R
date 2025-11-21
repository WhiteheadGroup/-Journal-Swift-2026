#Ultrasensitive-Switch Model
#SDS

#Set Working Directory
setwd("/Users/samdswift/Desktop")
library(ggplot2)

#Set Initiaal Conditions
K_d1 = 115.6 #nM
K_dT = 2.1 #M
K_d2 = 105 #nM
P_0 = 8 #nM
R_0 = 960 #nM

#Without Titrant
T_0 = 0
Y = NULL
for (C in c(seq(0.01,7.99,0.01))) { 
  x = (K_d2*C)/(P_0-C)
  L = (K_d1*x)/(R_0-C-x-(x*T_0)/(x+K_dT))
  PL = x
  S = (x*T_0)/(x+K_dT)  
  Ligand_norm = L + PL + C + S
  newrow = data.frame(C,Ligand_norm)
  Y = rbind(Y,newrow)
}

#With Titrant
T_0 = 800 #nM
Z = NULL
for (C in c(seq(0.01,7.99,0.01))) { 
  x = (K_d2*C)/(P_0-C)
  L = (K_d1*x)/(R_0-C-x-(x*T_0)/(x+K_dT))
  PL = x
  S = (x*T_0)/(x+K_dT)  
  Ligand_titrated = L + PL + C + S
  newrow = data.frame(Ligand_titrated)
  Z = rbind(Z,newrow)
}

Graph = cbind(Y,Z)

dev.new()
pdf(file = "Ultrasensitive-Switch Model.pdf")
ggplot(data = Graph) + 
  geom_path(mapping = aes(x = log10(Ligand_norm), y = C,colour= "- Titrant"), size = 1.5) +
  geom_path(mapping = aes(x = log10(Ligand_titrated), y = C,colour= "+ Titrant"), size = 1.5) +
  xlab("Ligand") +
  ylab("Complex") +
  ggtitle("Ultrasensitive-Switch Model") +
  xlim(0,5)+
  ylim(0,8)+
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, face="bold"),
    axis.title.x = element_text(size=15, face="bold"),
    axis.title.y = element_text(size=15, face="bold"))
dev.off()
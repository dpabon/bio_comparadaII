## Grafico de simulaciones
library(ggplot2)

datos <- read.csv("~/MEGAsync/bio_comparadaII/result/simulations/simulations_Graph.csv")

## Box Plot
caja <- ggplot(data = datos, aes(x= MODELO, y= AIC_D, fill=MODELO)) +
  geom_boxplot() +
  scale_x_discrete(name= "Modelo") +
  scale_y_continuous(name = "Delta AIC") +
  theme(axis.text.x = element_text(colour="grey20",size=11,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=0,face="plain"),
        title = element_text(size=16,angle=0,hjust=.5,vjust=.5,face="plain"),
        legend.text=element_text(size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        legend.position="none")

caja



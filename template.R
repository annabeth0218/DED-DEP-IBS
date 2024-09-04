# install.packages("tidyverse")
library(tidyverse)
read.xlsx("filename.xlsx")
X2sample_CoJo$n<-as.factor(X2sample_CoJo$n)
ggplot(data = X2sample_CoJo, aes(x = n , y = m ,fill = group))+
  geom_bar(stat = "identity", position = position_dodge())+
  geom_hline(yintercept = 3, color = "black")+
  labs(x="number of sample", y="-log(P value)")+
  scale_fill_manual(values = c( "#FFC900", "#FF9000", "#22D222", "#005800", "#D879FF", "#7100A4"))


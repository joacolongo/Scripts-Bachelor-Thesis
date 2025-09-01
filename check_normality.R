# Autor Joaquin Algorta Bove
# check_normality
# Funcion:  hace un test shapiro-wilk para comprobar la normalidad de los datos de color de la población
#           crea graficos de densidad para comprobar el aspecto de la distribución
#           en ellos pinta los valores de Endural y Aldura en nuestro caso

library(tidyverse)

####  cosecha 2019
x<-1:105
#datos fenotipicos de la ponblación
b_2019<-c(17.84,17.35,16.47,17.96,18.28,18.74,17.11,17.58,17.45,17.17,
          17.26,17.24,17.8,17.35,17.96,16.72,18.19,17.54,16.63,17.26,
          17.51,17.71,17.24,17.99,17.18,19.05,17.28,17.23,17.38,15.77,
          17.34,17.07,17.45,15.32,17.14,16.53,17.21,18,19.1,17.34,17.11,
          17.17,15.85,17,15.83,17.34,17.25,17.65,17.36,16.56,16.99,17.21,
          16.89,17.4,16.76,16.93,18.07,17.99,16.54,18.36,18.49,16.64,16.97,
          18.01,16.06,16.53,16.56,16.6,18.22,16.26,19.58,17.69,18.72,16.98,
          16.89,16.14,17.75,15.92,17.67,18.05,18.04,17.48,16.84,18.19,17.91,
          18.37,16.11,15.77,14.83,17.37,17.5,16.11,19.09,16.64,18.4,17.75,
          17.98,17.13,18.41,18.21,17.38,18.25,17.43,17.59,16.68)
matrix2019<-as.data.frame(cbind(x,b_2019))
#datos de los parentales
Endural2019<-16.47
Aldura2019<-19.38


#### cosecha 2020
x<-1:106
#datos fenotípicos de la población
b_2020<-c(17.38,16.85,15.5,16.47,16.75,15.93,15.64,15.89,15.32,16.83,
          17.33,17.57,16.02,16.93,17.16,15.96,15.75,17.14,15.59,15.37,
          16.56,16.03,16.03,16.97,17.43,15.69,15.67,14.89,15.26,14.93,
          16.13,16.34,18.3,15.25,15.87,14.7,17.04,16.07,16.69,16.31,
          16.66,16.83,15.22,15.92,16.18,16.7,15.71,16,14.85,14.85,15.22,
          16,15.09,17.87,17.21,16.18,16.34,15.29,15.17,14.94,16.66,16.89,
          16.13,15.58,16.48,14.96,14.49,16.04,15.39,16.34,16.94,17.06,14.94,
          16.45,16.09,15.08,14.92,15.7,17.48,17.56,16.89,16.32,15.46,16.03,
          16.26,16.24,17.17,16.23,14.93,14.54,16.97,14.39,16.25,16.43,
          14.56,16.82,17.39,16.78,15.44,16.82,16.26,16.97,15.79,16.48,16.34,16.66)
matrix2020<-as.data.frame(cbind(x,b_2020))
#datos de los parentales
Endural2020<-14.54
Aldura2020<-17.98

#test shapiro-wilk
shapiro.test(b_2020)
shapiro.test(b_2019)

#creación de gráficos de densidad
ggplot(matrix2019, aes(x = b_2019, y = after_stat(density)),color=) +
  geom_histogram(bins=15,fill="gray") +
  xlim(14, 20.5) +
  geom_vline(xintercept=Endural2019, color="blue")+
  geom_vline(xintercept=Aldura2019, color="red")+
  theme_bw() 

ggplot(matrix2020, aes(x = b_2020, y = after_stat(density))) +
  geom_histogram(bins=15,fill="gray") +
  xlim(13.5, 19)+
  geom_vline(xintercept=Endural2020, color="blue")+
  geom_vline(xintercept=Aldura2020, color="red")+
  theme_bw() 


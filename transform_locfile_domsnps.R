# Autor: Joaquin Algorta Bove
# transform_locfile_domsnps
# Funcion: transforma los marcadores tipo snps con problemas de homeologia y con una buena segregacion
#          a formato loc (requerido para JoinMap). Se codifican como dominantes
# Archivo de entrada: marcadores snps con problemas de homeologia y polimórficos entre Endural y Aldura -> Archivo "dom_snps_buenos.csv"
#         Filas- marcadores moleculares
#         Columnas- genotipo individuos  

library(tidyverse)

#cargar los datos
df <- read.csv("C:/RStudio/resultados/EnduralXAldura/snps/dom_snps_buenos.csv" ,header = TRUE, row.names = 1)

#numero de las columnas de endural y aldura
endural <- "X908623031001_A_1"
aldura <- "X908623031002_F_2"

col_endural <- as.numeric(which( colnames(df)== endural))
col_aldura <-as.numeric(which( colnames(df)== aldura ))

df_transformed<-df


##CAMBIOS PREVIOS - eliminar los 2s
#si ademas de 2s tenemos 0s entonces cambiamos todos los 2s como 1s
#si ademas de 2s tenemos 1s entonces cmabiamos los 2s comom 1s y los 1s
#como 0s para hacer que los 2s sean los dominantes

for (i in 1:nrow(df_transformed)){
  if (df_transformed[i,col_endural]==2) {
    if (df_transformed[i,col_aldura]==0){
      df_transformed[i,][df_transformed[i,] == 2] <- 1
    }
    else if (df_transformed[i,col_aldura]==1) {
      df_transformed[i,][df_transformed[i,] == 1] <- 0
      df_transformed[i,][df_transformed[i,] == 2] <- 1
    }
  }
  else if (df_transformed[i,col_aldura]==2) {
    if (df_transformed[i,col_endural]==0){
      df_transformed[i,][df_transformed[i,] == 2] <- 1
    }
    else if (df_transformed[i,col_endural]==1) {
      df_transformed[i,][df_transformed[i,] == 1] <- 0
      df_transformed[i,][df_transformed[i,] == 2] <- 1
    }
  }
}

#parental 1 - poner a si es un 0 y d si es un 1
#ir individuo a individuo poniendo a o d si coinciden 
for (i in 1:nrow(df_transformed)){
  if (df_transformed[i,col_endural]==1){
    df_transformed[i,col_endural]="d"
    for (j in 3:ncol(df_transformed)){
      if (df_transformed[i,j]==1 & !is.na(df_transformed[i,j])){df_transformed[i,j]="d" }
      else if (is.na(df_transformed[i,j])){next}
    }
  }
  else if (df_transformed[i,col_endural]==0){
    df_transformed[i,col_endural]="a"
    for (j in 3:ncol(df_transformed)){
      if (df_transformed[i,j]==0  & !is.na(df_transformed[i,j])){df_transformed[i,j]="a" }
      else if (is.na(df_transformed[i,j])){next}
    }
  }
}

### Transformar formato locfile 

#parental 2 - poner b si es un 0 y c si es un 1
#ir individuo a individuo poniendo a o d si coinciden 
for (i in 1:nrow(df_transformed)){
  if (df_transformed[i,col_aldura]==1){
    df_transformed[i,col_aldura]="c"
    for (j in 3:ncol(df_transformed)){
      if (df_transformed[i,j]==1 & !is.na(df_transformed[i,j])){df_transformed[i,j]="c" }
      else if (is.na(df_transformed[i,j])){next}
    }
  }
  else if (df_transformed[i,col_aldura]==0){
    df_transformed[i,col_aldura]="b"
    for (j in 3:ncol(df_transformed)){
      if (df_transformed[i,j]==0  & !is.na(df_transformed[i,j])){df_transformed[i,j]="b" }
      else if (is.na(df_transformed[i,j])){next}
    }
  }
}


#este archivo ya se puede meter como data set en JoinMap
#hay que introducir la información sobre cuantos marcadores hay y cuantas lineas se han analizado
write.csv(df_transformed, file="C:/RStudio/resultados/EnduralXAldura/snps/dom_snps_locfile.csv")


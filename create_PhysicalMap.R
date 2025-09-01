# Autor: Joaquin Algorta Bove
# create_PhysicalMap
# Funcion: a partir de un mapa genético busca cada marcador en el archivo original
#          y añade su posición fisica. Despues transforma esa posicion a medidas 
#          proporcionales a las del mapa genético para poder comparar a escala
# Archivo de entrada: mapa genetico de un cromosoma y archivo original 
#                     el nombre de marcador debe coincidir
 
library("readxl")
library(tidyverse)

darts_and_snps <- read_excel("C:/RStudio/resultados/EnduralXAldura/darts_and_snps/darts_and_snps_locfile.xlsx", sheet = "correspondencias")
mapa<-read_excel("C:/RStudio/resultados/EnduralXAldura/darts_and_snps/mapas.xlsx",sheet="7B")

#añadir la posicion fisica en una nueva columna
result <- mapa %>%
  add_column(PosFisica = "")


for (i in 1:nrow(result)){
  for (j in 1:nrow(darts_and_snps)){
    if (result[i,3]==darts_and_snps[j,3]){
      result[i,7]<- as.character(darts_and_snps[j,4])
    }
  }
}
#añadir la posicion fisica proporcional a la posicion genetica
result2 <- result %>%
  add_column(PosMapa = 0)

result2$PosFisica<-as.numeric(result2$PosFisica)
for (i in 1:nrow(result2)){
  pos<-as.numeric(result2[i,7])/as.numeric(max(result2$PosFisica))*as.numeric(result2[nrow(result2),6])
  result2[i,8]<-pos
}

write.csv(result2,file="C:/RStudio/resultados/EnduralXAldura/darts_and_snps/mapa7BconPos.csv",row.names=FALSE)


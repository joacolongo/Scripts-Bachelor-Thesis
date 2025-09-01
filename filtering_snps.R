# Autor: Joaquin Algorta Bove
# filtering_snps
# Funcion: filtra marcadores tipo dart y los separa en 3 grupos:
#          1- marcadores con genotipo conocido 0/1 y 1/0 en Endural/Aldura. Archivo salida: "codom_snps.csv"
#          2- marcadores con genotipo 0/2 o 1/2 o 2/0 o 2/1 en Endural/Aldura. Archivo salida: "dom_snps.csv"
#          3- marcadores con genotipo faltante en algun parental. A estos se aplica filtro MAF>0.3 y NA<10%
#          y una imputación del genotipo parental desconocido con el genotipo conocido del otro parental.
#         Archivo salida: "NAparents_imputed_snps.csv"
# Archivo de entrada: datos brutos de genotipado de la poblacion tipo SNP 
#         Filas- marcadores moleculares
#         Columnas- genotipo individuos + datos adicionales

library(tidyverse)

#cargo datos brutos
snps <- read.csv("C:/RStudio/datos/EnduralXAldura/Report_SW23-428_SNP_mapping.csv", na = "-")
names(snps) <- snps[5,]
snps<- snps[6:nrow(snps),]

#### filtro diferentes entre parentales sin tener en cuenta los NA
endural <- "908623031001_A_1"
aldura <- "908623031002_F_2"

#pongo las columnas de endural y aldura juntas para facilitar la visualizacion
snps <- snps %>%
  select(c(1:25,all_of(endural),all_of(aldura), everything()))

#numero de las columnas de endural y aldura
col_endural <- as.numeric(which( colnames(snps)== endural))
col_aldura <-as.numeric(which( colnames(snps)== aldura ))

######## GRUPO 1 - 0/1 y 1/0
#filtramos solo los que tengan un 0 en un parental y un 1 en el otro
snps_group1 <-snps %>%
  filter((.[[col_endural]]  == 0 & .[[col_aldura]]  == 1) | (.[[col_endural]]  == 1 & .[[col_aldura]]  == 0) )

write.csv(snps_group1,file="C:/RStudio/resultados/EnduralXAldura/snps/codom_snps.csv",row.names=FALSE)


######## GRUPO 2 - marcadores que sean 0/2 o 1/2 o 2/0 o 2/1
#filtramos solo los que tengan un 2 en un parental y un 1 o 0 en el otro
snps_group2 <-snps %>%
  filter((.[[col_endural]]  == 0 & .[[col_aldura]]  == 2) | (.[[col_endural]]  == 1 & .[[col_aldura]]  == 2) 
         | (.[[col_endural]]  == 2 & .[[col_aldura]]  == 0) | (.[[col_endural]]  == 2 & .[[col_aldura]]  == 1))

write.csv(snps_group2,file="C:/RStudio/resultados/EnduralXAldura/snps/dom_snps.csv",row.names=FALSE)

######## GRUPO 3 - marcadores que sean NA/algo

snps_group3 <-snps %>%
  filter((is.na(.[[col_endural]]) &  !is.na(.[[col_aldura]])) | (is.na(.[[col_aldura]]) &  !is.na(.[[col_endural]])))

#filtramos pNA<0.1
tsnps <- as.data.frame(t(snps_group3))#hay que trasponer
names(tsnps) <- tsnps[1,]

#nos quedamos solo con lo que sean datos de genotipado
tsnps<-tsnps[col_endural:nrow(tsnps),]
for (i in 1:ncol(tsnps)){
  tsnps[,i]<-as.numeric(tsnps[,i])
}
K = 0
prueba <- tsnps
col <- NULL
for (i in c(1:ncol(prueba))){
  num_na <-(sum(is.na(prueba[,i])))/109 #aquí va el número de variedades
  if (num_na <= 0.2){
    K = K + 1
    col <- c(col, i)
  }
}
prueba_final <- subset(prueba, select = col)
dim(prueba_final)

#filtramos por maf>0.3
skim100 <- prueba_final
dim(skim100)
resultado <- data.frame()
vec = NULL
col <- NULL
for (i in c(1:ncol(skim100))){
  datos <- as.data.frame(table(skim100[,i]))
  for(j in c(1:nrow(datos))){
    datos$Freq[j]<- (datos$Freq[j]/109)#poner el numero de variedades
  }
  if (datos[1,2] >= 0.3  && datos[2,2] >= 0.3){#poner la frecuencia por al que se filtra
    col <- c(col, i)
  }
  resultado[i,1] <- datos[1,2]
  resultado[i,2] <- datos[2,2]
  resultado[i,3] <- (1 - (datos[1,2] + datos[2,2]))
  rowname <- colnames(skim100[i])
  vec = c(vec,rowname)
}
skim100_final <- subset(skim100, select = col)
dim(skim100_final)

#guardamos
NAparents_snps<-as.data.frame(t(skim100_final))

write.csv(NAparents_snps,file="C:/RStudio/resultados/EnduralXAldura/snps/NAparents_snps.csv",row.names=TRUE)

#eliminar los que tienen 2/Na en los parentales porque no podemos imputar si es un 1 o un 0

for (i in 1:nrow(NAparents_snps)){
  if (is.na(NAparents_snps[i,1]) & NAparents_snps[i,2]==2){NAparents_snps<-NAparents_snps[-i,]}
  else if (is.na(NAparents_snps[i,2]) & NAparents_snps[i,1]==2){NAparents_snps<-NAparents_snps[-i,]}
}
#lo hago dos veces porque me dejaba algun valor
for (i in 1:nrow(NAparents_snps)){
  if (is.na(NAparents_snps[i,1]) & NAparents_snps[i,2]==2){NAparents_snps<-NAparents_snps[-i,]}
  else if (is.na(NAparents_snps[i,2]) & NAparents_snps[i,1]==2){NAparents_snps<-NAparents_snps[-i,]}
}


#imputar parental restante 
for (i in 1:nrow(NAparents_snps)){
  if (is.na(NAparents_snps[i,1]) & NAparents_snps[i,2]==1){NAparents_snps[i,1]<-0}
  else if (is.na(NAparents_snps[i,1]) & NAparents_snps[i,2]==0){NAparents_snps[i,1]<-1}
  else if (is.na(NAparents_snps[i,2]) & NAparents_snps[i,1]==1){NAparents_snps[i,2]<-0}
  else if (is.na(NAparents_snps[i,2]) & NAparents_snps[i,1]==0){NAparents_snps[i,2]<-1}
}

write.csv(NAparents_snps,file="C:/RStudio/resultados/EnduralXAldura/snps/NAparents_imputed_snps.csv",row.names=TRUE)
                
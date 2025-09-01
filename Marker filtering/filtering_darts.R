# Autor: Joaquin Algorta Bove
# filtering_darts
# Funcion: filtra aquellos marcadores dart polimorficos entre Endural y Aldura 
#          y recupera marcadores con NA en un parental, MAF>0.3 y NA<10% por medio de una
#          imputación del genotipo parental desconocido con el genotipo conocido del otro parental.
# Archivo de entrada: datos brutos de genotipado de la poblacion tipo DART 
#         Filas- marcadores moleculares
#         Columnas- genotipo individuos + datos adicionales

library(tidyverse)

#cargo datos brutos
darts <- read.csv("C:/RStudio/datos/EnduralXAldura/Report_SW23-428_SilicoDArT.csv")
names(darts) <- darts[5,]
darts<- darts[6:nrow(darts),]
darts<-darts[,c(1,25:ncol(darts))]

endural <- "908623031001_A_1"
aldura <- "908623031002_F_2"

#pongo las columnas de endural y aldura juntas para facilitar la visualizacion
darts <- darts %>%
  select(c(CloneID,all_of(endural),all_of(aldura), everything()))

#numero de las columnas de endural y aldura
col_endural <- as.numeric(which( colnames(darts)== endural))
col_aldura <-as.numeric(which( colnames(darts)== aldura ))


#### filtramos los marcadores diferentes entre parentales, sin tener en cuenta los NA
darts_poly <- darts %>%
  filter(.[[col_endural]] != "-" & .[[col_aldura]] != "-" & .[[col_endural]] != .[[col_aldura]])

#contamos cuantos NA hay en Endural y Aldura como info
datos_endural <- as.data.frame(table(darts[,col_endural]))
datos_aldura<-as.data.frame(table(darts[,col_aldura]))
NA_endural<-as.numeric(datos_endural[1,2])
NA_aldura<-as.numeric(datos_aldura[1,2])

#### filtramos los que tienen un NA y un numero 0 o 1
NA_darts <- darts %>%
  filter(.[[col_endural]] != .[[col_aldura]] & .[[col_endural]] == "-" |.[[col_endural]] != .[[col_aldura]] & .[[col_aldura]] == "-")

#aplicamos filtro datos faltanes<10%
tNA_darts <- as.data.frame(t(NA_darts))#hay que trasponer
names(tNA_darts) <- tNA_darts[1,]
tNA_darts<-tNA_darts[2:nrow(tNA_darts),]

#pasarlo a numeric para cambiar - por NA
for(i in 1:ncol(tNA_darts)){
  tNA_darts[,i]<-as.numeric(tNA_darts[,i]) 
}

K = 0
prueba <- tNA_darts
col <- NULL
for (i in c(1:ncol(prueba))){
  num_na <-(sum(is.na(prueba[,i])))/109 #aquí va el número de variedades
  if (num_na <= 0.1){
    K = K + 1
    col <- c(col, i)
  }
}
prueba_final <- subset(prueba, select = col)
dim(prueba_final)

#aplicamos filtro MAF>0.3
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
  #resultado[i,1] <- datos[1,2]
  #resultado[i,2] <- datos[2,2]
  #resultado[i,3] <- (1 - (datos[1,2] + datos[2,2]))
  #rowname <- colnames(skim100[i])
  #vec = c(vec,rowname)
}
skim100_final <- subset(skim100, select = col)
dim(skim100_final)
NA_darts_poly<-skim100_final

#Imputar valor 1 o 0 como el contrario al otro parental

for (i in 1:ncol(NA_darts_poly)){
  if (is.na(NA_darts_poly[1,i])){
    if (NA_darts_poly[2,i]==0){
      NA_darts_poly[1,i]<-1
    }
    else {
      NA_darts_poly[1,i]<-0
    }
  }
  if (is.na(NA_darts_poly[2,i])){
    if (NA_darts_poly[1,i]==0){
      NA_darts_poly[2,i]<-1
    }
    else {
      NA_darts_poly[2,i]<-0
    }
  }
}

#añadir estos valores al anterior dataframe
tNA_darts_poly<-as.data.frame(t(NA_darts_poly))

write.csv(tNA_darts_poly, "C:/RStudio/resultados/EnduralXAldura/darts/selected_darts.csv")
tNA_darts_poly <- read.csv("C:/RStudio/resultados/EnduralXAldura/darts/selected_darts.csv")

names(tNA_darts_poly)<-names(darts_poly)
tNA_darts_poly[1,]<-as.character(tNA_darts_poly[1,])
selected_darts <- full_join(darts_poly, tNA_darts_poly)
dim(selected_darts)

#pasarlo a numeric para cambiar - por NA
for(i in 1:ncol(selected_darts)){
  selected_darts[,i]<-as.numeric(selected_darts[,i]) 
}

#archivo definitivo para pasar a formato reconocido por JoinMap
write.csv(selected_darts, "C:/RStudio/resultados/EnduralXAldura/darts/selected_darts.csv",row.names = FALSE)

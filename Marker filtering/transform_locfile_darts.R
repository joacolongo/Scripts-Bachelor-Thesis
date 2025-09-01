# Autor: Joaquin Algorta Bove
# transform_locfile_darts
# Funcion: transforma los marcadores tipo darts a formato loc (requerido para JoinMap)
#          Además filtra marcadores de calidad MAF>0.05 datos faltantes<20%
#          Por ultimo, añade posiciones fisicas conocidas a los marcadores con ese dato
#          a partir del archivo original de darts
# Archivo de entrada: marcadores polimórficos tipo darts entre Endural y Aldura  
#         Filas- marcadores moleculares
#         Columnas- genotipo individuos  

library(tidyverse)

# cargar los datos
df <- read.csv("C:/RStudio/resultados/EnduralXAldura/darts/selected_darts.csv" ,header = TRUE, row.names = 1)

#numero de las columnas de endural y aldura
endural <- "X908623031001_A_1"
aldura <- "X908623031002_F_2"

col_endural <- as.numeric(which( colnames(df)== endural))
col_aldura <-as.numeric(which( colnames(df)== aldura ))

df_transformed<-df

#parental 1 - poner a si es un 0 y d si es un 1
#ir individuo a individuo poniendo a o d si coinciden 
for (i in 1:nrow(df_transformed)){
  if (df_transformed[i,col_endural]==1 ){
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
write.csv(df_transformed, file="C:/RStudio/resultados/EnduralXAldura/darts/darts_locfile.csv")

#### FILTRAMOS  MARCADORES DE BUENA CALIDAD
#eliminamos los parentales
darts<-df_transformed[,3:ncol(df_transformed)]

#Filtramos monomorficos
remove<-as.numeric()
for (i in 1:nrow(darts)){
  row<-as.character(darts[i,])
  as<-length(which(row=="a"))
  bs<-length(which(row=="b"))
  cs<-length(which(row=="c"))
  ds<-length(which(row=="d"))
  nas<-sum(is.na(row))
  
  if(as+nas==107 | bs+nas==107 | cs+nas==107| ds+nas==107){
    remove<-c(remove,as.numeric(i))
  }
}
darts<-darts[-remove,]
dim(darts)

#filtro por datos faltantes menor del 20% 
tdarts<-as.data.frame(t(darts))
K = 0
prueba <- tdarts
col <- NULL
for (i in c(1:ncol(prueba))){
  num_na <-(sum(is.na(prueba[,i])))/107 #aquí va el número de variedades
  if (num_na <= 0.2){
    K = K + 1
    col <- c(col, i)
  }
}
prueba_final <- subset(prueba, select = col)
dim(prueba_final)

#filtro MAF>0.05 (minimum allele frequency)
skim100 <- prueba_final
dim(skim100)
resultado <- data.frame()
vec = NULL
col <- NULL
for (i in c(1:ncol(skim100))){
  datos <- as.data.frame(table(skim100[,i]))
  for(j in c(1:nrow(datos))){
    datos$Freq[j]<- (datos$Freq[j]/107)#poner el numero de variedades
  }
  rowname <- colnames(skim100[i])
  if (datos[1,2] >= 0.05  && datos[2,2] >= 0.05){#poner la frecuencia por al que se filtra
    col <- c(col, i)
  }
  #resultado[i,1] <- datos[1,2]
  #resultado[i,2] <- datos[2,2]
  #resultado[i,3] <- (1 - (datos[1,2] + datos[2,2]))
  
  vec = c(vec,rowname)
}
skim100_final <- subset(skim100, select = col)
dim(skim100_final)

darts_FINAL<- as.data.frame(t(skim100_final))

#tramsformar NA a "-"
darts_FINAL[is.na(darts_FINAL)] <- "-"

#añadir el nombre del cromosoma al que pertenece el marcador
#para ello hay que leer el archivo codom_darts.csv
original_darts<-read.csv("C:/RStudio/datos/EnduralXAldura/Report_SW23-428_SilicoDArT.csv")
names(original_darts) <- original_darts[5,]

darts <- darts_FINAL %>%
  add_column(chrom = "")

col_chrom <-as.numeric(which( colnames(original_darts)== "Chrom_Wheat_Durum_Svevo_v2" ))
row_numbers<-which(  original_darts$CloneID %in% row.names(darts_FINAL))
for (i in 1:nrow(darts)){
  marker<-rownames(darts)[i]
  rowNUMBER<-which( original_darts$CloneID %in% marker)
  #print(rowNUMBER)
  #print(marker)
  darts[i,108]<-original_darts[(rowNUMBER),col_chrom][1]
}

#añadir la posicion física de cada marcador
darts <- darts %>%
  add_column(pos = "")

col_pos <-as.numeric(which( colnames(original_darts)== "ChromPos_Wheat_Durum_Svevo_v2" ))
row_numbers<-which(  original_darts$CloneID %in% row.names(darts_FINAL))
for (i in 1:nrow(darts)){
  marker<-rownames(darts)[i]
  rowNUMBER<-which(  original_darts$CloneID %in% marker)
  #print(rowNUMBER)
  #print(marker)
  darts[i,109]<-original_darts[(rowNUMBER),col_pos][1]
}

write.csv(darts,file="C:/RStudio/resultados/EnduralXAldura/darts/darts_locfile.csv",row.names=TRUE)

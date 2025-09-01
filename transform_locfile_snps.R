# Autor: Joaquin Algorta Bove
# transform_locfile_snps
# Funcion: transforma los marcadores tipo snps a formato loc (requerido para JoinMap)
#          Además filtra marcadores de calidad MAF>0.05 datos faltantes<20%
#          Por ultimo, añade posiciones fisicas conocidas a los marcadores con ese dato
#          a partir del archivo original de snps
# Archivo de entrada: marcadores tipo snps polimórficos entre Endural y Aldura  
#         Filas- marcadores moleculares
#         Columnas- genotipo individuos  

library(tidyverse)

# cargar los snps filtrados
df <- read.csv("C:/RStudio/resultados/EnduralXAldura/snps/codom_JUNTOS_snps.csv",row.names = 1)

#filtro diferentes entre parentales dejando las que tienen NAs
endural <- "X908623031001_A_1"
aldura <- "X908623031002_F_2"

col_endural <- as.numeric(which( colnames(df)== endural))
col_aldura <-as.numeric(which( colnames(df)== aldura ))

#transformar los 2s de los parentales por 0 o 1 dependiendo del que no esta en el otro parental
df_adapted<-df
for (i in 1:nrow(df_adapted)){
  parenta<-df_adapted[i,col_endural]
  parentb<-df_adapted[i,col_aldura]
  if (parenta==2){
    if (parentb==0){df_adapted[i,col_endural]=1}
    else if (parentb==1){df_adapted[i,col_endural]=0}
  }
  else if (parentb==2){
    if (parenta==0){df_adapted[i,col_aldura]=1}
    else if (parenta==1){df_adapted[i,col_aldura]=0}
  }
}

#Transformar formato locfile

#poner los 2 como h
df_adapted[df_adapted == 2] <- "h"

#poner a si igual al primer parental y b si igual al segundo
for (i in 1:nrow(df_adapted)){
  parenta=df_adapted[i,col_endural]
  parentb=df_adapted[i,col_aldura]
  for (j in (col_aldura+1):ncol(df_adapted)){
    if (df_adapted[i,j]==parenta & !is.na(df_adapted[i,j]) & df_adapted[i,j]!="h"){df_adapted[i,j]="a"}
    else if (df_adapted[i,j]==parentb & !is.na(df_adapted[i,j]) & df_adapted[i,j]!="h"){df_adapted[i,j]="b"}
  }
}

#columna parental 1 todos los valores igual a a y parental 2 a b
df_adapted[col_endural,]<-as.character(df_adapted[col_endural,])
df_adapted[col_aldura,]<-as.character(df_adapted[col_aldura,])
df_adapted <- df_adapted %>%
  mutate_at(col_endural,list(~str_replace(., ., "a")))
df_adapted <- df_adapted %>%
  mutate_at(col_aldura,list(~str_replace(., ., "b")))

#este archivo ya se puede meter como data set en JoinMap
#hay que introducir la información sobre cuantos marcadores hay y cuantas lineas se han analizado

#### FILTRAMOS  MARCADORES DE BUENA CALIDAD
#eliminamos los parentales
snps<-df_adapted[,3:ncol(df_adapted)]

#eliminar monomorficos
remove<-as.numeric()
for (i in 1:nrow(snps)){
  row<-as.character(snps[i,])
  unos<-length(which(row=="a"))
  ceros<-length(which(row=="b"))
  
  if(unos==0 | ceros==0){
    remove<-c(remove,as.numeric(i))
  }
}

snps<-snps[-remove,]


#filtro pNA<0.2
tsnps<-as.data.frame(t(snps))

K = 0
prueba <- tsnps
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

#filtro maf>0.05
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
  if (datos[1,2] >= 0.05  && datos[2,2] >= 0.05){#poner la frecuencia por al que se filtra
    col <- c(col, i)
  }
  #resultado[i,1] <- datos[1,2]
  #resultado[i,2] <- datos[2,2]
  #resultado[i,3] <- (1 - (datos[1,2] + datos[2,2]))
  rowname <- colnames(skim100[i])
  vec = c(vec,rowname)
}
skim100_final <- subset(skim100, select = col)
dim(skim100_final)

snps_FINAL<- as.data.frame(t(skim100_final))

#tramsformar NA a "-"
snps_FINAL[is.na(snps_FINAL)] <- "-"

#añadir el nombre del cromosoma al que pertenece el marcador
#para ello hay que leer el archivo codom_snps.csv
original_snps<-read.csv("C:/RStudio/datos/EnduralXAldura/Report_SW23-428_SNP_mapping.csv")
names(original_snps) <- original_snps[5,]

snps <- snps_FINAL %>%
  add_column(chrom = "")

col_chrom <-as.numeric(which( colnames(original_snps)== "Chrom_Wheat_Durum_Svevo_v2" ))
row_numbers<-which(  original_snps$AlleleID %in% row.names(snps_FINAL))
for (i in 1:nrow(snps)){
  marker<-rownames(snps)[i]
  rowNUMBER<-which(  original_snps$AlleleID %in% marker)
  #print(rowNUMBER)
  #print(marker)
  snps[i,108]<-original_snps[(rowNUMBER),col_chrom]
}

#añadir la posicion fisica en una nueva columna 
snps <- snps %>%
  add_column(pos = "")

col_pos <-as.numeric(which( colnames(original_snps)== "ChromPos_Wheat_Durum_Svevo_v2" ))
row_numbers<-which(  original_snps$AlleleID %in% row.names(snps_FINAL))
for (i in 1:nrow(snps)){
  marker<-rownames(snps)[i]
  rowNUMBER<-which(  original_snps$AlleleID %in% marker)
  #print(rowNUMBER)
  #print(marker)
  snps[i,109]<-original_snps[(rowNUMBER),col_pos]
}



write.csv(snps,file="C:/RStudio/resultados/EnduralXAldura/snps/codom_snps_locfile.csv",row.names=TRUE)

           
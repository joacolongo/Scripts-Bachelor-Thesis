# Autor: Joaquin Algorta Bove
# pre_analysis_darts
# Funcion: Analiza el numero de marcadores tipo dart polimorficos entre Endural y Aldura 
#          en cada cromosoma en posiciones geneticas y fisicas. Para ello crea archivos 
#          con histogramas y gráficos de densidad por cromosomas. Ademas, analiza distribución 
#          de marcadores calculando medias maximos y minimosde las distancias entre marcadores 
#          adyacentes y crea boxplots
# Archivo de entrada: archivo bruto de DARTs sin las 5 lineas primera que no aportan inf 
#         Filas- marcadores moleculares
#         Columnas- informacion de marcadores + genotipo individuos 

library(dplyr)
library(RColorBrewer)
library(writexl)
library(ggplot2)
library(xlsx)
library(readxl)
library(stringr)

setwd("C:/RStudio/")
# load / cleaned data
df <- read_xlsx("resultados/darts/datosLimpios1.xlsx")

#Especificar lineas parentales
genot1 <- 662 #Endural
genot2 <- 672 #Aldura

#SUBSET 
#creamos subset con las columnas que usaremos + genotipos parentales
df_subset <- df[,c("CloneID","AlleleSequence",  "Chrom_Wheat_Durum_Svevo_v2", 
                   "ChromPos_Wheat_Durum_Svevo_v2", "Chrom_Wheat_ConsensusMap_version_4", 
                   "position_ConsensusMap_version4", genot1, genot2)]
df_subset<-as.data.frame(df_subset)

#FILTER DATA
#filtramos marcadores polimorficos entre Endural y Aldura
marcadores_df <- df_subset %>%
  filter(.[[7]] != "-" & .[[8]] != "-" & .[[7]] != .[[8]])

#ponemos todos los marcadores pretenecientes al genomio D en un mismo grupo
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "1D", "D")))
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "2D", "D")))
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "3D", "D")))
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "4D", "D")))
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "5D", "D")))
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "6D", "D")))
marcadores_df <- marcadores_df %>%
  mutate_at("Chrom_Wheat_ConsensusMap_version_4",list(~str_replace(., "7D", "D")))
  
#Guardamos datos filtrados
write.xlsx(marcadores_df, file=paste("resultados/darts/datosLimpios",genot1, genot2 ,".xlsx", sep = "_"),col.names=TRUE)

#### Analisis
#Nº marcadores
num_marc <- marcadores_df %>% count(CloneID) %>% summarise(numero_marcadores = n())

#Nº marcadores por cromosoma
marc_CR <- marcadores_df %>% arrange(Chrom_Wheat_Durum_Svevo_v2) %>% count(Chrom_Wheat_Durum_Svevo_v2) 

marc_genet_CR <- marcadores_df %>% arrange(Chrom_Wheat_ConsensusMap_version_4) %>% count(Chrom_Wheat_ConsensusMap_version_4) 

#Guardamos en diferentes hojas del excel
write.xlsx(num_marc, file = paste("resultados/darts/analisis_marcadores", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "num_marcadores", col.names = TRUE, row.names = FALSE, append = FALSE)
write.xlsx(marc_genet_CR, file = paste("resultados/darts/analisis_marcadores", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "genet_mark", col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx(marc_CR, file = paste("resultados/darts/analisis_marcadores", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "CR_marc", col.names = TRUE, row.names = FALSE, append = TRUE)

#### Histogramas posicion fisica
#dataframe con columnas de Cromosoma y posicion fisica solo
data_histF <- marcadores_df[,c("Chrom_Wheat_Durum_Svevo_v2", "ChromPos_Wheat_Durum_Svevo_v2")]
data_histF <- na.omit(data_histF)

#crear histograma para cada cromosoma-> function(i)
data_histF$ChromPos_Wheat_Durum_Svevo_v2 <- as.numeric(data_histF$ChromPos_Wheat_Durum_Svevo_v2)
hist_F = lapply(sort(unique(data_histF$Chrom_Wheat_Durum_Svevo_v2)), function(i) {
  ggplot(data_histF[data_histF$Chrom_Wheat_Durum_Svevo_v2==i,], aes(x = ChromPos_Wheat_Durum_Svevo_v2, fill=Chrom_Wheat_Durum_Svevo_v2)) +
    geom_histogram(show.legend=FALSE, colour="grey20", fill = "#26783B", bins = 50) + 
    facet_wrap(~Chrom_Wheat_Durum_Svevo_v2, scales ="free") #separate by chr
})
hist_F


#density plot
hist_F2 = lapply(sort(unique(data_histF$Chrom_Wheat_Durum_Svevo_v2)), function(i) {
  ggplot(data_histF[data_histF$Chrom_Wheat_Durum_Svevo_v2==i,], aes(x = ChromPos_Wheat_Durum_Svevo_v2, fill=Chrom_Wheat_Durum_Svevo_v2)) +
    geom_density(show.legend=FALSE, fill = "#408D61", alpha = 1) +
    facet_wrap(~Chrom_Wheat_Durum_Svevo_v2, scales ="free") 
})
hist_F2

#guardamos histogramas en pdf
pdf(file = paste("resultados/darts/histograma_fisic", genot1, genot2 ,".pdf", sep = "_"))
print(hist_F) 
dev.off()

#save density plot
pdf(file = paste("resultados/darts/density_fisic", genot1, genot2 ,".pdf", sep = "_"))
print(hist_F2)     
dev.off()

#### Histogramas posicion genetica
#dataframe con cromosoma y posicion genetica
data_histG <- marcadores_df[,c("Chrom_Wheat_ConsensusMap_version_4", "position_ConsensusMap_version4")]
data_histG <- na.omit(data_histG)
data_histG$position_ConsensusMap_version4 <- as.numeric(data_histG$position_ConsensusMap_version4)


hist_G = lapply(sort(unique(data_histG$Chrom_Wheat_ConsensusMap_version_4)), function(i) {
  ggplot(data_histG[data_histG$Chrom_Wheat_ConsensusMap_version_4==i,], aes(x = position_ConsensusMap_version4, fill=Chrom_Wheat_ConsensusMap_version_4)) +
    geom_histogram(show.legend=FALSE, colour="grey20", fill = "#009999", bins = 50) +
    facet_wrap(~Chrom_Wheat_ConsensusMap_version_4, scales ="free") 
})

#guarmos en pdf
pdf(file = paste("resultados/darts/histograma_genetic", genot1, genot2 ,".pdf", sep = "_"))
print(hist_G)     
dev.off()

#density plot genet
dens_G2 = lapply(sort(unique(data_histG$Chrom_Wheat_ConsensusMap_version_4)), function(i) {
  ggplot(data_histG[data_histG$Chrom_Wheat_ConsensusMap_version_4==i,], aes(x = position_ConsensusMap_version4, fill=Chrom_Wheat_ConsensusMap_version_4)) +
    geom_density(show.legend=FALSE, fill = "#009999", alpha = 1) +
    facet_wrap(~Chrom_Wheat_ConsensusMap_version_4, scales ="free") 
})

#guaramos density plot
pdf(file = paste("resultados/darts/density_genet", genot1, genot2 ,".pdf", sep = "_"))
print(dens_G2)     
dev.off()

#### ANALISIS DISTRIBUCIÓN - POSICION FISICA

#ordenamos los marcadores por cromosoma en orden ascendente de posicion fisica
#calculamos distancias entre marcadores adyacentes (columna distance)
distance.F <- data_histF %>% 
  group_by(Chrom_Wheat_Durum_Svevo_v2) %>%
  arrange(Chrom_Wheat_Durum_Svevo_v2,ChromPos_Wheat_Durum_Svevo_v2) %>%
  mutate(difference = ChromPos_Wheat_Durum_Svevo_v2 - lag(ChromPos_Wheat_Durum_Svevo_v2, 1)) %>%
  as.data.frame()

     

write.xlsx(distance.F, file= paste("resultados/darts/distancias_F", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "distancias", col.names = TRUE, row.names = FALSE, append = FALSE)

#MEDIA, MIN & MAX de las distancias
sum_dist_F <- distance.F %>%                               
  group_by(Chrom_Wheat_Durum_Svevo_v2) %>% 
  mutate(average = mean(difference, na.rm = TRUE),
         min = min(difference, na.rm = TRUE),
         max = max(difference, na.rm = TRUE)) %>%
  select(Chrom_Wheat_Durum_Svevo_v2, average, min, max) %>%
  group_by(Chrom_Wheat_Durum_Svevo_v2, average, min, max) %>% summarise() %>%
  as.data.frame()
write.xlsx(sum_dist_F, file = paste("resultados/darts/distancias_F", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "calculos", col.names = TRUE, row.names = FALSE, append = TRUE)


#BOXPLOT 
png(file = paste("resultados/darts/boxplot_F", genot1, genot2,".png", sep = "_"),width = 1000, height = 500)
boxplot(distance.F$difference ~ distance.F$Chrom_Wheat_Durum_Svevo_v2, outline = FALSE, xlab = "CR", ylab = "eje y",
                      col = "#26783B", frame = FALSE)
dev.off()

#boxplot con outliers
png(file = paste("resultados/darts/PUNTOS-BOXPLOT-F", genot1, genot2 ,".png", sep = "_"), width = 1000, height = 500)
boxplot_PUNTOSF <- ggplot(distance.F, aes(x = Chrom_Wheat_Durum_Svevo_v2, y = difference, fill = Chrom_Wheat_Durum_Svevo_v2)) +
  geom_boxplot(show.legend = FALSE)
print(boxplot_PUNTOSF)
dev.off()

####GENETIC DISTANCIES
#mismo procedimiento
distancia.G <- data_histG %>%
  group_by(Chrom_Wheat_ConsensusMap_version_4) %>%
  arrange(Chrom_Wheat_ConsensusMap_version_4,position_ConsensusMap_version4) %>%
  mutate(differenceG = position_ConsensusMap_version4 - lag(position_ConsensusMap_version4, 1)) %>%
  as.data.frame()
write.xlsx(distancia.G, file = paste("resultados/darts/distancias_G", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "dist", col.names = TRUE, row.names = FALSE, append = FALSE)

#MEDIA, MIN & MAX OF DISTANCES
sum_dist_G <- distancia.G %>%                               
  group_by(Chrom_Wheat_ConsensusMap_version_4) %>% 
  mutate(average = mean(differenceG, na.rm = TRUE),
         min = min(differenceG, na.rm = TRUE),
         max = max(differenceG, na.rm = TRUE)) %>%
  select(Chrom_Wheat_ConsensusMap_version_4, average, min, max) %>%
  group_by(Chrom_Wheat_ConsensusMap_version_4, average, min, max) %>% summarise() %>%
  as.data.frame()
write.xlsx(sum_dist_G, file = paste("resultados/darts/distancias_G", genot1, genot2 ,".xlsx", sep = "_"), sheetName = "calculos", col.names = TRUE, row.names = FALSE, append = TRUE)

#BOXPLOT GENETIC

png(file = paste("resultados/darts//boxplot_G", genot1, genot2,".png", sep = "_"),width = 1000, height = 500)
boxplot(distancia.G$differenceG ~ distancia.G$Chrom_Wheat_ConsensusMap_version_4, outline = FALSE, xlab = "CR", ylab = "eje y",
        col = "#009999", frame = FALSE)
dev.off()

#boxplot con outliers
png(file = paste("resultados/darts/PUNTOS-BOXPLOT-G", genot1, genot2 ,".png", sep = "_"), width = 1000, height = 500)
boxplot_PUNTOSG <- ggplot(distancia.G, aes(x = Chrom_Wheat_ConsensusMap_version_4, y = differenceG, fill = Chrom_Wheat_ConsensusMap_version_4)) +
  geom_boxplot(show.legend = FALSE)
print(boxplot_PUNTOSG)
dev.off()






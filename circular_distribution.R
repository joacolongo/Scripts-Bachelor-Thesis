# Autor: Joaquin Algorta Bove
# circular_distribution
# Funcion: crea un gráfico circular que muestra el número de marcadores según posiciones
#          físicas y genéticas. Además, separa por sectores los distintos cromosomas
# Archivo de entrada: marcadores polimórficos entre Endural y Aldura  
#         Filas- marcadores moleculares
#         Columnas- datos de marcadores + genotipo individuos  

library(tidyverse)
library(readxl)
library(svglite)
setwd("C:/RStudio/")

#######   POR POSICIONES FÍSICAS

#cargar datos limpios y seleccionar datos de position
df<- read_xlsx("resultados/pre_estudio/snps/datosLimpios_662_672_.xlsx")
mapping <- df[,c("AlleleID",  "Chrom_Wheat_Durum_Svevo_v2", 
                   "ChromPos_Wheat_Durum_Svevo_v2")]

#cambiar los nombres de las columnas para que sean iguales al otro data frame
colnames(mapping) <- c("AlleleID","CHROM","POS")

#pasamos la col POS a numeric
mapping <-mapping %>% mutate_at("POS", as.numeric)

#ordenar de menor a mayor
mapping <- mapping %>%
  group_by(CHROM) %>%
  arrange(POS, .by_group=TRUE)

#crear dataframe con posición inicial y final de los cromosomas
chromo_maxpos =
  mapping %>%
  select(CHROM,POS) %>%
  group_by(CHROM) %>%
  summarise(max=max(POS))
chromo_maxpos$Chr = as.numeric(chromo_maxpos$CHROM)
chromo_maxpos = na.omit(chromo_maxpos)

delta = 0
chromo_maxpos$init = c(-delta,cumsum(chromo_maxpos$max))[1:14]+delta
chromo_maxpos$fin = c(cumsum(chromo_maxpos$max))

#quitar los chrUn
mapping.chr = filter(mapping, CHROM != "chrUn") %>%
  arrange(CHROM,POS)

#posiciones todas continuas a pesar de cambiar de cromosoma
for (i in 1:nrow(mapping.chr)){
  if(mapping.chr[i,]$CHROM==2){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[1,5])
  }
  else if(mapping.chr[i,]$CHROM==3){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[2,5])
  }
  else if(mapping.chr[i,]$CHROM==4){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[3,5])
  }
  else if(mapping.chr[i,]$CHROM==5){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[4,5])
  }
  else if(mapping.chr[i,]$CHROM==6){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[5,5])
  }
  else if(mapping.chr[i,]$CHROM==7){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[6,5])
  }
  else if(mapping.chr[i,]$CHROM==8){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[7,5])
  }
  else if(mapping.chr[i,]$CHROM==9){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[8,5])
  }
  else if(mapping.chr[i,]$CHROM==10){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[9,5])
  }
  else if(mapping.chr[i,]$CHROM==11){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[10,5])
  }
  else if(mapping.chr[i,]$CHROM==12){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[11,5])
  }
  else if(mapping.chr[i,]$CHROM==13){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[12,5])
  }
  else if(mapping.chr[i,]$CHROM==14){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[13,5])
  }
}
new.pos<-mapping.chr$POS


#parte específica de la densidad de marcadores, 
#histograma donde el tamaño del intervalo es mbp
#puse este valor para que haya 700 intervalos
mbp = 14166613
#redondea las posiciones para que se identifiquen con su intervalo
window = data.frame(new.pos = mbp*new.pos%/%(mbp))
#contar cuantos marcadores hay en cada intervalo
window = window %>% count(new.pos)

window$new.pos <- as.numeric(window$new.pos)
for (i in 1:700){
  valor=mbp*i
  if (!(valor %in% window$new.pos)){
    window = window %>% add_row(new.pos=valor,n=0)
  }
}

#añadimos columna con 0 para dibujar la linea basal
window <- window %>%
  add_column(min_val=0)

#nombre cromosomas
cromosomas=c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B")



#circulo1<-ggplot()+
ggplot()+
  geom_rect(data=chromo_maxpos,
            aes(xmin=init,xmax=fin,alpha=as.factor(as.numeric(CHROM)%%2),fill=as.factor(as.numeric(CHROM)%%2),
            ymin=log(1),ymax=5.5))+
  
  theme_classic()+
  
  scale_alpha_manual(values = c(0.65,0.6))+
  scale_fill_manual(values = alpha(c("#B1E4E1", "#639FAB")))+
  
  coord_polar()+
 
  guides(alpha="none")+
  theme(legend.position = "none")+   # comentar esto para leyenda
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank())+
  xlab("")+ylab("")+
  geom_text(data=chromo_maxpos,
           aes(label=cromosomas,x=(fin+init)/2,y=6.2))+
  geom_line(data=window,
            aes(x=new.pos,y=3.15+(n/50)),alpha=0.9,color="#65334D")+
  geom_line(data=window, aes(x=new.pos, y=3.15+(min_val/50)), linewidth = 0.2, alpha=0.7,color="black")

ggsave("resultados/pre_estudio/snps/circulos/circulofisico_snps.svg",plot=circulo1, width = 20, height = 20, units = "cm",device="svg")


##### POR POSICIONES GENETICAS

mapping <- df[,c("AlleleID",  "Chrom_Wheat_ConsensusMap_version_4", 
                      "position_ConsensusMap_version4")]

#cambiar los nombres de las columnas para que sean iguales al otro data frame
colnames(mapping) <- c("ID","CHROM","POS")

#pasamos la col CHROM a numeric
mapping <- mapping %>%
  mutate_at("CHROM",as.numeric)
mapping$CHROM<-as.numeric(mapping$CHROM) #los D te los cambia a NA
#pasamos la col POS a numeric
mapping <-mapping %>% mutate_at("POS", as.numeric)

#ordenar de menor a mayor
mapping <- mapping %>%
  group_by(CHROM) %>%
  arrange(POS, .by_group=TRUE)
mapping = na.omit(mapping)

chromo_maxpos =
  mapping %>%
  select(CHROM,POS) %>%
  group_by(CHROM) %>%
  summarise(max=max(POS))
chromo_maxpos$Chr = as.numeric(chromo_maxpos$CHROM)
chromo_maxpos = na.omit(chromo_maxpos)

delta = 0
chromo_maxpos$init = c(-delta,cumsum(chromo_maxpos$max))[1:14]+delta
chromo_maxpos$fin = c(cumsum(chromo_maxpos$max))

#quitar los chrUn y NA
mapping.chr = filter(mapping, CHROM != "chrUn") %>%
  arrange(CHROM,POS)

#posiciones todas continuas a pesar de cambiar de cromosoma
for (i in 1:nrow(mapping.chr)){
  if(mapping.chr[i,]$CHROM==2){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[1,5])
  }
  else if(mapping.chr[i,]$CHROM==3){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[2,5])
  }
  else if(mapping.chr[i,]$CHROM==4){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[3,5])
  }
  else if(mapping.chr[i,]$CHROM==5){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[4,5])
  }
  else if(mapping.chr[i,]$CHROM==6){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[5,5])
  }
  else if(mapping.chr[i,]$CHROM==7){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[6,5])
  }
  else if(mapping.chr[i,]$CHROM==8){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[7,5])
  }
  else if(mapping.chr[i,]$CHROM==9){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[8,5])
  }
  else if(mapping.chr[i,]$CHROM==10){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[9,5])
  }
  else if(mapping.chr[i,]$CHROM==11){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[10,5])
  }
  else if(mapping.chr[i,]$CHROM==12){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[11,5])
  }
  else if(mapping.chr[i,]$CHROM==13){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[12,5])
  }
  else if(mapping.chr[i,]$CHROM==14){
    mapping.chr[i,]$POS=as.numeric(mapping.chr[i,]$POS)+as.numeric(chromo_maxpos[13,5])
  }
}

new.pos<-mapping.chr$POS

#parte específica de la densidad de marcadores, 
#histograma donde el tamaño del intervalo es mbp
#puse este valor para que haya 700 intervalos
mbp = 2.9
#redondea las posiciones para que se identifiquen con su intervalo
window = data.frame(new.pos = mbp*new.pos%/%(mbp))
# contar cuantos marcadores hay en cada intervalo
window = window %>% count(new.pos)

window$new.pos <- as.numeric(window$new.pos)
for (i in 1:700){
  valor=mbp*i
  if (!(valor %in% window$new.pos)){
    window = window %>% add_row(new.pos=valor,n=0)
  }
}

window <- window %>%
  add_column(min_val=0)

cromosomas=c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B")


#circulo1<-ggplot()+
ggplot()+
  geom_rect(data=chromo_maxpos,
            aes(xmin=init,xmax=fin,alpha=as.factor(as.numeric(CHROM)%%2),fill=as.factor(as.numeric(CHROM)%%2),
                ymin=log(1),ymax=7.5))+
  
  theme_classic()+
  
  scale_alpha_manual(values = c(0.65,0.6))+
  scale_fill_manual(values = alpha(c( "#C1CDCD","#E0EEEE")))+
  
  coord_polar()+
  
  guides(alpha="none")+
  theme(legend.position = "none")+   
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank())+
  xlab("")+ylab("")+
  geom_text(data=chromo_maxpos,
            aes(label=cromosomas,x=(fin+init)/2,y=9))+
  geom_line(data=window,
            aes(x=new.pos,y=3.15+(n/50)),alpha=5,color="black",linewidth=0.65)+
  geom_line(data=window, aes(x=new.pos, y=3.15+(min_val/50)), linewidth = 0.2, alpha=0.7,color="red")

#guardamos
ggsave("resultados/pre_estudio/snps/circulos/circulogenetico_snps.svg",plot=circulo1, width = 20, height = 20, units = "cm",device="svg")


# 2	Códigos de R 
# 2.1	Calcular promedio a partir de las repeticiones
# Programa que agrupa todos los experimentos con iguales parámetros 
# y calcula el promedio de las variables calculadas
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())
require(plyr)
require(binom)
require(psych)
source("binom_sd.R")
source("beta_sd.R")

file_in<-"resultados_por_seed"

seed_data<-read.table(paste(file_in,".txt",sep=""),header=TRUE)
seed_data<-unique(seed_data)

variables_iguales<-names(seed_data)[c(2:19)]

data <- ddply(seed_data, variables_iguales, summarise,
              nMonteCarlo = sum(!is.na(percolado)),
              firstSeed = min(seed),
              percolado_temp = mean(percolado),
              numeroCNTPercolando_temp = mean(numeroCNTPercolando),
              numeroContactosCNTPercolando_temp = mean(numeroContactosCNTPercolando),
              PromedioNumeroContactosCNTPercolando_temp = mean(PromedioNumeroContactosCNTPercolando),
              porcentajeContactos_temp = mean(porcentajeContactos),
              numeroNodosA_temp = mean(numeroNodosA),
              numeroElementosA_temp = mean(numeroElementosA),
              promedioConexionesA_temp = mean(promedioConexionesA),
              R_temp = geometric.mean(R)
)
colnames(data)[21]<-"percolado"
colnames(data)[22]<-"numeroCNTPercolando"
colnames(data)[23]<-"numeroContactosCNTPercolando"
colnames(data)[24]<-"PromedioNumeroContactosCNTPercolando"
colnames(data)[25]<-"porcentajeContactos"
colnames(data)[26]<-"numeroNodosA"
colnames(data)[27]<-"numeroElementosA" 
colnames(data)[28]<-"promedioConexionesA" 
colnames(data)[29]<-"R"

file_name<-"resumen_seed.txt"
write.table(data, file_name,row.names = FALSE, sep=" ")


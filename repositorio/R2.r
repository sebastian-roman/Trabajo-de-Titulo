# 2.2	Determinar número de experimentos para cada conjunto de parámetros
# Programa que calcula las repiticiones realizadas para cada 
#  conjunto de parámetros. Se utiliza para hacer seguimiento
#  de los experimentos realizados y los faltantes
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
require(xlsx)
rm(list = ls())

file_in<-"resumen_seed"
data<-read.table(paste(file_in,".txt",sep=""),header=TRUE)
data<-unique(data)

variables<-c('n','pAglomerado','strain','nAglomerados')
n<-levels(factor(data$n))
pAglomerado<-levels(factor(data$pAglomerado))
strain<-levels(factor(data$strain))
nAglomerados<-levels(factor(data$nAglomerados))
N<-levels(factor(data$N))

lista<-list(n,pAglomerado,strain,nAglomerados)
names(lista)<-variables
fors<-expand.grid(lista,stringsAsFactors = TRUE)
fors<- fors[order(fors$nAglomerados,fors$n,fors$pAglomerado), ]

data_variables<-data[,variables]

id<-apply(fors, 1, function(x) {
    iguales<-apply(data_variables, 1, function(y) { all(x == y) }) 
    sapply(N,function(z){sum(data$nMonteCarlo[iguales & data$N == z])})
})
id<-aperm(id)

for (i in 1:length(N)){
    var_name<-paste("N",N[i],sep="_")
    fors[var_name]<-id[,i]
}

write.xlsx(x = fors, file = "combinaciones_for.xlsx", sheetName = "todos",row.names = FALSE)


# 2.3	Agrupar resultados según variables a graficar
# Programa que agrupa las variables segun un conjunto de parametros 
# indicados. Permite agrupar todas los conjuntos que tengan la misma 
# variable, a ser utilizada en el eje X de los gráficos y por otra 
# parte los divide según las variable que determina las distintas
# curvas a graficar
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())
file_in<-"resumen_seed"

ejeX<-'pAglomerado'
curvas<-'N'
var_in<-c('strain','n')

# ejeX<-'N'
# curvas<-'strain'
# var_in<-c('n','pAglomerado')

# ejeX<-'strain'
# curvas<-'N'
# var_in<-c('pAglomerado','n')

data<-read.table(paste(file_in,".txt",sep=""),header=TRUE)
data<-unique(data)
data<-data[order(data[[curvas]]),]

numero_var_in = length(var_in)

lista_valores_var<-list()
cantidad_var<-vector()

for (k in 1:numero_var_in){
    valores_var<-levels(factor(data[,var_in[k]]))    
    lista_valores_var<-append(lista_valores_var, list(valores_var))
    cantidad_var[k]<-length(valores_var)    
}
comb<-expand.grid(lista_valores_var,stringsAsFactors = FALSE)
comb <- comb[order(comb$Var1), ]

if (numero_var_in == 1)
    {comb <- as.data.frame(lista_valores_var)}

filas<-length(comb[,1])
columnas<-length(comb[1,])


for (i in 1:filas)
{
    var <- as.vector(as.matrix(comb[i, ]))
    var_data<-data
    for (j in 1:columnas)
    {
        var_data <- var_data[var_data[,var_in[j]] == var[j],]
    }
    
    var<-gsub("\\.", "_", var)
    
    var_name<-""
    for (j in 1:columnas-1)
    {
        var_name<-paste(var_name,var_in[j],var[j],"_",sep="")

    }
    var_name<-paste(var_name,var_in[columnas],var[columnas],sep = "")
    var_name<-sub("_", "", var_name)
    
    #Se escribe el archivo .txt
    file_name<-paste("ejeX_",ejeX,"_curvas_",curvas,".txt",sep = "",collapse="")
    
    #la primera vez se crea de nuevo el archivo
    if (i==1) 
        {cat(var_name,"=[\n",file=file_name,sep="")}
    else
        {cat(var_name,"=[\n",file=file_name,sep="",append = TRUE)}
    
    write.table(var_data, file = file_name,append = TRUE,row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("];",file=file_name,sep="\n",append=TRUE)

    
    #Se escribe los archivos resumen.m y run.m
    file_name_m<-paste("ejeX_",ejeX,"_curvas_",curvas,".m",sep = "",collapse="")
    file_name_m_run<-paste("run_ejeX_",ejeX,"_curvas_",curvas,".m",sep = "",collapse="")
    
    #la primera vez se crea de nuevo el archivo
    if (i==1) 
    {
        #archivo run.m
        cat("clear all\n",file=file_name_m_run,sep="")
        file_name_aux<-paste("ejeX_",ejeX,"_curvas_",curvas,sep = "",collapse="")
        cat("run('",file_name_aux,"');\n",file=file_name_m_run,sep="",append = TRUE)
        cat("guardar = 0;\n",file=file_name_m_run,sep="",append = TRUE)
        n_ejeX <- length(levels(factor(data[[ejeX]])))
        cat("n_ejeX = ",n_ejeX,";\n",file=file_name_m_run,sep="",append = TRUE)
        cat(ejeX," = 0;\n",file=file_name_m_run,sep="",append = TRUE)
        cat(curvas," = 0;\n",file=file_name_m_run,sep="",append = TRUE)
        cat("delay = 3;\n",file=file_name_m_run,sep="",append = TRUE)
        # cat("graficar_phi_c = 0;\n",file=file_name_m_run,sep="",append = TRUE)
        # cat("guardar_phi_c = 0;\n\n",file=file_name_m_run,sep="",append = TRUE)
        cat("if(exist('",var_name,"','var') && sum(size(",var_name,"))>0)\n",
file=file_name_m_run,sep="",append = TRUE)
        cat("graficar_ejeX_curvas(",var_name,",",ejeX,",",curvas,",n_ejeX,guardar);pause(delay);\n",
file=file_name_m_run,sep="",append = TRUE)
        cat("end\n",file=file_name_m_run,sep="",append = TRUE)
        #archivo resumen.m
        cat(var_name,"=[\n",file=file_name_m,sep="")
    }
    else
    {
        #archivo run.m
        cat("if(exist('",var_name,"','var') && sum(size(",var_name,"))>0)\n",
file=file_name_m_run,sep="",append = TRUE)
        cat("graficar_ejeX_curvas(",var_name,",",ejeX,",",curvas,",n_ejeX,guardar);pause(delay);\n",
file=file_name_m_run,sep="",append = TRUE)
        cat("end\n",file=file_name_m_run,sep="",append = TRUE)
        #archivo resumen.m
        cat(var_name,"=[\n",file=file_name_m,sep="",append = TRUE)
    }
    
    write.table(var_data, file = file_name_m,append = TRUE,row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("];",file=file_name_m,sep="\n",append=TRUE)
    
} #end for i 

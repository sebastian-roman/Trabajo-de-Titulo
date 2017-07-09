% 3.2	Graficar probabilidad de percolación y resistencia para todas las combinaciones de parámetros
% Funcion que grafica la probabilidad de percolación, numero de CNT, numero de contactos, numero promedio de 
% contactos, resistencia, y cambio en la resistencia, para todas las combinaciones posibles de parámetros, dados
% la variable que se grafica en el eje X y la variable que determina las distintas curvas

clear all
run('ejeX_strain_curvas_N');
guardar = 1;
n_ejeX = 5;
strain = 0;
N = 0;
delay = 3;

if(exist('pAglomerado0_n1','var') && sum(size(pAglomerado0_n1))>0)
graficar_ejeX_curvas(pAglomerado0_n1,strain,N,n_ejeX,guardar);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_n5','var') && sum(size(pAglomerado0_n5))>0)
graficar_ejeX_curvas(pAglomerado0_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_n9','var') && sum(size(pAglomerado0_n9))>0)
graficar_ejeX_curvas(pAglomerado0_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_1_n1','var') && sum(size(pAglomerado0_1_n1))>0)
graficar_ejeX_curvas(pAglomerado0_1_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_1_n5','var') && sum(size(pAglomerado0_1_n5))>0)
graficar_ejeX_curvas(pAglomerado0_1_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_1_n9','var') && sum(size(pAglomerado0_1_n9))>0)
graficar_ejeX_curvas(pAglomerado0_1_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_2_n1','var') && sum(size(pAglomerado0_2_n1))>0)
graficar_ejeX_curvas(pAglomerado0_2_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_2_n5','var') && sum(size(pAglomerado0_2_n5))>0)
graficar_ejeX_curvas(pAglomerado0_2_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_2_n9','var') && sum(size(pAglomerado0_2_n9))>0)
graficar_ejeX_curvas(pAglomerado0_2_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_3_n1','var') && sum(size(pAglomerado0_3_n1))>0)
graficar_ejeX_curvas(pAglomerado0_3_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_3_n5','var') && sum(size(pAglomerado0_3_n5))>0)
graficar_ejeX_curvas(pAglomerado0_3_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_3_n9','var') && sum(size(pAglomerado0_3_n9))>0)
graficar_ejeX_curvas(pAglomerado0_3_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_4_n1','var') && sum(size(pAglomerado0_4_n1))>0)
graficar_ejeX_curvas(pAglomerado0_4_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_4_n5','var') && sum(size(pAglomerado0_4_n5))>0)
graficar_ejeX_curvas(pAglomerado0_4_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_4_n9','var') && sum(size(pAglomerado0_4_n9))>0)
graficar_ejeX_curvas(pAglomerado0_4_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_6_n1','var') && sum(size(pAglomerado0_6_n1))>0)
graficar_ejeX_curvas(pAglomerado0_6_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_6_n5','var') && sum(size(pAglomerado0_6_n5))>0)
graficar_ejeX_curvas(pAglomerado0_6_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_6_n9','var') && sum(size(pAglomerado0_6_n9))>0)
graficar_ejeX_curvas(pAglomerado0_6_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_8_n1','var') && sum(size(pAglomerado0_8_n1))>0)
graficar_ejeX_curvas(pAglomerado0_8_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_8_n5','var') && sum(size(pAglomerado0_8_n5))>0)
graficar_ejeX_curvas(pAglomerado0_8_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado0_8_n9','var') && sum(size(pAglomerado0_8_n9))>0)
graficar_ejeX_curvas(pAglomerado0_8_n9,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado1_n1','var') && sum(size(pAglomerado1_n1))>0)
graficar_ejeX_curvas(pAglomerado1_n1,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado1_n5','var') && sum(size(pAglomerado1_n5))>0)
graficar_ejeX_curvas(pAglomerado1_n5,strain,N,n_ejeX,guardar);pause(delay);
end
if(exist('pAglomerado1_n9','var') && sum(size(pAglomerado1_n9))>0)
graficar_ejeX_curvas(pAglomerado1_n9,strain,N,n_ejeX,guardar);pause(delay);
end
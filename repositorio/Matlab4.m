% 3.4	 Graficar el cambio en la percolación al aplicar deformación, para todas las aglomeraciones y geometrías 
% Funcion que grafica el cambio en la pendiente de la porobabilidad de percolacion al variar la deformación 
% para cada concentracion utilizada y para todas las aglomeraciones y geometrias de nanotubos utilizadas
clear all
run('ejeX_strain_curvas_N');
guardar = 1;
n_ejeX = 5;
strain = 0;
N = 0;
delay = 1;

graficar_pendiente_ejeX_curvas(pAglomerado0_n1,strain,N,n_ejeX,0,0,1);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_1_n1,strain,N,n_ejeX,0,1,2);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_2_n1,strain,N,n_ejeX,0,1,3);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_3_n1,strain,N,n_ejeX,0,1,4);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_4_n1,strain,N,n_ejeX,0,1,5);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_6_n1,strain,N,n_ejeX,0,1,6);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_8_n1,strain,N,n_ejeX,0,1,7);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado1_n1,strain,N,n_ejeX,guardar,1,8);pause(delay);

graficar_pendiente_ejeX_curvas(pAglomerado0_n5,strain,N,n_ejeX,0,0,1);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_1_n5,strain,N,n_ejeX,0,1,2);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_2_n5,strain,N,n_ejeX,0,1,3);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_3_n5,strain,N,n_ejeX,0,1,4);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_4_n5,strain,N,n_ejeX,0,1,5);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_6_n5,strain,N,n_ejeX,0,1,6);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_8_n5,strain,N,n_ejeX,0,1,7);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado1_n5,strain,N,n_ejeX,guardar,1,8);pause(delay);

graficar_pendiente_ejeX_curvas(pAglomerado0_n9,strain,N,n_ejeX,0,0,1);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_1_n9,strain,N,n_ejeX,0,1,2);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_2_n9,strain,N,n_ejeX,0,1,3);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_3_n9,strain,N,n_ejeX,0,1,4);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_4_n9,strain,N,n_ejeX,0,1,5);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_6_n9,strain,N,n_ejeX,0,1,6);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado0_8_n9,strain,N,n_ejeX,0,1,7);pause(delay);
graficar_pendiente_ejeX_curvas(pAglomerado1_n9,strain,N,n_ejeX,guardar,1,8);pause(delay);

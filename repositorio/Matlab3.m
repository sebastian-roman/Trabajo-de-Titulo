% 3.3	Graficar el cambio en la percolación al aplicar deformación, para aglomeración y geometría fija
% Funcion que grafica el cambio en la pendiente de la porobabilidad de percolacion al variar la deformación 
% para cada concentracion utilizada y una aglomeración y geometría de nanotubo fija
function [] = graficar_pendiente_ejeX_curvas(r,ejeX,curvas,n_ejeX,guardar,mantener,color_i)

numero_curvas = size(r,1)/n_ejeX;
X_name = inputname(2);
curvas_name = inputname(3);
vars = {'N','L','n','n_aglomerados','n_aglomerados_x','n_aglomerados_y','n_aglomerados_z','pAglomerado',...
        'theta_max','d0','dvdw','dt','Lx','Ly','Lz','strain','poisson','R_pol','N_MonteCarlo','first_seed',...
        'percolado','numero_cnt_percolando','numero_contactos_cnt_percolando','promedio_numero_contactos_cnt_percolando','porcentaje_contactos',...
        'numeroNodosA','numeroElementosA','promedioConexionesA','R'};

for i = 1:size(r,2)
    eval([vars{i} '=  r(:,i);']);
end

texto1 ={[ ],['    L=' num2str(min(L)) '      n=' num2str(min(n)) '      N=' num2str(min(N)) ],[ ],...
    ['    p_a_g= ' num2str(min(pAglomerado),'%1.1f')   '       n_a_g= ' ...
 num2str(min(n_aglomerados)) ' / [' num2str(min(n_aglomerados_x)) ',' num2str(min(n_aglomerados_y)) ',' num2str(min(n_aglomerados_z)) ']' ],[ ],...
['    \epsilon=' num2str(min(strain)) '      \nu=' num2str(min(poisson)) '      \theta_m_a_x=' num2str(min(theta_max),3) ],[ ],...
    ['    L_x/L=' num2str(min(Lx./L)) '    L_y/L=' num2str(min(Ly./L)) '    L_z/L=' num2str(min(Lz./L))   ],[ ],...
    ['    d_0= ' num2str(min(d0)) '     d_v_d_w= ' num2str(min(dvdw)) '     d_t= ' num2str(min(dt))],[ ],...
    ['    R_p_o_l= ' num2str(min(R_pol)) '     N_m_c= ' num2str(min(N_MonteCarlo))  ]};

if (mantener == 0)
    h = figure();
    set(gcf, 'Position', get(0,'Screensize'));
    annotation('textbox', [.14 .53 .20 .43], 'String',texto1,'FontSize',13);
end

for (i = 1:numero_curvas)
%     valor = eval(['num2str(' curvas_name '(i*n_ejeX));']);
    valor = eval([ curvas_name '(i*n_ejeX);']);
    if strcmp(curvas_name,'N')
        v_cnt = pi*(min(d0)^2/4)*min(L);
        v_vol = (min(Lx)*min(Ly)*min(Lz));
        valor = valor.*(v_cnt/v_vol)*100;
    end
    valor = num2str(valor);
    if (length(curvas_name)==11) %pAglomerado
        labels{i}=[ 'p_a_g = ', valor ];
    elseif (length(curvas_name)==6) %strain
        labels{i}=[ '\epsilon = ', valor ];
    elseif strcmp(curvas_name,'N') %N se pasa a phi
        labels{i}=[ '\phi_v_o_l = ', valor(1:4) ];
    else % n
        labels{i}=[ curvas_name ' = ', valor ];
    end

end

vectors = {[0 1], [0 0.5], [0 1]};
colores = combvec(vectors{:})';

eval(['X_plot = ' X_name '(((i-1)*n_ejeX+1):(i*n_ejeX));']);
if strcmp(X_name,'N')
    v_cnt = pi*(min(d0)^2/4)*min(L);
    v_vol = (min(Lx)*min(Ly)*min(Lz));
    X_plot = X_plot.*v_cnt./v_vol*100;
    X_name = '\phi_v_o_l';
end

for i = 1:numero_curvas
    eval(['curvas_plot(i) = ' curvas_name '(i*n_ejeX);']);
end    
if strcmp(curvas_name,'N')
    v_cnt = pi*(min(d0)^2/4)*min(L);
    v_vol = (min(Lx)*min(Ly)*min(Lz));
    curvas_plot = curvas_plot.*v_cnt./v_vol*100;
    curvas_name = '\phi_v_o_l';
end

for i = 1:numero_curvas
%     subplot(2,3,2); hold on;
    percolado_plot = percolado(((i-1)*n_ejeX+1):(i*n_ejeX));
    promedio_pendiente_percolado(i) = 0;
    suma_pendiente_percolado = 0;
    for j = 1:(n_ejeX-1)
        suma_pendiente_percolado = suma_pendiente_percolado + (percolado_plot(j+1) - percolado_plot(j));
    end
    promedio_pendiente_percolado(i) = -suma_pendiente_percolado;
end
    subplot(2,2,2); 
    hold on;
    plot(curvas_plot,promedio_pendiente_percolado,'color',colores(color_i,1:3));
    
    if (color_i == 8)
        title(['Cambio en la percolacion          n = ' num2str(min(n))]);xlabel(curvas_name);
     ylabel('\Delta probabilidad percolacion');
        label = {'p_a_g=0','p_a_g=0.1','p_a_g=0.2','p_a_g=0.3','p_a_g=0.4','p_a_g=0.6','p_a_g=0.8','p_a_g=1.0',};
        legend(label,'Fontsize',8,'Location','NorthWest');grid on;
    end

if (guardar)
    if strcmp(X_name,'\phi_v_o_l')
        X_name = 'phi';
    end
    if strcmp(curvas_name,'N')
        curvas_name = 'phi';
    end
    if strcmp(curvas_name,'\phi_v_o_l')
        curvas_name = 'phi';
    end

    rname = ['cambio_pendiente_ejeX_' X_name '_curvas_' curvas_name '_parametros_' inputname(1)];
    nombre_archivo = [pwd '\cambio_percolacion\' rname];
    f=gcf; %f is the handle of the figure you want to export
    figpos=getpixelposition(f); %dont need to change anything here
    resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
    set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]);
    
    saveas(f,nombre_archivo ,'fig');
    saveas(f,nombre_archivo ,'png');
    print('-dpdf','-r1200',nombre_archivo);
 
end %if guardar

end %function
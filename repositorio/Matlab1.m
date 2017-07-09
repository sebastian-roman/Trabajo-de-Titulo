% 3	Códigos MATLAB 
% 3.1	Graficar probabilidad de percolación y resistencia para un conjunto de parámetros fijos
% Funcion que grafica la probabilidad de percolación, numero de CNT, numero de contactos, numero promedio de 
% contactos, resistencia, y cambio en la resistencia, para un experimento con iguales parámetros.
% Se entrega la variable que se grafica en el eje X y la variable que determina las distintas curvas
function [] = graficar_ejeX_curvas(r,ejeX,curvas,n_ejeX,guardar)

numero_curvas = size(r,1)/n_ejeX;
X_name = inputname(2);
curvas_name = inputname(3);
vars = {'N','L','n','n_aglomerados','n_aglomerados_x','n_aglomerados_y','n_aglomerados_z','pAglomerado',...
'theta_max','d0','dvdw','dt','Lx','Ly','Lz','strain','poisson','R_pol','N_MonteCarlo','first_seed',...       'percolado','numero_cnt_percolando','numero_contactos_cnt_percolando',
'promedio_numero_contactos_cnt_percolando','porcentaje_contactos',...
        'numeroNodosA','numeroElementosA','promedioConexionesA','R'};

for i = 1:size(r,2)
    eval([vars{i} '=  r(:,i);']);
end

texto1 ={[ ],['    L=' num2str(min(L)) '      n=' num2str(min(n)) '      N=' num2str(min(N)) ],[ ],...
    ['    p_a_g= ' num2str(min(pAglomerado),'%1.1f')   '       n_a_g= ' ...
 num2str(min(n_aglomerados)) ' / [' num2str(min(n_aglomerados_x)) ',' num2str(min(n_aglomerados_y)) ','       num2str(min(n_aglomerados_z)) ']' ],[ ],...
['    \epsilon=' num2str(min(strain)) '      \nu=' num2str(min(poisson)) '      \theta_m_a_x=' num2str(min(theta_max),3) ],[ ],...
    ['    L_x/L=' num2str(min(Lx./L)) '    L_y/L=' num2str(min(Ly./L)) '    L_z/L=' num2str(min(Lz./L))   ],[ ],...
    ['    d_0= ' num2str(min(d0)) '     d_v_d_w= ' num2str(min(dvdw)) '     d_t= ' num2str(min(dt))],[ ],...
    ['    R_p_o_l= ' num2str(min(R_pol)) '     N_m_c= ' num2str(min(N_MonteCarlo))  ]};

h = figure();
set(gcf, 'Position', get(0,'Screensize'));
annotation('textbox', [.14 .53 .20 .43], 'String',texto1,'FontSize',13);

for (i = 1:numero_curvas)

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
labels

vectors = {[0 1], [0 0.5], [0 1]};
colores = combvec(vectors{:})';

eval(['X_plot = ' X_name '(((i-1)*n_ejeX+1):(i*n_ejeX));']);
if strcmp(X_name,'N')
    v_cnt = pi*(min(d0)^2/4)*min(L);
    v_vol = (min(Lx)*min(Ly)*min(Lz));
    X_plot = X_plot.*v_cnt./v_vol*100;
    X_name = '\phi_v_o_l [%]';
end

for i = 1:numero_curvas
    subplot(2,3,2); hold on;
    percolado_plot = percolado(((i-1)*n_ejeX+1):(i*n_ejeX));
    plot(X_plot,percolado_plot,'color',colores(i,1:3));
    if (i == numero_curvas)
        title('Porcentaje percolacion');xlabel(X_name);ylabel('% perc');
     legend(labels,'Fontsize',8,'Location','NorthEast');grid on;
    end
end

for i = 1:numero_curvas
    subplot(2,3,3); hold on; 
    plot(X_plot,numero_cnt_percolando(((i-1)*n_ejeX+1):(i*n_ejeX)),'color',colores(i,1:3));
    if (i == numero_curvas)
        title('Numero CNT percolando');
     xlabel(X_name);ylabel('# CNT percolando');legend(labels,'Fontsize',8,'Location','Best');grid on;
    end
end

for i = 1:numero_curvas
    subplot(2,3,4); hold on;
    promedio_numero_contactos_cnt_percolando_plot=
promedio_numero_contactos_cnt_percolando(((i-1)*n_ejeX+1):(i*n_ejeX));
    plot(X_plot,promedio_numero_contactos_cnt_percolando_plot,'color',colores(i,1:3));
    if (i == numero_curvas)
        title('Promedio Numero Contactos CNT percolando');
        xlabel(X_name);
        ylabel('promedio contactos CNT percolando');legend(labels,'Fontsize',8,'Location','Best');grid on;
    end
end

for i = 1:numero_curvas
    subplot(2,3,5); hold on;
    plot(X_plot,R(((i-1)*n_ejeX+1):(i*n_ejeX)),'color',colores(i,1:3));
    if (i == numero_curvas)
        title('Resistencia');xlabel(X_name);
     ylabel('R [\Omega]');legend(labels,'Fontsize',8,'Location','SouthEast');grid ON; set(gca,'YScale','log');
    end
end

for i = 1:numero_curvas
    subplot(2,3,6); hold on;
    R_plot = R(((i-1)*n_ejeX+1):(i*n_ejeX));
    RR0_plot = (R_plot-R_plot(1))./R_plot(1);
    plot(X_plot,RR0_plot,'color',colores(i,1:3));
    if (i == numero_curvas)
    	title('Cambio en la Resistencia');
xlabel(X_name);ylabel('\DeltaR');legend(labels,'Fontsize',8,'Location','NorthEast');grid ON; 
    end
end


if (guardar)
    if strcmp(X_name,'\phi_v_o_l')
        X_name = 'phi';
    end
    if strcmp(curvas_name,'N')
        curvas_name = 'phi';
    end

    rname = ['ejeX_' X_name '_curvas_' curvas_name '_parametros_' inputname(1)];
    nombre_archivo = [pwd '\graficos\' rname];
    f=gcf; %f is the handle of the figure you want to export
    figpos=getpixelposition(f); %dont need to change anything here
    resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
    set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]);
    
    saveas(h,nombre_archivo ,'fig');
    print(h,'-dpdf','-r1200',nombre_archivo);
    % print(h,'-djpeg','-r600',nombre_archivo);
end %if guardar

end %function

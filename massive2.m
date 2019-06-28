%algorithm to generate list of cluster images to be sent to cluster
%detection
% select only relevants
%------------------------------------
% corregir 21 de febrero
% los datos deben ser tomados de variable factora no segmented, ver linea
% linea 250 CorreccionsCh

%MARZO 2019 SE ELIMINA LO SIGUIENTE
%L=logical(erodedBWch1);
%box=regionprops(L,'Area', 'BoundingBox', 'centroid','Perimeter'); 
%areasas = cat(1, box.Area);
%Centroidc = cat(1, box.Centroid);
%bbox = cat(1, box.BoundingBox);
%Pericc = cat(1, box.Perimeter);
%bboxint=int16(bbox);


cuentaclusters=0;
clear dG;
%nana=imcomplement(erodedBW);
%nana=erodedBWch1;  
%nana1=label2rgb(nana*512); %imagen que muestra a colores clusters
logics=0; % bandera que indica que el cluster actual no contiene celulas
egray=rgb2gray(he);
testi=uint8(Itrue).*egray;
%Ccu = bwconncomp(erodedBWch1);
%statsh5=regionprops(Ccu, 'Area');
%areash5 = cat(1, statsh5.Area);

%numPixels = cellfun(@numel,Ccu.PixelIdxList);

%while (min(numPixels)< Th)
%    [biggest,idx] = min(numPixels);
%    erodedBWch1(Ccu.PixelIdxList{idx}) = 0;
%    Ccu = bwconncomp(erodedBWch1);
%    numPixels = cellfun(@numel,Ccu.PixelIdxList);
%end

%centroidf=0;

%dic2018
%figure
%hold on
%imshow(he)
vectordetected=0;  % este vector contiene el listado de indices validos de bbox, del total de elementos
clear yfitC;
clear nam;
flag1=0;
clear crHv;
clear crSv;
clear crVv;
% lee la dimension de tabla anterior
%if c1==1 % no hubo ninguna celula antes
%     c1=1;
%   else  %lee tabla de salida existente
%     c1=c1-1; % Tab_clust se genero para colocar celulas del background durante back def process
     
%end

     
%conteo matriz debe iniciar desde donde se quedo external cells
if (flagelo==0)
    conteo_matrix=1;
else
    conteo_matrix=height(Tab_clust(:,1));
end
if (back_density==0)
  Areafin=0;
  Eccenfin=0;
  Centroidfin=0;
  Mayor_ellipsefin=0;
  Minor_ellipsefin=0;
  Orien_Anglefin=0;
  Perimeterfin=0;
  
end

cluster_watern;
%[Lbw,n] = bwlabeln(erodedBWch4);
[La,numa] = bwlabel(segCluster);
[sx,sy]=size(segCluster);
%nana5=zeros(sx,sy);
a=La;
element=zeros(sx,sy,numa);
Tot=numa;
Xa = sprintf('Total de clusters %d ',numa);
disp(Xa)

x=1;

disp('inicia cuenta de Clusters')
while(x<(Tot)+1)


    
    % AQUI SE DEBERA CREAR LA IMAGEN DE CLUSTER DILATADA
    %L2=Lbw;
    %L2(L2~=x)=0;
    %se2 = strel('sphere',16);
    %L2 = imdilate(L2,se2);
    %factorb=L2; 
    
    a(La~=x)=0;
    a(La==x)=1;
    element(:,:,x)=a(:,:);
    %a=La;
    factorb=element(:,:,x);
    segmentedb=factorb.*Itrue;
    %nana5=label2rgb(L2*512);
    %imshow(nana5)
    %imshow(segmentedb)
    
    % CREAR IMAGEN CON PRIMER CLUSTER Y LOS DEMAS ELIMINADOS
    L=logical(factorb);
    box=regionprops(L,'Area', 'BoundingBox', 'centroid','Perimeter'); 
    areasas = cat(1, box.Area);
    Centroidc = cat(1, box.Centroid);
    bbox = cat(1, box.BoundingBox);
    Pericc = cat(1, box.Perimeter);
    bboxint=int16(bbox);
    %figure(11)
    %imshow(segmentedb)
    %title('actual cluster proicessing')
    % DILATAR
    
    % OBTENER PROPIEDADES DE RESPECTIVO CLUSTER
    
    
    %****************************
    % Calibrar
    if(areasas(1)>(8*153)) % solo procesa un tamaño definido por 14 veces una celula
           
        %Izzz = imcrop(erodedBW,bbox(x,:));
        %Izvv = imcrop(eBWcluster,bbox(x,:));
        Izzz = La(bboxint(1,2):bboxint(1,2)+bboxint(1,4)-1, bboxint(1,1):bboxint(1,1)+bboxint(1,3)-1);
        Izvv = segmentedb(bboxint(1,2):bboxint(1,2)+bboxint(1,4)-1, bboxint(1,1):bboxint(1,1)+bboxint(1,3)-1);
        compa=bwconncomp(Izvv);
        test1=regionprops(compa,'Area');  
        area_cel_c = cat(1, test1.Area);
        [axax,ayay]=size(area_cel_c);
        
        if(max(max(Izvv))~=0)
            [logics,eBW2,L4,La] = tacut(Izzz,Izvv);
            
            nana5(bboxint(1,2):bboxint(1,2)+bboxint(1,4)-1, bboxint(1,1):bboxint(1,1)+bboxint(1,3)-1,:)=L4;
            
            % logics ==1 indica que no hay celulas adentro
            if(logics==1)
              if(axax>1)
              %Diciembre2018  
              
              %rectangle('Position', [bbox(1,1),bbox(1,2),bbox(1,3),bbox(1,4)],...
              %   'EdgeColor','r','LineWidth',1 )
              cuentaclusters=cuentaclusters+1;
              centroidf(cuentaclusters,:)=Centroidc(1,:);
              pericluster(cuentaclusters,1)=Pericc(1);
              areacluster(cuentaclusters,1)=areasas(1);
              % figure
              % subplot(3,1,1)
              % imshow(Izzz)
              % subplot(3,1,2)
              % imshow(Izvv)
              % subplot(3,1,3)
              % imshow(L4)
              % pause
              ccn = bwconncomp(segmentedb);
              statsh3=regionprops(ccn, 'Area');
              areash3 = cat(1, statsh3.Area);
              %Th=Mina*cell_s;
              numPixels3 = cellfun(@numel,ccn.PixelIdxList);

                while (min(numPixels3)< Th)
                  [biggest,idx] = min(numPixels3);
                   segmentedb(ccn.PixelIdxList{idx}) = 0;
                   ccn = bwconncomp(segmentedb);
                   numPixels3 = cellfun(@numel,ccn.PixelIdxList);
                end
              cc = bwconncomp(segmentedb);
              t_cells(cuentaclusters,1)  = cc.NumObjects;          
              % subplot(1,2,2)
              %imshow(L4)
              % pause
              
              % a continuacion se analizan las celulas dentro de
              % correspondiente cluster
              
              

               %AQUI CORREGIR
              
              C = bwconncomp(segmentedb);
              statsh3=regionprops(C, 'Area', 'Eccentricity','centroid', 'MajorAxisLength','MinorAxisLength','Orientation','Perimeter','BoundingBox');
              Area = cat(1, statsh3.Area);
              Eccen = cat(1, statsh3.Eccentricity);
              Centroid = cat(1, statsh3.Centroid);
              Mayor_ellipse = cat(1, statsh3.MajorAxisLength);
              Minor_ellipse = cat(1, statsh3.MinorAxisLength);
              Orien_Angle = cat(1, statsh3.Orientation);
              Perimeter = cat(1, statsh3.Perimeter);
              bboxclu = cat(1, statsh3.BoundingBox);
              bboxintclu=int16(bboxclu);
              
              vectordetected(cuentaclusters)=x;
               % este indice es un cluster valido
              
              % Definido en challenge2
              % Th=Mina*cell_s;
              %Th=Mina*mean(areash2);
              % aqui se determina el valor promedio de las celulas
              % IMPORTANTE -------------------------------------
               
              B=size(Area);
              
              clear Round;
              clear yfitC;
              count=1;
              flag1=0;
              disp('inicia cuenta de celulas dentro de cluster')
              if (modoe==0)
                    
                    figure(1)
                    subplot(3,2,1)
                    imshow(densitymatrix5_h);
                    title('Hue region density near cell');
                    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                    %   'EdgeColor','r','LineWidth',1 )
                    subplot(3,2,2)
                    imshow(densitymatrix5_v);
                    title('V region density near cell');
                    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                    %   'EdgeColor','r','LineWidth',1 )
                    subplot(3,2,3)
                    imshow(densitymatrix);
                    title('Back_density');
                    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                    %   'EdgeColor','r','LineWidth',1 )
          
                    subplot(3,2,4)
                    imshow(densitymatrix2);
                    title('Positive_density');
                    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                    %  'EdgeColor','r','LineWidth',1 )   
                    subplot(3,2,5)
                    imshow(densitymatrix3);
                    title('BBack');
                    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                    %   'EdgeColor','r','LineWidth',1 )    
                    subplot(3,2,6)
                    imshow(densitymatrix4);
                    title('BaD');
                    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                    %  'EdgeColor','r','LineWidth',1 ) 
                    disp('Processing clusters')
              end
              while(count<(B(1,1)+1))
                
                    
                  
                Xa = sprintf('Counting Total: %d  actual cell: %d  actual cluster: %d',B(1,1), count,cuentaclusters);
                disp(Xa)
                
                if (c1>9)
                   if(c1>99)
                     if(c1>999)
                        filenamecc = ['Cell_',num2str(c1),'.tif'];
                        filenamec2 = ['Cellmask_',num2str(c1),'.tif'];
                     else
                        filenamecc = ['Cell_0',num2str(c1),'.tif'];
                        filenamec2 = ['Cellmask_0',num2str(c1),'.tif'];
                     end
                   else
                      filenamecc = ['Cell_00',num2str(c1),'.tif'];
                      filenamec2 = ['Cellmask_00',num2str(c1),'.tif'];
                   end
                else
                   filenamecc = ['Cell_000',num2str(c1),'.tif'];
                   filenamec2 = ['Cellmask_000',num2str(c1),'.tif'];
                end
                
                % compute the roundness metric
                if(Perimeter(count)==0)
                    Perimeter(count)=1;
                end
                Round(count) = 4*pi*Area(count)/Perimeter(count)^2;
                
                % *******************************************
                % A continuacion se calcula textura de celulas dentro de
                % cluster    y se calcula la densidad background
                Itexturec = testi(bboxintclu(count,2):bboxintclu(count,2)+bboxintclu(count,4)-1, bboxintclu(count,1):bboxintclu(count,1)+bboxintclu(count,3)-1);
                
                H1 = hImage(bboxintclu(count,2):bboxintclu(count,2)+bboxintclu(count,4)-1, bboxintclu(count,1):bboxintclu(count,1)+bboxintclu(count,3)-1);
                S1 = sImage(bboxintclu(count,2):bboxintclu(count,2)+bboxintclu(count,4)-1, bboxintclu(count,1):bboxintclu(count,1)+bboxintclu(count,3)-1);
                V1 = vImage(bboxintclu(count,2):bboxintclu(count,2)+bboxintclu(count,4)-1, bboxintclu(count,1):bboxintclu(count,1)+bboxintclu(count,3)-1);
                
                Hei = he(bboxintclu(count,2):bboxintclu(count,2)+bboxintclu(count,4)-1, bboxintclu(count,1):bboxintclu(count,1)+bboxintclu(count,3)-1,:);
                He_bw= Itrue(bboxintclu(count,2):bboxintclu(count,2)+bboxintclu(count,4)-1, bboxintclu(count,1):bboxintclu(count,1)+bboxintclu(count,3)-1,:);
                
                
                hz=size(nonzeros(H1));
                sz=size(nonzeros(S1));
                vz=size(nonzeros(V1));
                
                % These are cell colour parameters
                Ha1(count)=sum(sum(H1))/hz(1,1);
                Sa1(count)=sum(sum(S1))/sz(1,1);
                Va1(count)=sum(sum(V1))/vz(1,1);
                
                glcmc = graycomatrix(Itexturec,'Offset',[2 0;0 2]);
                statsc = graycoprops(glcmc,{'contrast','homogeneity'});
                
                
                % These are texture cell parameters
                text_contrastc(count,:)=statsc.Contrast;
                text_homogenetyc(count,:)=statsc.Homogeneity;
                
                index1c=uint16(Centroid(count,1));
                index2c=uint16(Centroid(count,2));
                back_densityc(count) = densitymatrix(index2c,index1c);
                positive_densityc(count) = densitymatrix2(index2c,index1c);
                Bbackc(count)= densitymatrix3(index2c,index1c);
                BaDc(count)= densitymatrix4(index2c,index1c);
                
                % Histogram regional parameter is a 4 dimension vector
                crHv(count,:)=density_histogramH(index2c,index1c,:);
                crSv(count,:)=density_histogramS(index2c,index1c,:);
                crVv(count,:)=density_histogramV(index2c,index1c,:);
                
                % these are intercell colour parameters
                intercellH(count)=densitymatrix5_h(index2c,index1c);
                intercellS(count)=densitymatrix5_s(index2c,index1c);
                intercellV(count)=densitymatrix5_v(index2c,index1c);
                
                if (modoe==0)
                
                    figure(2)
                    hold off
                    subplot(2,2,1)
                    imshow(Hei)
                    title('Selected cell');
                    subplot(2,2,2)
                    hold off
                    imshow(He_bw)
                    title('Cell mask')
                    
                figure(3)
                imshow(he);
                hold on;
                rectangle('Position', [bboxclu(count,1),bboxclu(count,2),bboxclu(count,3),bboxclu(count,4)],...
                   'EdgeColor','r','LineWidth',1 )    
                   
                                                       
                   prompt = 'Indicate lymphocyte=1, Malignant=2, Normal=3, Not a Cell=4, Skip=5, Quit=6? ';
                   dG(count,1) = input(prompt);
                   if(dG(count,1)==1)
                     yfitC(count,:)='lymphocyt';
                     %filenamecc = ['C_lym_',num2str(c1),'.tif']; 
                     nam(count,:)=filenamecc
                   end
                   if(dG(count,1)==2)
                     yfitC(count,:)='Malignant';
                     %filenamecc = ['C_mal_',num2str(c1),'.tif'];
                     nam(count,:)=filenamecc;
                   end
                   if (dG(count,1)==3)
                     yfitC(count,:)='Normalepy';
                     %filenamecc = ['C_nor_',num2str(c1),'.tif'];
                     nam(count,:)=filenamecc;
                   end
                   if (dG(count,1)==4)
                     yfitC(count,:)='Not_cell_';
                     %filenamecc = ['C_not_',num2str(c1),'.tif'];
                     nam(count,:)=filenamecc;
                   end
                   if (dG(count,1)==5)
                     flag5=5;
                     flag1=1;
                     c1=c1-1;
                     %nam(count,:)=filenamecc;
                   end
                   if (dG(count,1)==6)
                     flag5=6;
                     flag1=1;
                     c1=c1-1;
                     %nam(count,:)=filenamecc;
                   end
                     count=count+1;
                     c1=c1+1;
                     imwrite(Hei,filenamecc)
                     imwrite(He_bw,filenamec2)
                   if (flag1==1)
                     flag1=count-2;
                     count=B(1,1)+1;
                     disp('checando')
                   end
                else
                   count=count+1;  
                end
                
              end
              
              if(flag1==0)
                  flag1=count-1;
              end
              
              Xa = sprintf('Cluster completed, cells included: %d ',flag1);
              disp(Xa)
              
              if(B(1,1)>0)
              clear dG;
              
              Round=Round';
              back_densityc=back_densityc';
              positive_densityc=positive_densityc';
              Bbackc=Bbackc';
              BaDc=BaDc';
              intercellH=intercellH';
              intercellS=intercellS';
              intercellV=intercellV';
              Ha1=Ha1';
              Sa1=Sa1';
              Va1=Va1';
              
              
              
              % *******************************************
              % *****************************************
              % calcular disytancia entre Centroid y centroidf
              %centr_dist=
              Clus_cen=ones(B(1,1),1).*centroidf(cuentaclusters,:);
              centr_dist=lvdt(Clus_cen(:,1),Clus_cen(:,2),Centroid(:,1),Centroid(:,2));
              T_clusters=cuentaclusters.*ones(B(1,1),1);
              A_clusters=areasas.*ones(B(1,1),1);
              
              
              cell_unit=ones(B(1,1),1);
 
              meanIpar_H=cell_unit*meanItotH;
              meanIpar_S=cell_unit*meanItotS;
              meanIpar_V=cell_unit*meanItotV;
              Ppar=cell_unit*denp;
              Bpar=cell_unit*denb; 
              Bapar=cell_unit*denba;
              Epar=cell_unit*denE;
              
              % solo muestra el promedio de los datos de cada cluster
              Mean_area(cuentaclusters,1)=mean(Area);
              Mean_perimeter(cuentaclusters,1)=mean(Perimeter);
              Mean_round(cuentaclusters,1)=mean(Round);
              Mean_eccentricity(cuentaclusters,1)=mean(Eccen);
              Mean_major_ellipse(cuentaclusters,1)=mean(Mayor_ellipse);
              Mean_minor_ellipse(cuentaclusters,1)=mean(Minor_ellipse);
              Mean_orientation_angle(cuentaclusters,1)=mean(Orien_Angle);
              %size(Tab_clust)
              % aqui ya se tienen los datos completos de el cluster n
              if(flagelo==0)
                  conteo_matrix=conteo_matrix-1;
              end
              if (modoe)
                Tab_clust(conteo_matrix+1:conteo_matrix+t_cells(cuentaclusters,1),:)=table(back_densityc, positive_densityc, Bbackc, BaDc, intercellH, intercellS, intercellV,crHv(:,1), crHv(:,2), crHv(:,3), crHv(:,4), crSv(:,1), crSv(:,2), crSv(:,3), crSv(:,4), crVv(:,1), crVv(:,2),crVv(:,3),crVv(:,4),Area, Eccen, Round, Centroid(:,1), Centroid(:,2), Mayor_ellipse, Minor_ellipse, Orien_Angle, Perimeter,text_contrastc(:,1),text_contrastc(:,2),text_homogenetyc(:,1), text_homogenetyc(:,2),Ha1,Va1,Sa1,T_clusters,A_clusters,centr_dist,meanIpar_H, meanIpar_S,meanIpar_V,Ppar,Bpar,Bapar,Epar);
              else
                  flag1=flag1+1;
                  
                  Tab_clust(conteo_matrix+1:conteo_matrix+flag1-1,:)=table(back_densityc(1:flag1-1), positive_densityc(1:flag1-1), Bbackc(1:flag1-1), BaDc(1:flag1-1), intercellH(1:flag1-1), intercellS(1:flag1-1), intercellV(1:flag1-1), crHv(1:flag1-1,1), crHv(1:flag1-1,2), crHv(1:flag1-1,3), crHv(1:flag1-1,4), crSv(1:flag1-1,1), crSv(1:flag1-1,2), crSv(1:flag1-1,3), crSv(1:flag1-1,4), crVv(1:flag1-1,1), crVv(1:flag1-1,2),crVv(1:flag1-1,3),crVv(1:flag1-1,4),Area(1:flag1-1), Eccen(1:flag1-1), Round(1:flag1-1), Centroid(1:flag1-1,1), Centroid(1:flag1-1,2), Mayor_ellipse(1:flag1-1), Minor_ellipse(1:flag1-1), Orien_Angle(1:flag1-1), Perimeter(1:flag1-1),text_contrastc(1:flag1-1,1),text_contrastc(1:flag1-1,2),text_homogenetyc(1:flag1-1,1),text_homogenetyc(1:flag1-1,2),Ha1(1:flag1-1),Va1(1:flag1-1),Sa1(1:flag1-1),T_clusters(1:flag1-1),A_clusters(1:flag1-1),centr_dist(1:flag1-1),meanIpar_H(1:flag1-1), meanIpar_S(1:flag1-1),meanIpar_V(1:flag1-1),Ppar(1:flag1-1),Bpar(1:flag1-1),Bapar(1:flag1-1),Epar(1:flag1-1),yfitC,nam);
                  
              end 
              
              % aqui llena con ceros 
              %conteo_matrix=conteo_matrix+t_cells(cuentaclusters,1);
              conteo_matrix=conteo_matrix+flag1-1;
              end
              end
            end
        end
        
        

    end
    %limpiar variables
    clear Round;
    clear back_densityc;
    clear positive_densityc;
    clear Bbackc;
    clear BaDc;
    clear intercellH;
    clear intercellS;
    clear intercellV;
    clear text_contrastc;
    clear text_homogenetyc;
    clear Ha1;
    clear Sa1;
    clear Va1;
    clear Area;
    clear Eccen;
    clear Centroid;
    clear Mayor_ellipse;
    clear Minor_ellipse;
    clear Orien_Angle;
    clear Perimeter;
    clear nam;
    clear crHv;
    clear crSv;
    clear crVv;
    [La,numa] = bwlabel(segCluster);
    clear a;
    a=La;

x=x+1;
if (flag5==6) % en este modo sale definitivamente del programa
    x=Tot+1;
end
end
distancea=0;
B=size(vectordetected);
for x=1:B(1,2)
    for y=1:B(1,2)
        if(y~=x)
            dataa(1,:)=centroidf(x,:);
            dataa(2,:)=centroidf(y,:);
            distancea(y,1)=pdist(dataa);
            %line(dataa(1,:), dataa(2,:));
            %pause
        end
    end
    [min_dist(x,1) min_dist(x,2)]=min(distancea);
   
end
min_d=min_dist(:,1);


%DIC2018
%hold off

tot_cluster=cuentaclusters;
if(tot_cluster<0)
  if(length(areasas)>0)
    %ToTab2_3=table(centroidf, min_d, Mean_area,Mean_perimeter,Mean_round,Mean_eccentricity,Mean_major_ellipse,Mean_minor_ellipse,Mean_orientation_angle, pericluster, areacluster );
  end
end

if (modoe==0)  
    hold off;
end

 %Tab_allcells=table(Areafin, back_density, positive_density, Bback, BaD, Eccenfin, Roundfin, Centroidfin, Mayor_ellipsefin, Minor_ellipsefin, Orien_Anglefin, Perimeterfin,text_contrast,text_homogenety,hValue,sValue,vValue);
 %writetable(Tab_allcells,'tablaclusters.txt');
 %save Tab_allcells Tab_allcells
 save Tab_clust Tab_clust;

 disp('Table Tab_allcells completed...PAUSED')
 






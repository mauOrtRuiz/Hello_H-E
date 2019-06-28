
% routine to process density
% Must be run before massive 2
% Background is BWHe this is Eosin Region obtained from IMGa
% eeBW is the total positive cells including cells out of positive DAB
%region
% Bluebackground density is determined in Bback table data
% Background density is determined in BaD table data

[J1, K1] = size(BWHe);
size_cell=30;
load modoe modoe;
xn=1;
yn=1;
imagen_intercellh=hImage.*imcomplement(Itrue);
imagen_intercells=sImage.*imcomplement(Itrue);
imagen_intercellv=vImage.*imcomplement(Itrue);

if flag3==1 % primera vez que genera tabla de salida
     c1=1;
   else  %lee tabla de salida existente
     Ta=readtable('All_trained_data.txt');
     c1=height(Ta(:,1))+1;
end
  

%figure
%imshow(BWHe);
%hold on
for x1=1:size_cell:J1
    if((x1+size_cell-1)<J1)
        limit1=x1+size_cell-1;
    else
        limit1=J1;
    end
    for y1=1:size_cell:K1
        if((y1+size_cell-1)<K1)
                limit2=y1+size_cell-1;
        else
                limit2=K1;
        end
         
         % Inter-cell region is calculated by mean of non-zero elements
         numPixelsh= sum(sum(imagen_intercellh(x1:limit1,y1:limit2)));
         %s=size(imagen_intercellh(x1:limit1,y1:limit2));
         hh=size(nonzeros(imagen_intercellh(x1:limit1,y1:limit2)));
         AVIh=numPixelsh/hh(1,1);
         
         numPixelss= sum(sum(imagen_intercells(x1:limit1,y1:limit2)));
         %s=size(imagen_intercells(x1:limit1,y1:limit2));
         hs=size(nonzeros(imagen_intercells(x1:limit1,y1:limit2)));
         AVIs=numPixelss/hs(1,1);
         
         numPixelsv= sum(sum(imagen_intercellv(x1:limit1,y1:limit2)));
         %s=size(imagen_intercellv(x1:limit1,y1:limit2));
         hv=size(nonzeros(imagen_intercellv(x1:limit1,y1:limit2)));
         AVIv=numPixelsv/hv(1,1);
         
         
         %total number of non zero
         %nonzeros(imagen_intercellh(x1:limit1,y1:limit2))
         
         % PINK REGION
         numWhitePixels=sum(sum(BWHe(x1:limit1,y1:limit2)));
         numBlackPixels=sum(sum(~BWHe(x1:limit1,y1:limit2)));
         
         % BLUE CELLS REGION
         numWhitePixels2=sum(sum(Itrue(x1:limit1,y1:limit2)));
         numBlackPixels2=sum(sum(~Itrue(x1:limit1,y1:limit2)));
         
         %BLUE REGION
         numWhitePixels3=sum(sum(erodedBWch1(x1:limit1,y1:limit2)));
         numBlackPixels3=sum(sum(~erodedBWch1(x1:limit1,y1:limit2)));
         
         % Background region
         numWhitePixels4=sum(sum(BWDsub(x1:limit1,y1:limit2)));
         numBlackPixels4=sum(sum(~BWDsub(x1:limit1,y1:limit2)));
         
         % Histogram of Regional Area in HSV parameters
         heH=Itot(x1:limit1,y1:limit2,1);
         heS=Itot(x1:limit1,y1:limit2,2);
         heV=Itot(x1:limit1,y1:limit2,3);
         [Jreg, Kreg, Mreg] = size(heH);
         
         Tarreg=Jreg*Kreg;
         [crH,brH] = imhist(heH,4);
         [crS,brS] = imhist(heS,4);
         [crV,brV] = imhist(heV,4);
         crH=crH/Tarreg;
         crS=crS/Tarreg;
         crV=crV/Tarreg;
         
         %rectangle('Position', [y1,x1,size_cell,size_cell],...
         %        'EdgeColor','r','LineWidth',2 )
         
         AVRG=numWhitePixels/(size_cell*size_cell);
         densitymatrix(x1:limit1,y1:limit2)=AVRG;
         
         AVRG2=numWhitePixels2/(size_cell*size_cell);
         densitymatrix2(x1:limit1,y1:limit2)=AVRG2;
         
         AVRG3=numWhitePixels3/(size_cell*size_cell);
         densitymatrix3(x1:limit1,y1:limit2)=AVRG3;
         
         AVRG4=numWhitePixels4/(size_cell*size_cell);
         densitymatrix4(x1:limit1,y1:limit2)=AVRG4;
         
         densitymatrix5_h(x1:limit1,y1:limit2)=AVIh;
         densitymatrix5_s(x1:limit1,y1:limit2)=AVIs;
         densitymatrix5_v(x1:limit1,y1:limit2)=AVIv;
         
         % Density Histogram is a regional histogram nearby the selected
         % Cell
         density_histogramH(x1:limit1,y1:limit2,1)=crH(1);
         density_histogramH(x1:limit1,y1:limit2,2)=crH(2);
         density_histogramH(x1:limit1,y1:limit2,3)=crH(3);
         density_histogramH(x1:limit1,y1:limit2,4)=crH(4);
         
         density_histogramS(x1:limit1,y1:limit2,1)=crS(1);
         density_histogramS(x1:limit1,y1:limit2,2)=crS(2);
         density_histogramS(x1:limit1,y1:limit2,3)=crS(3);
         density_histogramS(x1:limit1,y1:limit2,4)=crS(4);
         
         density_histogramV(x1:limit1,y1:limit2,1)=crV(1);
         density_histogramV(x1:limit1,y1:limit2,2)=crV(2);
         density_histogramV(x1:limit1,y1:limit2,3)=crV(3);
         density_histogramV(x1:limit1,y1:limit2,4)=crV(4);
         
         AV_matriz(xn,yn)=AVRG;
         AV_matriz2(xn,yn)=AVRG2;
         AV_matriz3(xn,yn)=AVRG3;
         AV_matriz4(xn,yn)=AVRG4;
         
         %AV_matrix5_s(nx,yn)=AVIs;
         yn=yn+1;
         if(limit2==K1)
             yn=1;
         end
          
    end
   
    xn=xn+1;
end

%figure
%DIC 18
%imshow(densitymatrix)
%hold on
back_density=0;
positive_density=0;
Bback=0;
BaD=0;
if (modoe==0)
    figure(4)
    
    subplot(1,2,1)
    imshow(densitymatrix5_h);
    title('Hue region density near cell');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )
    subplot(1,2,2)
    imshow(imagen_intercellh)
    title('Hue region density near cell');
    
    figure(5)
    subplot(1,2,1)
    imshow(densitymatrix5_v);
    title('V region density near cell');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )
    subplot(1,2,2)
    imshow(imagen_intercellv)
    title('V region density near cell');
    
    figure(6)
    subplot(1,2,1)
    imshow(densitymatrix5_s);
    title('S region density near cell');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )
    subplot(1,2,2)
    imshow(imagen_intercells)
    title('S region density near cell');
    
    figure(7)
    subplot(1,2,1)
    imshow(densitymatrix);
    title('Back_density (Pink - eosin)');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )
    subplot(1,2,2)
    imshow(BWHe)
    title('Background image eosin')
    
    
    figure(8)
    subplot(1,2,1)
    imshow(densitymatrix2);
    title('Positive_density');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )   
    subplot(1,2,2)
    imshow(Itrue)
    title('Positive cells image');
    
    figure(9)
    subplot(1,2,1)
    imshow(densitymatrix3);
    title('BBack');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )    
    subplot(1,2,2)
    imshow(erodedBWch1)
    title(' Cluster Region');
    
    figure(10)
    subplot(1,2,1)
    imshow(densitymatrix4);
    title('BaD');
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 ) 
    subplot(1,2,2)
    imshow(BWDsub)
    title('Background Region')
    
    %figure(2)
    %imshow(he);
    %hold on;
    %rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
    %             'EdgeColor','r','LineWidth',1 )
end
flag1=0;
B1=size(Centroid); %corregir Centroid
count=1;
flagelo=1;
if (B1(1,1)==0)
    flagelo=0;
end
 while(count<(B1(1,1)+1))
     
     % compute the density metric
     index1=uint16(Centroid(count,1));
     index2=uint16(Centroid(count,2));
     back_density(count) = densitymatrix(index2,index1);
     positive_density(count) = densitymatrix2(index2,index1);
     Bback(count)= densitymatrix3(index2,index1);
     BaD(count)= densitymatrix4(index2,index1);
     
     % Histogram regional parameter is a 4 dimension vector
     crHv(count,:)=density_histogramH(index2,index1,:);
     crSv(count,:)=density_histogramS(index2,index1,:);
     crVv(count,:)=density_histogramV(index2,index1,:);
     
     intercellHnc(count)=densitymatrix5_h(index2,index1);
     intercellSnc(count)=densitymatrix5_s(index2,index1);
     intercellVnc(count)=densitymatrix5_v(index2,index1);
     
     Hei = he(bbbox(count,2):bbbox(count,2)+bbbox(count,4)-1, bbbox(count,1):bbbox(count,1)+bbbox(count,3)-1,:);
     He_bw= Itrue(bbbox(count,2):bbbox(count,2)+bbbox(count,4)-1, bbbox(count,1):bbbox(count,1)+bbbox(count,3)-1,:);
     
     if (modoe==0)
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
                    
          figure(2)
          hold off
          subplot(2,2,1)
          imshow(Hei)
          title('Selected cell');
          hold off
          subplot(2,2,2)
          imshow(He_bw)
          title('Cell mask')
                             
          figure(4)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          
          figure(5)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          figure(6)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          figure(7)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          figure(8)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          figure(9)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          figure(10)
          subplot(1,2,1)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          subplot(1,2,2)
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          hold off
          
          
          figure(3)
          clf
          imshow(he);
          hold on;
          rectangle('Position', [bbbox(count,1),bbbox(count,2),bbbox(count,3),bbbox(count,4)],...
                 'EdgeColor','r','LineWidth',1 )
          %Area(count)
          %back_density(count)
          %positive_density(count)
          %Bback(count)
          %BaD(count)
          %intercellHnc(count)
          %intercellSnc(count)
          %intercellVnc(count)
          
          prompt = 'Indicate lymphocyte=1, Malignant=2, Normal=3?, Not a Cell=4, Skip=5, Quit=6? ';
          dG(count,1) = input(prompt);
          if(dG(count,1)==1)
            yfitC(count,:)='lymphocyt';
            %filenamecc = ['C_lym_',num2str(c1),'.tif']; 
            nam(count,:)=filenamecc;
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
           if (dG(count,1)>6)
            yfitC(count,:)='Not_cell_';
            %filenamecc = ['C_not_',num2str(c1),'.tif'];
            nam(count,:)=filenamecc;
          end
          if (dG(count,1)==5)
            flag5=5;
            flag1=1; 
            disp('cinco')
          end
          if (dG(count,1)==6)
            flag5=6;
            flag1=1;  
          end
          count
          count=count+1;
          c1=c1+1;
          imwrite(Hei,filenamecc)
          imwrite(He_bw,filenamec2)
          
          if (flag1==1)
            flag1=count-1;
            c1=c1-1;
            count=B1(1,1)+1;
            disp('checando')
          end
     else
         count=count+1;
     end
     %  Aqui se debera generar tabla de celulas externas
     
     
     %Tab_clust(count,:)=table(Area, back_density, positive_density, Bback, BaD, intercellHnc, intercellSnc, intercellVnc, Eccen, Round, Centroid, Mayor_ellipse, Minor_ellipse, Orien_Angle, Perimeter,text_contrastc,text_homogenetyc,Ha1,Va1,Sa1,T_clusters,A_clusters,centr_dist);
     
    % pause
 
 end
 if(flag1==0)
   flag1=count;
 end
 
 if(flagelo)
 
 %hold off;
 %figure
 %imshow(densitymatrix2)
 back_density=back_density';
 positive_density=positive_density';
 Bback=Bback';
 BaD=BaD';
 intercellHnc=intercellHnc';
 intercellSnc=intercellSnc';
 intercellVnc=intercellVnc';
 %These are final cell parameters
 % falta poner ceros a parametros de clusters que noaplica
 cell_zer=zeros(B1(1,1),1);
 cell_unit=ones(B1(1,1),1);
 
 meanIpar_H=cell_unit*meanItotH;
 meanIpar_S=cell_unit*meanItotS;
 meanIpar_V=cell_unit*meanItotV;
 Ppar=cell_unit*denp;
 Bpar=cell_unit*denb; 
 Bapar=cell_unit*denba;
 Epar=cell_unit*denE;
 
 T_clusters=cell_zer;
 A_clusters=cell_zer;
 centr_dist=cell_zer;
 if (modoe)
  Tab_clust(:,:)=table(back_density, positive_density, Bback, BaD, intercellHnc, intercellSnc, intercellVnc, crHv(:,1), crHv(:,2), crHv(:,3), crHv(:,4), crSv(:,1), crSv(:,2), crSv(:,3), crSv(:,4), crVv(:,1),  crVv(:,2), crVv(:,3), crVv(:,4), Area, Eccen, Round, Centroid(:,1), Centroid(:,2), Mayor_ellipse, Minor_ellipse, Orien_Angle, Perimeter,text_contrast(:,1),text_contrast(:,2), text_homogenety(:,1),text_homogenety(:,2), hValue,vValue,sValue,T_clusters,A_clusters,centr_dist,meanIpar_H, meanIpar_S,meanIpar_V,Ppar,Bpar,Bapar,Epar);
 else
   %flag1=flag1+1; 
   Tab_clust(1:flag1-1,:)=table(back_density(1:flag1-1), positive_density(1:flag1-1), Bback(1:flag1-1), BaD(1:flag1-1), intercellHnc(1:flag1-1), intercellSnc(1:flag1-1), intercellVnc(1:flag1-1), crHv(1:flag1-1,1), crHv(1:flag1-1,2), crHv(1:flag1-1,3), crHv(1:flag1-1,4), crSv(1:flag1-1,1), crSv(1:flag1-1,2), crSv(1:flag1-1,3), crSv(1:flag1-1,4), crVv(1:flag1-1,1), crVv(1:flag1-1,2),crVv(1:flag1-1,3),crVv(1:flag1-1,4), Area(1:flag1-1), Eccen(1:flag1-1), Round(1:flag1-1), Centroid(1:flag1-1,1), Centroid(1:flag1-1,2), Mayor_ellipse(1:flag1-1), Minor_ellipse(1:flag1-1), Orien_Angle(1:flag1-1), Perimeter(1:flag1-1),text_contrast(1:flag1-1,1),text_contrast(1:flag1-1,2),text_homogenety(1:flag1-1,1),text_homogenety(1:flag1-1,2),hValue(1:flag1-1),vValue(1:flag1-1),sValue(1:flag1-1),T_clusters(1:flag1-1),A_clusters(1:flag1-1),centr_dist(1:flag1-1),meanIpar_H(1:flag1-1), meanIpar_S(1:flag1-1),meanIpar_V(1:flag1-1),Ppar(1:flag1-1),Bpar(1:flag1-1),Bapar(1:flag1-1),Epar(1:flag1-1),yfitC,nam);
 end
%These are final cell parameters
Areafin=Area;
Eccenfin=Eccen;
Roundfin=Round;
Centroidfin=Centroid;
Mayor_ellipsefin=Mayor_ellipse;
Minor_ellipsefin=Minor_ellipse;
Orien_Anglefin=Orien_Angle;
Perimeterfin=Perimeter;
% text_contrast
% text_homogenety
% back_density
% positive_density
%Tab_allcells=table(Areafin, back_density, positive_density, Bback, BaD, Eccenfin, Roundfin, Centroidfin, Mayor_ellipsefin, Minor_ellipsefin, Orien_Anglefin, Perimeterfin,text_contrast,text_homogenety,hValue,sValue,vValue);

 end
disp('Background parameters completed...')





% generate segmentation Purple
[J, K] = size(he);
Igc=Ib;
for x=1:J
    for y=1:K
       if(Igc(x,y)==0)
          Igc(x,y)=255; 
       end
        
    end
end



%this is a big enough mask
III=rgb2gray(Igc);
II=imbinarize(III,0.55);

%nueva tecnica de binarizacion blue aqui empieza:
IIIm=medfilt2(III);
IIIm=medfilt2(IIIm);
[J, K] = size(IIIm);
IIIm(J,K)=255;
%IIIm is a good gray imjge with reinforcement in cells density
IIIa=rgb2gray(he);
IIIm2=imcomplement(double(IIIm)/double(255));
IIIa2=imcomplement(double(IIIa)/double(255));

% Here operations with image density can be done
IIIne=imcomplement(IIIm2.*IIIa2*.8+IIIa2*0.30);

IIIgama = imadjust(IIIm,[],[],0.1);

%IIIm2
IIIg2=imcomplement(double(IIIgama)/double(255));
IIIne3=imcomplement(IIIne);
IIIne4=imcomplement(IIIne3+IIIg2);

grayKmeans;


final=imbinarize(IIIne4,(umbpi+umbpi*.2));

Itrue=imbinarize(IIIne4,(umbpi+umbpi*.2));
%figure
%imshowpair(he, Itrue, 'montage')
%pause
% AQUI HAY QUE PROPONER UN UMBRAL MAS ESTRICTO
Itrueout=imbinarize(IIIne4,(umbpi+umbpi*.20));


%IIIf=imcomplement(double(final));
%IIIfinew=imcomplement(IIIf*.5+IIIm2);
%IIIfinew=medfilt2(IIIfinew);
% nueva tecnica de binarizacion aqui termina 


IIc=imcomplement(II);
eBW=imfill(IIc,'holes');
 % -------------------------------
 % ------------------------------
 % en esta parte se trabaja con la imagen binarizada de la mask del color
 % azul esto se debe corregir para usar una mejor imagen binarizada
 % ---------------------------------
 % ----------------------------------

C2 = bwconncomp(eBW);
statsh2=regionprops(C2, 'Area');
areash2 = cat(1, statsh2.Area);
Th=Mina*cell_s;
%Th=Mina*mean(areash2);
% aqui se determina el valor promedio de las celulas



numPixels = cellfun(@numel,C2.PixelIdxList);

while (min(numPixels)< Th)
    [biggest,idx] = min(numPixels);
    eBW(C2.PixelIdxList{idx}) = 0;
    C2 = bwconncomp(eBW);
    numPixels = cellfun(@numel,C2.PixelIdxList);
end
%figure
eeBW=imfill(eBW,'holes');
%imshow(eeBW)
%title('SEGMENTED before watershed')

%now is filtered again this is only a mask:

%AQUI SE PROCESA NUEVA SERIE DADA EN ENERO 2019  eeBWag
se2 = strel('sphere',1); %AQUI ESTOY X3
eeBWag=imdilate(eBWcluster,se2);
eeBWag=imfill(eeBWag,'holes');
eeBWag=imerode(eeBWag,se2);
eeBWag=imdilate(eeBWag,se2);
eeBWag=imfill(eeBWag,'holes');
eeBWag=imerode(eeBWag,se2);

trueBW=rgb2gray(he).*uint8(eeBWag);

[J, K] = size(trueBW);
for x=1:J
    for y=1:K
       if(trueBW(x,y)==0)
          trueBW(x,y)=255; 
       end
        
    end
end






% Modifica 2019
%Itrue=imbinarize(vImage,0.40);
%Itrueout=imbinarize(vImage,0.55);
%Itrue=imbinarize(trueBW,0.6);
Itrue=imcomplement(Itrue);
Itrueout=imcomplement(Itrueout);


se = strel('sphere',1);
Itrue2 = imerode(Itrue,se);
Itrue2o = imerode(Itrueout,se);

Itrue=Itrue2;
Itrueout=Itrue2o;

C2 = bwconncomp(Itrue);
statsh2=regionprops(C2, 'Area');
areash2 = cat(1, statsh2.Area);
numPixels = cellfun(@numel,C2.PixelIdxList);

while (min(numPixels)< Th)
    [biggest,idx] = min(numPixels);
    Itrue(C2.PixelIdxList{idx}) = 0;
    C2 = bwconncomp(Itrue);
    numPixels = cellfun(@numel,C2.PixelIdxList);
end



C2 = bwconncomp(Itrueout);
statsh2=regionprops(C2, 'Area');
areash2 = cat(1, statsh2.Area);
numPixels = cellfun(@numel,C2.PixelIdxList);

while (min(numPixels)< Th)
    [biggest,idx] = min(numPixels);
    Itrueout(C2.PixelIdxList{idx}) = 0;
    C2 = bwconncomp(Itrueout);
    numPixels = cellfun(@numel,C2.PixelIdxList);
end


factor=erodedBWch1;

% ***************************
%
%  AQUI HAY UN ERROR 21 DE FEBRE
%  eeBW es la segmentacion de V value, no del nuevo metodo adaptable
segmentedclus=erodedBWch1.*eeBW;

%THIS IS ENERO 2019
%Itrue=eeBWag;
%Generates DATA from He positive cells

C = bwconncomp(Itrue);
statsh3=regionprops(C, 'Area', 'BoundingBox');
Area = cat(1, statsh3.Area);
bbbox = cat(1, statsh3.BoundingBox);
bboxint=int16(bbbox);

 % obtain the area calculation corresponding to label 'k'
 B=size(Area);
 
 for count=1:B(1,1)
     % compute the roundness metric
     
     if (Area(count)>3*cell_s)
         Iwater = Itrue(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
         
         D = -bwdist(~Iwater);
         Ld = watershed(D);
         bw2 = Iwater;
         bw2(Ld == 0) = 0;
         mask = imextendedmin(D,1);
         %imshowpair(Iwater,mask,'blend')
         D2 = imimposemin(D,mask);
         Ld2 = watershed(D2);
         bw3 = Iwater;
         bw3(Ld2 == 0) = 0;
         Itrue(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1)=bw3;
         %figure
         %imshow(bw3)
         %title('Here applies watershed')
         %pause
     end
 end
 
 
% Aqui aplica watershed de externas a cluster
C = bwconncomp(Itrueout);
statsh3=regionprops(C, 'Area', 'BoundingBox');
Area = cat(1, statsh3.Area);
bbbox = cat(1, statsh3.BoundingBox);
bboxint=int16(bbbox);

 % obtain the area calculation corresponding to label 'k'
 B=size(Area);
 
 for count=1:B(1,1)
     % compute the roundness metric
     
     if (Area(count)>3*cell_s)
         Iwater = Itrueout(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
         
         D = -bwdist(~Iwater);
         Ld = watershed(D);
         bw2 = Iwater;
         bw2(Ld == 0) = 0;
         mask = imextendedmin(D,1);
         %imshowpair(Iwater,mask,'blend')
         D2 = imimposemin(D,mask);
         Ld2 = watershed(D2);
         bw3 = Iwater;
         bw3(Ld2 == 0) = 0;
         Itrueout(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1)=bw3;
         %figure
         %imshow(bw3)
         %title('Here applies watershed')
         %pause
     end
 end
 
 Itrue=imfill(Itrue,'holes');
 
 %figure
 %subplot(1,2,1)
 %imshow(Itrue)
 %subplot(1,2,2)
 %imshow(eeBW)
 %title('comparison of watershed applied')

 
factora=erodedBWch1; 
segmented=factora.*Itrue;

restante=imcomplement(factora).*Itrueout;

eeBW=restante;

%Generates DATA from He positive cells FUERA DE LOS CLUSTERS
C = bwconncomp(eeBW);
statsh3=regionprops(C, 'Area', 'Eccentricity','centroid', 'MajorAxisLength','MinorAxisLength','Orientation','Perimeter','BoundingBox');
Area = cat(1, statsh3.Area);
Eccen = cat(1, statsh3.Eccentricity);
Centroid = cat(1, statsh3.Centroid);
Mayor_ellipse = cat(1, statsh3.MajorAxisLength);
Minor_ellipse = cat(1, statsh3.MinorAxisLength);
Orien_Angle = cat(1, statsh3.Orientation);
Perimeter = cat(1, statsh3.Perimeter);
bbbox = cat(1, statsh3.BoundingBox);
bboxint=int16(bbbox);

egray=rgb2gray(he);


% para calcular texturas de centros de celula se esta usando eeBW, esta es
% solo cell outside cluster
% corregir la siguiente linea Febrero 22/19
[J, K] = size(eeBW);
for x=1:J
    for y=1:K
       if(eeBW(x,y)==0)
          egray(x,y)=0; 
       end
        
    end
end

 hValue=0;
 vValue=0;
 sValue=0;
 Round=0;
 text_contrast(1,1)=0;
 text_contrast(1,2)=0;
 text_homogenety(1,1)=0;
 text_homogenety(1,2)=0;
 % obtain the area calculation corresponding to label 'k' and obtains
 % texture of every cell
 B=size(Area);
 for count=1:B(1,1)
     % compute the roundness metric
     if(Perimeter(count)==0)
        Perimeter(count)=1;
     end
     Round(count) = 4*pi*Area(count)/Perimeter(count)^2;
     
     Itexture = egray(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
     glcm = graycomatrix(Itexture,'Offset',[2 0;0 2]);
     stats = graycoprops(glcm,{'contrast','homogeneity'});
     
     H1 = hImage(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
     S1 = sImage(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
     V1 = vImage(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
             
     hz=size(nonzeros(H1));
     sz=size(nonzeros(S1));
     vz=size(nonzeros(V1));
                
                
     hValue(count)=sum(sum(H1))/hz(1,1);
     vValue(count)=sum(sum(S1))/sz(1,1);
     sValue(count)=sum(sum(V1))/vz(1,1);
     
     text_contrast(count,:)=stats.Contrast;
     text_homogenety(count,:)=stats.Homogeneity;
     
     
 end
 
 hValue=hValue';
 vValue=vValue';
 sValue=sValue';
 text_contrast;
 text_homogenety;
 
 
 Round=Round';
 disp('Morphology parameters completed...')
%ESTO SE EJECUTARA HASTA EL FINAL
%se = strel('sphere',11);
%clBW = imdilate(eeBW,se);
%figure
%imshow(clBW)
%title('Estimated Hematoxyline region dilated 11 pixels')

%Generates Data from Hematoxilyn nuclei
% First separate overlapp nuclei

% 

% Esta seccion se deberà moves hasta calcular background parameters
% Table to predict data  use: yfit = trainedClassifier.predictFcn(Tabn);
% Tabn=table(Area(1:B), Eccen(1:B), Round(1:B), Mayor_ellipse(1:B), Minor_ellipse(1:B), Orien_Angle(1:B), Perimeter(1:B));
%  load Tab;
%  [trainedClassifier, validationAccuracy] = trainClassifier(Tab);
% yfit = trainedClassifier.predictFcn(Tabn);
  %Tab2=table(AreaH, EccenH, RoundH, CentroidH, Mayor_ellipseH, Minor_ellipseH, Orien_AngleH, PerimeterH);
  %Tab1=table(Area, Eccen, Round, Centroid, Mayor_ellipse, Minor_ellipse, Orien_Angle, Perimeter, yfit);
  %writetable(Tab1,'DABdata.txt')
  %writetable(Tab2,'HeData.txt')
  %Summary of Statistics GLOBAL data
  %Ro=0;
  %Elo=0;
  %Sur=0;
  %Nu=0;
  %for x=1:B
  %    if (yfit(x)=='Rounded')
  %        Ro=Ro+1;
  %    end
  %    if (yfit(x)=='Elongated')
  %        Elo=Elo+1;
  %    end
  %    if (yfit(x)=='Surrounding')
  %        Sur=Sur+1;
  %    end
  %    if (yfit(x)=='Nuclei_visible')
  %        Nu=Nu+1;
  %    end
  %end
  
  % Final table  DAB Positive
  % Total DAB+Cells | Mean Area | Mean Perimeter |Mean Round | Mean Eccen |
  % Mean Major ellipse | Min Minor | Mean Angle |Total Round Cells | Total Elongated
  % cells | Total Surrounding | Total Nuclei Visible 
  %Total_DAB=B(1,1);
  %Mean_area=mean(Area);
  %Mean_perimeter=mean(Perimeter);
  %Mean_round=mean(Round);
  %Mean_eccentricity=mean(Eccen);
  %Mean_major_ellipse=mean(Mayor_ellipse);
  %Mean_minor_ellipse=mean(Minor_ellipse);
  %Mean_orientation_angle=mean(Orien_Angle);
  %Total_round=Ro;
  %Total_elongated=Elo;
  %Total_surrounding=Sur;
  %Total_nuclei_visible=Nu;
  
  %ToTab=table(Total_DAB,Mean_area, Mean_perimeter,Mean_round,Mean_eccentricity,Mean_major_ellipse,Mean_minor_ellipse,Mean_orientation_angle,Total_round,Total_elongated,Total_surrounding,Total_nuclei_visible);
  
  %
  %figure
  %subplot(2,2,1)
  %imshow(he)
  %title('original')
  %subplot(2,2,2)
  %imshow(restante)
  %title('not clustered')
  %subplot(2,2,3)
  %imshow(segmented)
  %title(' segmented based on original he image')
  %subplot(2,2,4)
  %imshow(Itrue)
  %title(' all segmented cells')
  %pause
  
  
  %figure
  %subplot(2,2,3)
  %imshow(erodedBWch1)
  %title('DAB region after dilate')
  %subplot(2,2,4)
  %imshow(restante)
  %title('not clustered')
  %subplot(2,2,2)
  %imshow(segmented)
  %title('Actual cluster cells')
  %subplot(2,2,1)
  %imshow(he)
  %title('original')
  disp('CorrectionsCh Finish');
  
 
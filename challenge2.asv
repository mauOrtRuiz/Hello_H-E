% Algorithm to alanyse a full WSI image by sections
% visualization of HSV histograms of first region
% Segmentation algorithms for DAB & HE sections (brown - blue)
% Calibration Parameters:
% V value threshold  0.3 for clusters,  0.5 for cells
% Cell size = 153  considere to recalibrate
% Outputs:
%  File HeData.txt    classification od DAB positive cells (macrophages)
%  Table  Tab1    Positive cells detected (bown) summary
%  Table  Tab2    Negative cells detected (blue) summary
%  Table  ToTab   This is a summary of all cells (note: blue obtained from
%  QuPath)


Vthr=0.7;
Vthrclust=0.55; % Umbral para celulas fuera de cluster
Vthrclust2=0.55; %umbral para celulas dentro de cluster
cell_s=150;
flag5=0;

% Minimum area size
% Cells < than 0.15 mean area are discarded
Mina=0.15;

% El resultado del todo el procesamiento son imagenes desplegadas
% Ademas se presentan las tablas:
%  DABdata.txt
%  HeData.txt
%

% This line is only for GT analysys
%load GT
% next three lines should be eliminated


he = train1;
[J1, K1, M1] = size(he);
%eliminate yellow border
Tar=J1*K1;

imshow(he)
% Next line only for ground truth
%he=GT;

%figure


hsvImage = rgb2hsv(he);
hImage = hsvImage(:, :, 1);
sImage = hsvImage(:, :, 2);
vImage = hsvImage(:, :, 3);
%subplot(2,2,1);
%hHist = histogram(hImage);
%subplot(2,2,2);
%sHist = histogram(sImage);
%subplot(2,2,3);
%vHist = histogram(vImage);

[J, K] = size(hImage);
khImage=reshape(hImage,[1,J*K]);
[idx,C] = kmeans(khImage',2);
binI=reshape(idx,[J,K]);
JJ=imbinarize(he);
I1=JJ(:,:,1);
I2=I1;
I3=I1;
I4=I1;
%repeat for each cluster for K-means
for x=1:J
    for y=1:K
        I1(x,y)=0;
        I2(x,y)=0;
        I3(x,y)=0;
        I4(x,y)=0;
        if (binI(x,y)==1)
            I1(x,y)=1;
        end
        if (binI(x,y)==2)
            I2(x,y)=1;
        end
        if (binI(x,y)==3)
            I3(x,y)=1;
        end
        
    end
end
% aqui se muestra segmentacion por HSV + Kmeans clustering
        

%figure
%imshow(I1)
%title('Colour segmentation by HSV transformation - K-means')
%figure
%imshow(I2)
%title('Colour segmentation by HSV transformation - K-means')
%figure
%imshow(I3)
%title('Colour segmentation by HSV transformation - K-means')


% Now segmentation by *lab transformation + K-means


%cform = makecform('srgb2lab');
%lab_he = applycform(he,cform);
%ab = double(lab_he(:,:,2:3));
%nrows = size(ab,1);
%ncols = size(ab,2);
%ab = reshape(ab,nrows*ncols,2);

%nColors = 3;
% repeat the clustering 3 times to avoid local minima
%[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
%                                      'Replicates',3);
                                                                  
%pixel_labels = reshape(cluster_idx,nrows,ncols);
% Thhis is new
lab_he = rgb2lab(he);
ab = lab_he(:,:,2:3); 
nrows = size(ab,1); 
ncols = size(ab,2); 
ab = reshape(ab,nrows*ncols,2);  
nColors = 3;
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...                                       
    'Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols); 
%imshow(pixel_labels,[]), title('image labeled by cluster index');
segmented_images = cell(1,3); 
rgb_label = repmat(pixel_labels,[1 1 3]);  
for k = 1:nColors     
    color = he;     
    color(rgb_label ~= k) = 0;     
    segmented_images{k} = color; 
end
I1=segmented_images{1};
I2=segmented_images{2};
I3=segmented_images{3};

%[BWblue1,maskedRGBImageblue1] = createMaskbb(I1);
%[BWblue2,maskedRGBImageblue2] = createMaskbb(I2);
%[BWblue3,maskedRGBImageblue3] = createMaskbb(I3);

[BWblue1,maskedRGBImageblue1] = createMaskblue2(I1);
[BWblue2,maskedRGBImageblue2] = createMaskblue2(I2);
[BWblue3,maskedRGBImageblue3] = createMaskblue2(I3);

[i,j,s] = find(BWblue1);
den1=(length(i)/Tar);

[i,j,s] = find(BWblue2);
den2=(length(i)/Tar);

[i,j,s] = find(BWblue3);
den3=(length(i)/Tar);

V1=[den1 den2 den3];
[Ma, Ia]=max(V1);
%figure

Ib=segmented_images{Ia};

%imshow(Iblue)
%title('Blue image detected june 2019')
%pause
% I is blue image then remaining two images are back and pink
if (Ia==1)
    Ipp=I2;
    Ipba=I3;
else
    if(Ia==2)
        Ipp=I1;
        Ipba=I3;
    else
        if(Ia==3)
            Ipp=I1;
            Ipba=I2;
        end
    end
end


[BWpink1,maskedRGBImagepink2] = createMaskpink(Ipp);
[BWpink2,maskedRGBImagepink3] = createMaskpink(Ipba);

[i,j,s] = find(BWpink1);
den1=(length(i)/Tar);

[i,j,s] = find(BWpink2);
den2=(length(i)/Tar);

V1=[den1 den2];
[Ma, Ia]=max(V1);
%figure

if (Ia==1)
    Ip=Ipp;
    Iba=Ipba;
else
    Ip=Ipba;
    Iba=Ipp;
end



%Iphsv=rgb2hsv(Ipp);

% Pink image H value
%meanIph=mean(mean(Iphsv(:,:,1)));
%[i,j,s] = find(Iphsv(:,:,1));
%denph=(length(i)/Tar);
%menph=mean(s);

% Pink image S value
%meanIps=mean(mean(Iphsv(:,:,2)));
%[i,j,s] = find(Iphsv(:,:,2));
%denps=(length(i)/Tar);
%menps=mean(s);


%%%%
%Ibahsv=rgb2hsv(Ipba);

% Back image H value
%meanIbh=mean(mean(Ibahsv(:,:,1)));
%[i,j,s] = find(Ibahsv(:,:,1));
%denph=(length(i)/Tar);
%menph=mean(s);

% Back image S value
%meanIbs=mean(mean(Ibahsv(:,:,2)));
%[i,j,s] = find(Ibahsv(:,:,2));
%denps=(length(i)/Tar);
%menps=mean(s);

%%%%


%V=nonzeros(Ibahsv(:,:,2));
%n11=size(V);
%h1=histogram(V,4);
%hDatab=h1.Values/n11(1,1);

%V=nonzeros(Iphsv);
%n11=size(V);
%h1=histogram(V,4);
%hDatap=h1.Values/n11(1,1);

%if (hDatap(1)<hDatab(1))
%    Iba=Ipp;
%    Ip=Ipba;
%else
%    Ip=Ipp;
%    Iba=Ipba;
%end

%V2=[meanIps meanIbs];
%[M2, I2a]=min(V2);
%if (I2a==1)
%    Iba=Ipp;
%    Ip=Ipba;
%else
%    if(I2a==2)
%        Ip=Ipp;
%        Iba=Ipba;
%    end
%end
        
%%

%mx=mean(mean(I1(:,:,2)));
%my=mean(mean(I2(:,:,2)));
%mz=mean(mean(I3(:,:,2)));
%load x
%load y
%load z
%x1=mean(x);
%y1=mean(y);
%z1=mean(z);
%if (abs(mx-x1)<abs(my-x1) & abs(mx-x1)<abs(mz-x1))
    % This is Pink%
%    Ip=I1;
%else
%    if (abs(my-x1)<abs(mx-x1) & abs(my-x1)<abs(mz-x1))
        % This is Pink
%        Ip=I2;
%    else
%        Ip=I3;
%    end
%end

% Next Comparision
%if (abs(mx-y1)<abs(my-y1) & abs(mx-y1)<abs(mz-y1))
    % This is blue
%    Ib=I1;
%else
%    if (abs(my-y1)<abs(mx-y1) & abs(my-y1)<abs(mz-y1))
        % This is blue
%        Ib=I2;
%    else
%        Ib=I3;
%    end
%end

% Next Comparision
%if (abs(mx-z1)<abs(my-z1) & abs(mx-z1)<abs(mz-z1))
    % This is back
%    Iba=I1;
%else
%    if (abs(my-z1)<abs(mx-z1) & abs(my-z1)<abs(mz-z1))
        % This is back
%        Iba=I2;
%    else
%        Iba=I3;
%    end
%end
%%


figure(1)
subplot(1,4,1)
imshow(Ip)
title('PINK Componenet');
subplot(1,4,2)
imshow(Ib)
title('BLUE Componenet');
subplot(1,4,3)
imshow(Iba)
title('Background Componenet');
% AQUI ES DONDE SE ELIMINA LA PREGUNTA MANUAL E INICIA EL AUTOMATISMO
% -------------------------------------------------------------------
%prompt = 'Indicar si esta satisfecho? teclear 1 si requiere ajuste manual';
%ddd = input(prompt)
%if(ddd==1)

%IMGa=segmented_images{1};
%figure
%imshow(IMGa)
%pause
%prompt = 'Indicate Blue=1, Brown=2, white background=3 ? ';
%d1 = input(prompt)
%if(d1==1)
%    BWHe = imbinarize(IMGa(:,:,1));
%    Ip=IMGa;
    
%end
%if(d1==2)
%    BWDa = imbinarize(IMGa(:,:,1));
%    Ib=IMGa;
%    Igc=IMGa;
%end
%if(d1==3)
%    BWDsub = imbinarize(IMGa(:,:,1));
%    BWDr=imcomplement(BWDsub);
%    Iba=IMGa;
    
%end

%figure
%IMGb=segmented_images{2};
%imshow(IMGb)
%pause

%prompt = 'Indicate Blue=1, Brown=2, white background=3 ? ';
%d2 = input(prompt)
%if(d2==1)
%    BWHe = imbinarize(IMGb(:,:,1));
%    Ip=IMGb;
    
%end
%if(d2==2)
%    BWDa = imbinarize(IMGb(:,:,1));
%    Ib=IMGb;
%    Igc=IMGb;
%end
%if(d2==3)
%    BWDsub = imbinarize(IMGb(:,:,1));
%    BWDr=imcomplement(BWDsub);
%    Iba=IMGb;
%end

%figure
%IMGc=segmented_images{3};
%imshow(IMGc)
%pause

%prompt = 'Indicate Blue=1, Brown=2, white background=3 ? ';
%d3 = input(prompt)
%if(d3==1)
%    BWHe = imbinarize(IMGc(:,:,1));
%    Ip=IMGc;
%end
%if(d3==2)
%    BWDa = imbinarize(IMGc(:,:,1));
%    Ib=IMGc;
%    Igc=IMGc;
%end
%if(d3==3)
%    BWDsub = imbinarize(IMGc(:,:,1));
%    BWDr=imcomplement(BWDsub);
%    title('Background component')
%    Iba=IMGc;
%end

%end
%  AQUI COPNCLUYO LA RUTINA PREVIA AL AUTOMATISMO
%disp('Magenta componenets')
%[mean(mean(Ipink(:,:,2)))  mean(mean(Iblue(:,:,2)))  mean(mean(Iback(:,:,2)))]

%disp('Green componenets')
%[mean(mean(Ipink(:,:,3)))  mean(mean(Iblue(:,:,3)))  mean(mean(Iback(:,:,3)))]


%densities of color sections
% february 2019

Iphsv=rgb2hsv(Ip);
meanIphsv=mean(mean(Iphsv(:,:,1)));
[i,j,s] = find(Iphsv(:,:,1));
denp=(length(i)/Tar);
menp=mean(s);


Ibhsv=rgb2hsv(Ib);
meanIbhsv=mean(mean(Ibhsv(:,:,1)));
[i,j,s] = find(Ibhsv(:,:,1));
denb=(length(i)/Tar);
meanb=mean(s);

Ibahsv=rgb2hsv(Iba);
meanIbahsv=mean(mean(Ibahsv(:,:,1)));
[i,j,s] = find(Ibahsv(:,:,1));
denba=(length(i)/Tar);


Itot=rgb2hsv(he);
meanItotH=mean(mean(Itot(:,:,1)));
meanItotS=mean(mean(Itot(:,:,2)));
meanItotV=mean(mean(Itot(:,:,3)));

% histogram values
[coR,bnR] = imhist(he(:,:,1),16);
[coG,bnG] = imhist(he(:,:,2),16);
[coB,bnB] = imhist(he(:,:,3),16);
coR=coR/Tar;
coG=coG/Tar;
coB=coB/Tar;

% Histogram of Eosin elements and background elements
%[ip,jp,sp] = find(Ip(:,:,1));
%[countsp1,centersp1] = hist(sp,16);


%Eosin color filter
[BW1m, maskRGB]=createMask(he);
denE=sum(sum(BW1m))/Tar;

% Variables globales:
%meanItotH   -> this is the mean H component of original image
%meanItotS   -> this is the mean S component of original image
%meanItotV   -> this is the mean V component of original image
%denp        -> density of only total  elements of pink image(from Kmeans)
%denb        -> density of only total  elements of blue image(from Kmeans)
%denba       -> density of only total  elements of back image(from Kmeans)
%denE        -> density of eosin color after a filter mask created

%coR         -> 16 elements vector with R image historgram
%coG         -> 16 elements vector with G image historgram
%coB         -> 16 elements vector with B image historgram

%crHv        -> 4 elements regional vector (dimension 4)of HSV histogram  
%csSv        -> 4 elements regional vector (dimension 4)of HSV histogram
%crVv        -> 4 elements regional vector (dimension 4)of HSV histogram


BWDa=imbinarize(Ib(:,:,1));   % Ib
BWHe = imbinarize(Ip(:,:,1)); % Ip
BWDsub = imbinarize(Iba(:,:,1)); %Iba
BWDr=imcomplement(BWDsub);


%imshow(pixel_labels,[]), title('image labeled by cluster index');
%pause
%segmented_images = cell(1,3);
%rgb_label = repmat(pixel_labels,[1 1 3]);
% here ends new

%for k = 1:nColors
 %   color = he;
 %   color(rgb_label ~= k) = 0;
 %   segmented_images{k} = color;
%end
%la imagen de interes es numero 2


[BWaa,maskedaa] = createMaskpk2(he);
Im1=imfill(BWaa,'holes');
se = strel('disk',5);
closeBW = imclose(Im1,se);
closeBW=imcomplement(closeBW);
 


%blue substraction


%IMGnew=IMGc;



for x=1:J
    for y=1:K
       if(BWHe(x,y)==1)
          BWDr(x,y)=0; 
       end
        
    end
end


SE = strel('cube',3);
% CORRIGE DIC 19
%erodedBW = imdilate(BWDr,SE);
erodedBW=BWDr;

C2 = bwconncomp(erodedBW);
numPixels = cellfun(@numel,C2.PixelIdxList);

while (min(numPixels)< 10*153) %%%correct this  IMPORTANT !!!!!!!!!!
    [biggest,idx] = min(numPixels);
    erodedBW(C2.PixelIdxList{idx}) = 0;
    C2 = bwconncomp(erodedBW);
    numPixels = cellfun(@numel,C2.PixelIdxList);
end

erodedBWch=imfill(erodedBW,'holes');
se = strel('sphere',2);
erodedBWch1 = imdilate(erodedBWch,se);
erodedBWch2 = imerode(erodedBWch,se);


%CORREGIR ESTA PARTE CONTIENE LOS CLUSTERS EROSIONADOS HASTA SEPARARLOS
%ENTRE SI. SE DEBERAN DILATAR INDIVIDUALMENTE PARA RECONOCER CADA REGION
se2 = strel('sphere',16);
erodedBWch4 = imerode(erodedBWch1,se2);
Cclean=bwconncomp(erodedBWch4);
sc=regionprops(Cclean, 'Area', 'BoundingBox');
Ac = cat(1, sc.Area);
bc = cat(1, sc.BoundingBox);
bbc=int16(bc);



%figure
%imshow(erodedBWch1)
%pause
%erodeBW contiene la imagen segmentada por metodo de substraccion de azul
%
%C = bwconncomp(erodedBW);
%statsh1=regionprops(C, 'Area')
%areash1 = cat(1, statsh1.Area);
%redonh1 = cat(1, statsh1.Eccentricity);
%mean(areash1)
%mean(redonh1)
 

%segment by V value  el valor de calibracion de esta binarizacion en 0.5
% Imagen segmentada esta en erodedBW2 y eBW contiene los mayores al 15%

BWvaluea = imbinarize(vImage, Vthrclust2);
BWvaluecluster=imbinarize(vImage,Vthrclust);

BWvalue2a=imcomplement(BWvaluea);
BWvalue2=imfill(BWvalue2a,'holes');
BWvalue2cluster=imcomplement(BWvaluecluster);
erodedBW2 = imerode(BWvalue2,SE);
erodedBW2cluster = imerode(BWvalue2cluster,SE);



% Primer filtro, se eliminan los elementos menores a Tam=20 pixeles
eBWin=erodedBW2;
eBWcluster=erodedBW2cluster;


C2 = bwconncomp(eBWcluster);
statsh2=regionprops(C2, 'Area');
areash2 = cat(1, statsh2.Area);
Th=Mina*cell_s;
%Th=Mina*mean(areash2);
% aqui se determina el valor promedio de las celulas
% IMPORTANTE -------------------------------------


numPixels = cellfun(@numel,C2.PixelIdxList);

while (min(numPixels)< Th)
    [biggest,idx] = min(numPixels);
    eBWcluster(C2.PixelIdxList{idx}) = 0;
    C2 = bwconncomp(eBWcluster);
    numPixels = cellfun(@numel,C2.PixelIdxList);
end


%Results
%eBW         -> Segmented Image cells without small cells T<15% pixels
%erodedBW2   -> Segmented by thresholding V value from HSV transformation
%erodedBWch1 -> Segmented DAB region adjusted for cluster minimum
%closeBW     -> Colour separation by colour filtering in HSV space using CreateMask Function

disp('First variables completed...')



 
 
%Classification table to train a learning machine algorithm
% It analyse train data and adds a classificaction for the image
%BBBox contains the bounding fot all the detected cells
%eeBW contains all cells detected

%Longi=size(Areafin)
%To=Longi(1,1);

%L=logical(eBW);
%box=regionprops(L,'Area', 'BoundingBox', 'centroid','Perimeter'); 
%areasas = cat(1, box.Area);
%Centroidc = cat(1, box.Centroid);
%bbox = cat(1, box.BoundingBox);
%Pericc = cat(1, box.Perimeter);
%ref = imread(nombre3);
bboxint=int16(bbbox);
%cuentacells=0;
%nana=imcomplement(erodedBW);
%figure
%subplot(2,1,2)
%title('Original GT image')
%imshow(ref)

vectordetected=0;  % este vector contiene el listado de indices validos de bbox, del total de elementos
[xa,ya]=size(bbbox);
for x=1:xa
     % solo procesa un tamaño definido por 14 veces una celula
           
        %Izzz = imcrop(erodedBW,bbox(x,:));
        %Izvv = imcrop(eBWcluster,bbox(x,:));
        subplot(2,1,1)
        title('Original clean image w detected cell')
        imshow(he)
        rectangle('Position', [bbbox(x,1),bbbox(x,2),bbbox(x,3),bbbox(x,4)],...
                 'EdgeColor','r','LineWidth',1 )
        Izzz = eeBW(bboxint(x,2):bboxint(x,2)+bboxint(x,4)-1, bboxint(x,1):bboxint(x,1)+bboxint(x,3)-1);
        
        
        %subplot(1,3,3)
        %title('Magnified cell')
        %imshow(Izzz)
        
        prompt = 'Indicate lymphocyte=1, Malignant=2, Normal=3? ';
           dG(x,1) = input(prompt)
        if(dG(x,1)==1)
            yfitC(x,:)='lymphocyt';
        end
        if(dG(x,1)==2)
            yfitC(x,:)='Malignant';
        end
        if (dG(x,1)==3)
            yfitC(x,:)='Normalepy';
        end
        if (dG(x,1)==4)
            yfitC(x,:)='Notselect';
        end     
               
      
end
disp('Train data generated ok...')
 Tab_allcells_kn(n,:)=table(Areafin, back_density, Bback, BaD, positive_density, Eccenfin, Roundfin, Centroidfin, Mayor_ellipsefin, Minor_ellipsefin, Orien_Anglefin, Perimeterfin,text_contrast,text_homogenety,hValue,sValue,vValue, yfitC);
 save Tab_allcells_kn Tab_allcells_kn;
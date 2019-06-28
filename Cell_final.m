%if (back_density~=0)

BWcelularity=Itrue;
 
 A=table2array(Tab_clust);
 Tclass=table(A(:,1),A(:,2),A(:,3),A(:,4),A(:,5),A(:,6),A(:,7),A(:,8),A(:,9),A(:,10),A(:,11),A(:,12),A(:,13),A(:,14),A(:,15),A(:,16),A(:,17),A(:,18),A(:,19),A(:,20),A(:,21),A(:,22),A(:,23),A(:,24),A(:,25),A(:,26),A(:,27),A(:,28),A(:,29),A(:,30),A(:,31),A(:,32),A(:,33),A(:,34),A(:,35),A(:,36),A(:,37),A(:,38),A(:,39),A(:,40),A(:,41),A(:,42),A(:,43),A(:,44),A(:,45));
 
 % la siguiente linea solo se ejecuta el 18-junio borrar despues
 writetable(Tab_clust,'Datos_tab.txt');

 
 [l a]=size(he(:,:,1));
 BWcelularity = zeros(l,a);
 load trainedClassifier_tree;
 %load trainedClassifier
 yfitcc=trainedClassifier_tree.predictFcn(Tclass);
 
 xt=length(yfitcc);
 %pause
 %subplot(1,3,1)
 %imshow(he)
 %hold on
 bbboxint=int16(bbbox);
 S1='Malignant';
 totmal=0;
 figure(2)
 hold on
 for na=1:xt
    tf = strcmp(S1,yfitcc(na,1));
    if (tf)
        
        totmal=totmal+1;
       BWcelularity = BWcelularity+bwselect(Itrue,table2array(Tab_clust(na,23)),table2array(Tab_clust(na,24)));
       
       %rectangle('Position', [bbbox(na,1),bbbox(na,2),bbbox(na,3),bbbox(na,4)],...
       %         'EdgeColor','r','LineWidth',1 )
       % imshow(BWcelularity)
       
    end
    
 end
 totmal
 figure(2)
 subplot(1,2,1)
 imshow(he)
 subplot(1,2,2)
 imshow(BWcelularity)
 pause(2)
imwrite(BWcelularity,filename3); 
se = strel('sphere',11);
BWcel = imdilate(BWcelularity,se);
figure(1)
subplot(1,4,4)
imshow(BWcel)
title('Malignant cells body');
C2 = bwconncomp(BWcel);
stat=regionprops(C2, 'Area');
are = cat(1, stat.Area);
[l a]=size(he(:,:,1));
Celllularity=sum(are)/(l*a);
%pause
%pause
%title('Malignant cells detected')
%subplot(1,3,2)
%imshow(BWcel)
%subplot(1,3,3)
%imshow(eeBW)
%pause
%end
%Celllularity=0;
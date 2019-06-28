prot=1;
prot2=1;
while(prot==1)
prot2=prot2+1;
C = bwconncomp(Itrue);
statsh3=regionprops(C, 'Area', 'BoundingBox');
Area = cat(1, statsh3.Area);
bbbox = cat(1, statsh3.BoundingBox);
bboxint=int16(bbbox);
prot=0;
if (sum(Area>3*cell_s))
    prot=1;
end
if(prot2>2)
    prot=0;
end
 % obtain the area calculation corresponding to label 'k'
 B=size(Area);
 
 for count=1:B(1,1)
     % compute the roundness metric
     
     if (Area(count)>3*cell_s)
         Iwater = Itrue(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1);
         figure
         imshow(Iwater)
         D = -bwdist(~Iwater);
         
         Ld = watershed(D);
         
         bw2 = Iwater;
         bw2(Ld == 0) = 0;
        
         
         mask = imextendedmin(D,1);
         imshowpair(Iwater,mask,'blend')
         pause
         D2 = imimposemin(D,mask);
         Ld2 = watershed(D2);
         bw3 = Iwater;
         bw3(Ld2 == 0) = 0;
         figure
         imshow(bw3)
         Itrue(bboxint(count,2):bboxint(count,2)+bboxint(count,4)-1, bboxint(count,1):bboxint(count,1)+bboxint(count,3)-1)=bw3;
         %figure
         %imshow(bw3)
         %title('Here applies watershed')
         pause
     end
 end
 clear Area;
 clear C;
 clear bbbox;
 clear bboxint;
 
 
 end
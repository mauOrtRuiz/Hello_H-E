function [logic,eBW2,L4,L] = tacut(BW,eBW1)
%cluster detection function
%input parameters are: BW binary image with DAB region detected
%input parameter: Centroid - centroids of every elemt of BW image obtained
%from regionprops function
% input parameter  eBW1 binary image with positive cells detected
% input paratemer AREA corresponds to regionpropg and is sincronized with
% centroid1 

% output L : total cell whithin cluster
% output L4 : image showing cluster
% output eBW2 cluster'c positive cells
% now will be a function for only a single cluster not matrix-matrix



eBW2=eBW1;
color1=254;  %positive cell
color2=185;  %fondo es naranja
logic=0;
% calibrating window size
%T=50X50
[JJ,KK]=size(eBW1);
clusterMatrix1=zeros(JJ,KK);
clustergrafic=zeros(JJ,KK);


     % analize every cell and clasiffies it according to its region
     %px=int16(Centroid(count,1));
     %py=int16(Centroid(count,2));
     
     %cleaning
     for x1=1:JJ
         for y1=1:KK
             if(BW(x1,y1)~=0)
             clustergrafic(x1,y1)=color2;                    
                 if (eBW1(x1,y1)~=0)
                     clustergrafic(x1,y1)=color1; 
                     logic=1;
                 end
             end
             if(BW(x1,y1)==0)
               if (eBW1(x1,y1)~=0) 
                  eBW2(x1,y1)==0;
               end
             end
         end
     end
  if(logic==0)
      clustergrafic(1,1)=color1;
  end
 L = bwlabel(eBW2);
 
     
 L4=label2rgb(clustergrafic);
 %imshow(L4)
 
 
end

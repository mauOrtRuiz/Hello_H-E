%Segment clusters bu watershed applying Marker Controlled


IIu=bwulterode(erodedBWch4);

se = strel('sphere',15);
Imark=imdilate(IIu,se);


% Gradient of binary image
gmag = imgradient(erodedBWch1);


D = bwdist(erodedBWch1);
DL = watershed(D);
bgm = DL == 0;


% these are the marked points
fgm4=Imark;


gmag2 = imimposemin(gmag, bgm | fgm4);
L = watershed(gmag2);
la=logical(L);
segCluster=logical(uint8(erodedBWch1).*uint8(la));

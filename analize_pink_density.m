

%10
%analize good images Pink colour
train1=imread('99858_Region_27_crop.tif');
[BW1, maskRGB]=createMask(train1);
den=sum(sum(BW1))/(512*512)
A4(11,55)=den





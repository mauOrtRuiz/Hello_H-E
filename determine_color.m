[BWblue1,maskedRGBImageblue1] = createMaskbb(Ib);
[BWblue2,maskedRGBImageblue2] = createMaskbb(Iba);
[BWblue3,maskedRGBImageblue3] = createMaskbb(Ip);

[i,j,s] = find(BWblue1);
den1=(length(i)/Tar);

[i,j,s] = find(BWblue2);
den2=(length(i)/Tar);

[i,j,s] = find(BWblue3);
den3=(length(i)/Tar);

V1=[den1 den2 den3];
[M, I]=max(V1);

figure
Iblue=segmented_images{I};
imshow(Iblue)
title('Blue image detected june 2019')
pause


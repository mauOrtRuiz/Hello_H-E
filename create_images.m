%programa para graficar el analisis de parametros obtenidos por mau
% load A2 matriz que contiene todos los parametros
load AAmeanf;
figure(1)
clf
A2=AAmeanf;
hold off
plot(A2(:,46),A2(:,20),'-+')
xlabel('Cellularity')
ylabel('Area ')
title('Morphology parameters vs Cellularity')
dim = [.4 .5 .4 .4];
str = 'Area ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
pause


figure(2)
plot(A2(:,46),A2(:,21),'-*')
hold on
plot(A2(:,46),A2(:,22),'-+')
hold on
xlabel('Cellularity')
ylabel('Eccentricity & Roundness')
title('Morphology parameters vs Cellularity')
dim = [.4 .5 .4 .4];
str = 'Round(+) Eccen(*) ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
pause


figure(3)
clf
plot(A2(:,46),A2(:,1),'-+')
hold on
plot(A2(:,46),A2(:,2),'-*')
plot(A2(:,46),A2(:,3),'-o')
plot(A2(:,46),A2(:,4),'-x')

xlabel('Cellularity')
ylabel('Eosin, He, cluster and background regions')
title('Regional Colour segmentation components vs Cellularity')
dim = [.42 .5 .4 .4];
str = 'Eosin(+) He(*) Cluster(o) Background(x)';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
axis([0 1 0 1])
pause



figure(4)
clf
hold on
plot(A2(:,46),A2(:,5),'-+')
plot(A2(:,46),A2(:,6),'-*')
plot(A2(:,46),A2(:,7),'-o')

xlabel('Cellularity')
ylabel('H,S,V colour Inter-cell region ')
title('Colour value for HSV in the cell neighborhood')
dim = [.48 .5 .2 .2];
str = 'H(+) S(*) V(o) ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
pause


figure(5)
clf
plot(A2(:,46),A2(:,29),'-+')
hold on
plot(A2(:,46),A2(:,30),'-*')
%plot(A2(:,27),A2(:,19),'-o')
%plot(A2(:,27),A2(:,20),'-x')

xlabel('Cellularity')
ylabel('Texture Contrast ')
title('  Texture Contrast vs Cellularity')
dim = [.42 .15 .15 .15];
str = 'Contrast 1(+) Contrast 2(*) ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
pause

figure(6)
clf
plot(A2(:,46),A2(:,31),'-+')
hold on
plot(A2(:,46),A2(:,32),'-*')
%plot(A2(:,27),A2(:,19),'-o')
%plot(A2(:,27),A2(:,20),'-x')

xlabel('Cellularity')
ylabel('Texture Homogenity ')
title('  Texture Homogenity vs Cellularity')
dim = [.42 .5 .4 .4];
str = 'Homogenity 1(+) Homogenity 2(*) ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

figure(7)
clf
plot(A2(:,46),A2(:,33),'-+')
hold on
plot(A2(:,46),A2(:,34),'-*')
plot(A2(:,46),A2(:,35),'-o')

xlabel('Cellularity')
ylabel('HSV cell value')
title(' Inside Cell HSV colour vs Cellularity')
dim = [.42 .5 .3 .3];
str = 'H(+) S(*) V(o)';
annotation('textbox',dim,'String',str,'FitBoxToText','on');


figure(8)
clf
plot(A2(:,46),A2(:,37),'-+')
hold on
plot(A2(:,46),A2(:,38)*200,'-*')
plot(A2(:,46),A2(:,39)*4000,'-o')

xlabel('Cellularity')
ylabel('Clusters parameters')
title(' Cluster region vs Cellularity')
dim = [.4 .5 .4 .4];
str = 'Total Clusters(+) Area(*) Radial Distance(o)';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
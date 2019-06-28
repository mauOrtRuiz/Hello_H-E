
A3=A2;
figure(1)
clf
plot(A3(:,46),A3(:,39),'-+')
hold on
plot(A3(:,46),A3(:,40),'-*')
plot(A3(:,46),A3(:,41),'-x')
xlabel('Cellularity')
ylabel('Global mean values of H, S, V')
title('Global image color values vs Cellularity')
dim = [.48 .12 .5 .4];
str = 'Global mean values: H(+) S(*) V(x) ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
pause

figure(2)
clf
plot(A3(:,46),A3(:,42),'-+')
hold on
plot(A3(:,46),A3(:,43),'-*')
plot(A3(:,46),A3(:,44),'-x')
%plot(A3(:,46),A3(:,45),'-o')
xlabel('Cellularity')
ylabel('Global density values of Eosin, He and Back')
title('Global image density E, H and back densities')
dim = [.3 .5 .4 .4];
str = 'E:+, H:*, Back:x, ';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
axis([0 1 0 1])
pause



figure(3)
Z2=A3(:,8:8+3);
surf(Z2,'EdgeColor','None');
title('H component from hsv Histogram vs Celullarity')

figure(4)
Z2=A3(:,12:12+3);
surf(Z2,'EdgeColor','None');
title('S component from hsv Histogram vs Celullarity')

figure(5)
Z2=A3(:,16:16+3);
surf(Z2,'EdgeColor','None');
title('V component from hsv Histogram vs Celullarity')

figure(6)
clf
plot(A3(:,46),A3(:,45),'-+')

xlabel('Cellularity')
ylabel('Eosin filter by a mask')
title('Global image density color values vs Cellularity')


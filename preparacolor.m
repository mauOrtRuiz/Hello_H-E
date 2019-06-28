
Iphsv=rgb2hsv(Ip);


% Pink image H value
meanIph=mean(mean(Iphsv(:,:,1)));
[i,j,s] = find(Iphsv(:,:,1));
denph=(length(i)/Tar);
menph=mean(s);

% Pink image S value
meanIps=mean(mean(Iphsv(:,:,2)));
[i,j,s] = find(Iphsv(:,:,2));
denps=(length(i)/Tar);
menps=mean(s);

% Back image H value
Ibahsv=rgb2hsv(Iba);
meanIbah=mean(mean(Ibahsv(:,:,1)));
[i,j,s] = find(Ibahsv(:,:,1));
denbah=(length(i)/Tar);
meanbah=mean(s);
% Back image S value
meanIbas=mean(mean(Ibahsv(:,:,2)));
[i,j,s] = find(Ibahsv(:,:,2));
denbas=(length(i)/Tar);
menbas=mean(s);

% Blue image H value
Ibhsv=rgb2hsv(Ib);
meanIbh=mean(mean(Ibhsv(:,:,1)));
[i,j,s] = find(Ibhsv(:,:,1));
denbh=(length(i)/Tar);
meanbh=mean(s);
% Blue image S value
meanIbs=mean(mean(Ibhsv(:,:,2)));
[i,j,s] = find(Ibhsv(:,:,2));
denbs=(length(i)/Tar);
menbs=mean(s);
[meanIph meanIbah  meanIbh] 
[meanIps meanIbas  meanIbs] 
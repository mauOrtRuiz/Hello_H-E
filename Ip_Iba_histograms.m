%Pink and background histograms

[ip1,jp1,sp1] = find(Iphsv(:,:,1));
[counts,centers] = hist(sp1,16);
figure(4)
subplot(1,3,1)
plot(centers,counts/length(ip1))
Chist1(4,:)=counts/length(ip1);

[ip2,jp2,sp2] = find(Iphsv(:,:,2));
[counts,centers] = hist(sp2,16);

subplot(1,3,2)
plot(centers,counts/length(ip2))
Chist2(4,:)=counts/length(ip2);

[ip3,jp3,sp3] = find(Iphsv(:,:,3));
[counts,centers] = hist(sp3,16);

subplot(1,3,3)
plot(centers,counts/length(ip3))
Chist3(4,:)=counts/length(ip3);
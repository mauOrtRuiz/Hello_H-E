meanItotH2=meanItotH;
meanItotS2=meanItotS;
meanItotV2=meanItotV;
denp2=denp;        
denb2=denb;        
denba2=denba;       
coR2=coR;         
coG2=coG;    
coB2=coB;         



globaltable(1,:)=table(meanItotH, meanItotS, meanItotV, denp, denb, denba, coR', coG', coB');
globaltable(2,:)=table(meanItotH2, meanItotS2, meanItotV2, denp2, denb2, denba2, coR2', coG2', coB2');
A = mean(table2array(globaltable));
%This will generate row by row elements of a global matrix for global
%partameters vs cellularity

kk=0.0;
glob_I(11,:)=table(A(1),A(2),A(3),A(4),A(5),A(6),A(7),A(8),A(9),A(10),A(11),A(12),A(13),A(14),A(15),A(16),A(17),A(18),A(19),A(20),A(21),A(22),A(23),A(24),A(25),A(26),A(27),A(28),A(29),A(30),A(31),A(32),A(33),A(34),A(35),A(36),A(37),A(38),A(39),A(40),A(41),A(42),A(43),A(44),A(45),A(46),A(47),A(48),A(49),A(50),A(51),A(52),A(53),A(54),kk);

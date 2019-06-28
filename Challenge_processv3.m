%This is the method to validate and process series of images1
% This works with 1035 first images
% directory
clear 
modoe=1; % if modoe==0 se trata de modo de entrenamiento
         % if mode==1  se trata de analisis automatizado
save modoe modoe;
prompt = 'Indicar si es el primer analisis o se continua  1=INICIO, 2=CONTINUA ? ';
flag3 = input(prompt);
save flag3 flag3
   if(flag3==1) %Primera vez que genera tabla
       nn=1;
       %Para procesar el primer lote de 1000 imagenes');
       filename = fullfile('train_mauricio_1.csv');
       
       % Para entrenar una serie desde celullarity =1 hasta cero
       %filename = fullfile('train_mauricio_anal.csv');
       
       % Use this line to process the full data set
       %filename=fullfile('train_mauricio_1.csv');
       
       % to train use file 2:
       %filename=fullfile('train_mauricio_2.csv');
       
       T=readtable(filename);
       
       save T T;
   else  % lee tabla que ya existe
       load csvtable_ch;
       load T;
       nn=height(csvtable_ch(:,1))+1;
   end
flag2=1;

[xa ya]=size(T);

while(flag2==1)
if(flag2==1)
    disp('sigue');
end
n=nn;
for n=nn:xa;
    if (modoe)
    filename2 = [num2str(T.slide(n)),'_',num2str(T.rid(n)),'.tif']
    filename3 = [num2str(T.slide(n)),'_BW',num2str(T.rid(n)),'.tif'];
    flag5=0;
    else
    %en modo de entrernamiento las imagenes deben llamarse 1_Region_1_Crop.tif   
    filename2 = [num2str(T.slide(n)),'_Region_',num2str(T.rid(n)),'_crop.tif'];
    filename3 = [num2str(T.slide(n)),'_BW_',num2str(T.rid(n)),'.tif'];
    end
    
    train1 = imread(filename2);
    challenge2;
    correctionsCh;
    back_def_process;
    if (flag5~=6)
      massive2;
    end
    
    load flag3;
    if (modoe)
        Cell_final;
    else
        imwrite(Itrue,filename3);
        if flag3==1 % primera vez que genera tabla de salida
            writetable(Tab_clust,'All_trained_data.txt');
        else  %lee tabla de salida existente
            Ta=readtable('All_trained_data.txt');
            writetable(Tab_clust,'tempor_trained_data.txt');
            T2=readtable('tempor_trained_data.txt');
            c1=height(Tab_clust(:,1));
            c2=height(Ta(:,1))+1;
            Ta(c2:c2+c1-1,:)=T2;
            writetable(Ta,'All_trained_data.txt');
                       
            
                 
        end
     
      %entrena;
    end
    
    if (modoe)
      csvtable_ch(n,:)=table(T.slide(n),T.rid(n),Celllularity, T.y(n));
      %n
      [T.y(n)  Celllularity]
      save csvtable_ch csvtable_ch;
    else
      csvtable_ch(n,:)=table(T.slide(n),T.rid(n),1);
      save csvtable_ch csvtable_ch;
    end
    save n n
    save nn nn
    save modoe modoe
    save flag5 flag5
    save flag3 flag3
    save flag2 flag2;
    clear 
    load T;
    [xa ya]=size(T);
    load n;
    load nn;
    load flag5
    load flag2
    load flag3
    load modoe
    nRange = nn:xa;
    load csvtable_ch
    if(flag2==1)
      disp('sigue');
    end
      if (flag5==6) % este es el caso que sale totalmente del programa
        n=xa;
        flag2=0;
        if(flag2==1)
           disp('sigue');
        end
        if(flag2==0)
           disp('fin ahora cerrar aplicacion');
           pause
        end
      end
      
 
 

end
if(n==xa)
    flag2=0;
end
end
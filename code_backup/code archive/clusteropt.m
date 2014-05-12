clear all
close all
clc

% disp([num2str(k)])
load winAllSubAllSeed_partialCor_645_session1.mat

sizesum=size(partWinAllSubAllSeed);
tmp(1:sizesum(1,2))=0;
for i=1:sizesum(1,1)
    for j=1:sizesum(1,2)
        tmp(j)=tmp(j)+partWinAllSubAllSeed(i,j);
    end
end

% center of all the points
for j=1:sizesum(1,2)
    tmp(j)=tmp(j)/sizesum(1,1);
end

% T----(the SUM distance between the points to the center)^2
T=0;
for i=1:sizesum(1,1)
    for j=1:sizesum(1,2)
        T=T+(partWinAllSubAllSeed(i,j)-tmp(j))^2;
    end
end

PG(1:16)=0;
R2(1:16)=0;     %R2 method
F(1:16)=0;      %pseudo-F statistic method

for i=2:2:16
    load(['indx_',num2str(i),'clusters.mat']);
    Wi(1:i,1:sizesum(1,2))=0;
    Numi(1:i)=0;
    for j=1:sizesum(1,1)
        for k=1:sizesum(1,2)
            Wi(indx(j,1),k)=Wi(indx(j,1),k)+partWinAllSubAllSeed(j,k);            
        end
        Numi(indx(j))=Numi(indx(j))+1;
    end
        
    %Wi--the center of the cluster I
    for j=1:i
        for k=1:sizesum(1,2)
            Wi(j,k)=Wi(j,k)/Numi(j);
        end
    end
    
    %PG--(the SUM distance btween the points to their cluster I's center)^2
    for j=1:sizesum(1,1)
        for k=1:sizesum(1,2)
            PG(i)=PG(i)+(partWinAllSubAllSeed(j,k)-Wi(indx(j),k))^2;
        end
    end
    
    %R^2 method    
    R2(i)=1-PG(i)/T;
    
    %pseudo-F statistic method
    F(i)=(T-PG(i))/(i-1)/(PG(i)/(sizesum(1,1)-i));
    
    
    clear Wi Numi Wi    
end

R2
F

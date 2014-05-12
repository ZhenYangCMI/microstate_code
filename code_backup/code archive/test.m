clear all
close all
clc

for i=1:10
    for j=1:20
        for k=1:6
           pop_temp(i,j,k)= (i+k)/(i+j+k);
        end
    end
end

for k=1:6
    for i=1:10
        for j=1:20
            pop(i,j)=pop_temp(i,j,k);
        end
    end
    [wList, thetaList, lambdaList, errors] = GraphicalLassoPath(pop)
    k
end
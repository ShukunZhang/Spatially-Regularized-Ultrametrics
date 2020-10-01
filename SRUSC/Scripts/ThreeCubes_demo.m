
close all;
clear all;
num=20;
X=zeros(num^3,3);

for i=1:num
    for j=1:num
        for k=1:num
            X((i-1)*num^2+(j-1)*num+k,1) = i/num;
            X((i-1)*num^2+(j-1)*num+k,2) = j/num;
            X((i-1)*num^2+(j-1)*num+k,3) = k/num;
        end
    end
end
X=vertcat(X,X,X);
scatter3(X(:,1),X(:,2),X(:,3),[],'filled');


X([num^3+1:2*(num^3)],3) = X([num^3+1:2*(num^3)],3)+1.5;
X([2*num^3+1:3* num^3],3) = X([2*num^3+1:3* num^3],3)+3;

scatter3(X(:,1),X(:,2),X(:,3),[],'filled');
X=horzcat(X,zeros(3*num^3,196));
 
 [Q,R] = qr(randn(199));
 X=X*Q;
 X=horzcat(X,ones(3*num^3,1));
 for i=1:3
     for j=1:num^3
         X((i-1)*num^3 + j,200) = i/10; 
     end
 end
  
  C_1=randi([num^3/2 num^3/2+100],10,1);
  C_2=randi([2*num^3+num^3/2 2*num^3+num^3/2+100],10,1);
  
   temp =  X(C_1,:);
 
X(C_1,:) = X(C_2,:);
   X(C_2,:) = temp;


addpath(genpath('../../SRUSC'));
load('TC_GT_100x240.mat');

SetDefaultParameters 
DenoisingOpts.Method='None';
 
SRUSCopts.KNN = 9500; 
SpectralOpts.SigmaScaling = 'Manual';
SpectralOpts.NumEig = 10;
SpectralOpts.SigmaValues = 1.12;



SpatialReg.UseReg = 1;
SpatialReg.Width = 100;
SpatialReg.Height = 240;

SpatialReg.r=95;
MajorV.Use = 0;

 
ComparisonOpts.RunEucSC = 1;
ComparisonOpts.Kmeans = 0;
ComparisonOpts.EucSCSigmaScaling = 'Manual'; 
ComparisonOpts.EucSCSigmaValues=1.12;

GeneralScript_SRUSC

figure; imagesc(reshape(Labels_SRUSC_FullData,100,240));


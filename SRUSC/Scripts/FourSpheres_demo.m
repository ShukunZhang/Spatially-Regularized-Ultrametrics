clear all;
n=4900; % number of points that you want

 radius = 1.7; % radius of the circle
 angle1 = 2*pi*rand(n,1);
 angle2 = 2*pi*rand(n,1);
 A=rand(n,2);
 B=rand(n,2);
 C=rand(n,2);
 D=rand(n,2);

 for i=1:99
    center1 = [1,3];
    r1 = radius*sqrt(rand(n,1));
    X_1 = r1.*cos(angle1)+ center1(1);
    Y_1 = r1.*sin(angle1)+ center1(2);

    center2 = [1,5];
    r2 = radius*sqrt(rand(n,1));
    X_2 = r2.*cos(angle1)+ center2(1);
    Y_2 = r2.*sin(angle1)+ center2(2);
 
    center3 = [1,7];
    r3 = radius*sqrt(rand(n,1));
    X_3 = r3.*cos(angle1)+ center3(1);
    Y_3 = r3.*sin(angle1)+ center3(2);
    
    
    center4 = [5,5];
    r4 = radius*sqrt(rand(n,1));
    X_4 = r4.*cos(angle2)+ center4(1);
    Y_4 = r4.*sin(angle2)+ center4(2);
    
    
    A=horzcat(A,X_1,Y_1);
    B=horzcat(B,X_2,Y_2);
    C=horzcat(C,X_3,Y_3);
    D=horzcat(D,X_4,Y_4);
    
end

X=vertcat(A,B,C,D);

 
 %%
 close all;
 addpath(genpath( '../../SRUSC'));
 load('FS_GT_140x140.mat'); %For large dataset 140 x 140, please load this 
 SetDefaultParameters
 
 DenoisingOpts.Method='None';
 SRUSCopts.KNN = 5000; %Number of nearest neighbors in underlying graph
 SpectralOpts.SigmaScaling = 'Manual';
 SpectralOpts.SigmaValues = 8;
 
 SpatialReg.UseReg = 0;
 SpatialReg.Width = 140;
 SpatialReg.Height = 140;

 SpatialReg.r=65;
 MajorV.Use = 0;

 
 ComparisonOpts.RunEucSC = 0;
 ComparisonOpts.Kmeans = 0;
 ComparisonOpts.EucSCSigmaScaling = 'Automatic'; 

 
 GeneralScript_SRUSC
 scatter(X(Idx_Retain,3),X(Idx_Retain,4),[],Labels_SRUSC_FullData,'filled')
 figure; imagesc(reshape(Labels_SRUSC_FullData,140,140));


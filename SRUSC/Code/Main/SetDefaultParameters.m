% This script sets default values for SRUSCopts, DenoisingOpts, SpectralOpts, and ComparisonOpts, structures required to run GeneralScript_SRUSC_SC.
% SRUSCopts, DenoisingOpts, and SpectralOpts are optional input to the
% SRUSC_SpectralClustering function

% SRUSCopts: structure with options for SRUSC approximation
%     SRUSCopts.ThresholdMethod: 'LogScaling' or 'PercentileScaling' (method for computing sequence of scales used for SRUSC approximation; default: 'LogScaling')
%     SRUSCopts.UseFixedNumScales: 0 or 1 (if 0, number of scales is determined by SRUSCopts.PercentileSize or SRUSCopts.LogRatio; if 1, number of scales is fixed by SRUSCopts.NumScales; default: 1)
%     SRUSCopts.NumScales: number of scales to use in SRUSC approximation (high values yield low SRUSC approximation error; default: 20)
%     SRUSCopts.FullyConnected: 0 or 1 (if 1, when KNN graph is disconnected the algorithm adds small edges until the graph becomes connected; default: 1; SRUSCopts.FullyConnected=0 is not recommended) 
%     SRUSCopts.KNN: integer specifying number of NN's for the KNN graph construction used to approximate SRUSC (default:20)
%     SRUSCopts.NumberofSamplePoints = integer specifying how many points to sample from each connected components when determining which small edges to add to a disconnected KNN graph (required if SRUSCopts.FullyConnected=1; default: 10)
% 
% DenoisingOpts: structure controlling data denoising
%     DenoisingOpts.KNN: number of NN's to use during denoising (denoising is based on SRUSC to DenoisingOpts.KNN-th NN; default: 20) 
%     DenoisingOpts.Method: 'Cutoff', 'ExamineFigure', or 'None' ('Cutoff' removes all points whose SRUSC to their DenoisingOpts.KNN nearest neighbor exceeds the cutoff specified in DenoisingOpts.Cutoff; 'Examine Figure' allows user to select the cut-off after examining a plot of KNN SRUSC's; 'None' keeps all data points; default: 'ExamineFigure')  
%     DenoisingOpts.Cutoff: denoising cutoff (required if DenoisingOpts.Method='Cutoff'; algorithm removes all points whose SRUSC to their DenoisingOpts.KNN nearest neighbor exceeds this cutoff)
%         
% SpectralOpts: structure with SRUSC spectral clustering options
%     SpectralOpts.NumEig: integer number of Laplacian eigenvalues and eigenvectors to compute (default: 30)
%     SpectralOpts.Laplacian: 'Symmetric_Laplacian' or 'RandomWalk_Laplacian' (Laplacian to use for SRUSC spectral clustering; default: 'Symmetric_Laplacian')
%     SpectralOpts.RowNormalization: 0 or 1 (if 1, normalizes rows of eigenvector matrix before clustering; default: 0)
%     SpectralOpts.SigmaScaling: 'Automatic' or 'Manual' (method for computing sigma scaling for SRUSC SC; default: 'Automatic')
%     SpectralOpts.NumReplicates: integer specifying number of Kmeans replicates to use for clustering the SRUSC spectral embedding (default: 5)        
% 
% ComparisonOpts: structure specifying comparisons with other methods
%     ComparisonOpts.RunEucSC: 0 or 1 (if 1, algorithm returns Euclidean spectral clustering results; default: 0)
%     ComparisonOpts.RunKmeans: 0 or 1 (if 1, algorithm returns Kmeans clustering results; default: 1)
%     ComparisonOpts.EucSCNumEig: integer specifying how many Laplacian eigenvalues and eigenvectors to compute (default: 30)
%     ComparisonOpts.EucSCSigmaScaling: 'Automatic' or 'Manual' (method for computing sigma scaling for Euclidean SC; default: 'Automatic')
%     ComparisonOpts.EucSCNumScales: integer specifying how many scales to use when computing Euclidean SC (default: 20; required when ComparisonOpts.EucSCSigmaScaling='Automatic')
%     ComparisonOpts.EucSCSigmaValues: vector of sigma values to use when computing Euclidean SC (required when ComparisonOpts.EucSCSigmaScaling='Manual')
%     ComparisonOpts.EucSCSparsity: 'None' or 'KNN' ('None' uses dense Laplacian, 'KNN' uses sparse Laplacian constructed from SRUSCopts.KNN distances only; default: 'None') 
%     ComparisonOpts.EucSCNumReplicates: integer specifying number of Kmeans replicates to use for clustering the Euclidean spectral embedding (default: 5)
%


%% Set SRUSCopts
SRUSCopts.ThresholdMethod = 'LogScaling'; %How to compute approximate SRUSC distances
SRUSCopts.UseFixedNumScales = 1; %Forces a fixed number of scales in SRUSC approximation
SRUSCopts.NumScales = 20; %Sets number of scales if SRUSCopts.UseFixedNumScales==1
SRUSCopts.FullyConnected = 1; %Makes the underlying K-NN graph fully connected in case K is too small
SRUSCopts.KNN = 20; %Number of nearest neighbors in underlying graph
SRUSCopts.NumberofSamplePoints = 10; %If SRUSCopts.FullyConnected==1, how many samples are used to connect disconnected components

%% Set DenoisingOpts
DenoisingOpts.KNN = 20; %Which NN to look at for denoising
DenoisingOpts.Method='ExamineFigure'; %Look at figure to denoise, or set a cutoff?

% Set SpectralOpts
SpectralOpts.NumEig = 30; %How many SRUSC Laplacian eigenvectors
SpectralOpts.Laplacian = 'Symmetric_Laplacian'; %Which Laplacian normalization
SpectralOpts.RowNormalization = 0; %if 1, normalizes rows of eigenvector matrix before clustering; if 0, no row normalization is done.
SpectralOpts.SigmaScaling = 'Automatic'; %Set the range of sigma automatically or with a specified range?
SpectralOpts.NumReplicates=5; %How many replicates to use for K-means

%% Set ComparisonOpts
ComparisonOpts.RunEucSC = 0; %0 for no, 1 for yes
ComparisonOpts.RunKmeans = 0; %0 for no, 1 for yes
ComparisonOpts.EucSCNumEig = 30; %How many Euclidean Laplacian eigenvectors
ComparisonOpts.EucSCSigmaScaling = 'Automatic'; %Set the range of sigma automatically or with a specified range?
ComparisonOpts.EucSCNumScales = 20; %How Euclidean scales to use
ComparisonOpts.EucSCSparsity = 'None'; %Use dense Laplacian;
ComparisonOpts.EucSCNumReplicates=5; %How many replicates to use for K-means

%% Set Spatial Regularity
SpatialReg.UseReg = 0;
SpatialReg.Width = 10;
SpatialReg.Height = 10;
SpatialReg.Sigma = 10;

%% Set Majority Voting 
MajorV.Use = 0;
MajorV.Radius = 1;
MajorV.VoteRate = 1;
MajorV.Width = 0;
MajorV.Height = 0;


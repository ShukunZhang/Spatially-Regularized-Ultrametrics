function [Labels, K_EstimateLLPD, Sigma_EstimateLLPD, Sigma_LLPD, EigVals_LLPD, EigVecs_LLPD] = LLPD_SpectralClustering(X,varargin)

% Input:
%
% -X: n x D data matrix
% -LLPDopts (optional): Structure of options for computing LLPD
% -DenoisingOpts (optional): Structure of options for denoising with LLPD
% -SpectralOpts (optional): Structure of options for running LLPD spectral
% clustering
%
% Output: 
%
% -Labels: Labels learned from LLPD spectral clustering
% -K_EstimateLLPD: Estimated number of clusters
% -Sigma_EstimateLLPD: Estimated optimal choice of sigma
% -Sigma_LLPD: All sigma values considered
% -EigVals_LLPD: LLPD eigenvalues
% -EigVecs_LLPD: LLPD eigenvectors

LLPDopts = varargin{1};
DenoisingOpts = varargin{2};
SpectralOpts = varargin{3};

% Set Default Parameter Values:

if isempty(LLPDopts)
    LLPDopts.ThresholdMethod = 'LogScaling';
    LLPDopts.UseFixedNumScales = 1;
    LLPDopts.NumScales = 20;
    LLPDopts.FullyConnected = 1;
    LLPDopts.KNN = 20;
    LLPDopts.NumberofSamplePoints = 10;
end

if isempty(DenoisingOpts)
    DenoisingOpts.KNN = 20;
    DenoisingOpts.Method='ExamineFigure';
end

if isempty(SpectralOpts)
    SpectralOpts.NumEig = 30;
    SpectralOpts.Laplacian = 'Symmetric_Laplacian';
    SpectralOpts.RowNormalization = 0; %if 1, normalizes rows of eigenvector matrix before clustering; if 0, no row normalization is done.
    SpectralOpts.SigmaScaling = 'Automatic';
    SpectralOpts.NumReplicates=5;
end


ComparisonOpts.RunEucSC = 0; 
ComparisonOpts.RunKmeans = 0;

GeneralScript_LLPD_SC
Labels = Labels_LLPD_SC_FullData;

end




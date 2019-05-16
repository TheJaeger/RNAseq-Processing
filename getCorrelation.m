function corrMat = getCorrelation(xlsPath,R,outPath)
% getCorrelation computes correlation of all possible gene pairs in an
% excel sheet containing RNA sequence triplicates
%--------------------------------------------------------------------------
%
% Usage:
%-------
% corrMat = corrMat(path_to_xls)
%
% Required input:
%----------------
% 1. xlsPath: path to Excel sheet
%    [Gene C_1 C_2 R1_1 R1_2 R1_3 R2_1 R2_2 R2_3],
%     where C --> condition,
%     Ri_N --> RNAseq observation N, N = 1:3 for triplicates
%
% Optional input:
%----------------
% 2. T: number of observations per gene (default: 3 for triplicates)
%
% Author: Siddhartha Dhiman
% Email: sdhiman@buffalo.edu & dhiman@musc.edu
% Created with MATLAB 2019a

%% Perform Checks
%  Check for file existence
if ~exist(xlsPath); error('Input file does not exist'); end

%  Check whether file can be opened
try
    tmp = xlsfinfo(xlsPath);
    disp('Input is a valid Excel file');
catch
    error('Excel file cannot be read. Ensure compatibility with your MATLAB version');
end

%  Check whether output exists. If not, save current working directory
if ~exist('outPath','var') || isempty(outPath)
    outPath = pwd;
end

%% Load Data File
[num,txt,raw] = xlsread(xlsPath);
[~,fn,~] = fileparts(xlsPath);

%  Perform check on number of triplicates specified. Try to use number of
%  columns divided by two. If fails, default to triplicates where R = 3.
[nN nR] = size(num);
if exist('R','var') || ~isempty(R)
    if rem(nR,R) == 0
        disp(sprintf('Specified R = %d matches Excel',R));
    else
        tmp = nR/2;
        if isreal(tmp) && rem(tmp,1)==0
            warning(sprintf('Specified R = %d does not match Excel. Determined R = %d.',R,tmp));
            R = tmp;
        else
            R = 3;
            warning(sprintf('Specified R = %d does not match Excel. Using default R = 3.',R));
        end
    end
end


header = {'Gene','G1O1','G1O2','G1O3','G2O1','G2O2','G2O3'};
geneList = txt(1:end,1);

%  Check for NaN and Remove them
idxNaN = ~any(isnan(num),2);
if numel(find(idxNaN == 0)) > 0
    warning(sprintf('Found %d genes with missing observations or replicates < T...discarding',numel(find(idxNaN == 0))));
else
    fprintf('All genes have replicates = T...data is excellent.');
end
geneList = geneList(idxNaN);
num = num(idxNaN,:);

%% Standardize Data (Entire Row)
num = zscore(num,0,2);

%% Form Independent Conditions
%  R1 --> First condition
%  R2 --> Second condition
%  S --> Starting columnwise index of replicate
%  E --> Ending columnwise index of replicate
R1_S = 1;
R1_E = R;
R2_S = R + 1;
R2_E = R * 2;

%% Find Correlation
corrMat1 = corr(transpose(num(:,R1_S:R1_E)));
corrMat2 = corr(transpose(num(:,R2_S:R2_E)));

%   Obtain Upper Triangle
corrMat1 = triu(corrMat1,1);
corrMat2 = triu(corrMat2,1);

%   Compress to Sparse Matrix
corrMat1 = sparse(corrMat1);
corrMat2 = sparse(corrMat2);






end
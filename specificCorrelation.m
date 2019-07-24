function corrMat = specificCorrelation(xlsPath,xlsTarget,varargin)
% specificCorrelation computes correlation of all possible gene pairs in an
% Excel sheet containing RNA sequence triplicates, provided that a target
% list is provided as well
%--------------------------------------------------------------------------
%
% Usage:
%-------
% corrMat = corrMat(xlsPath,R,outPath)
%
% Required input:
%----------------
% 1. xlsPath: path to Excel sheet
%    [Gene C_1, C_2, R1_1, R1_2, R1_3, R2_1, R2_2, R2_3],
%     where C --> condition,
%     Ri_N --> RNAseq observation N, N = 1:3 for triplicates
%
% Optional input:
%----------------
% 2. R: number of observations per gene (default: 3 for triplicates)
%
% 3. outPath: output directory to save CSV of pairwise correlations
%             (defaults: current working directory; same name as input file)
%
% Output:
%--------
% corrMat: table contianing gene-gene pairs and their correlations for the
% two conditions
%   [Gene_Source Gene_Target Correlation_C1 Correlation_C2]
%
% Author: Siddhartha Dhiman
% Email: sdhiman@buffalo.edu & dhiman@musc.edu
% Created with MATLAB 2019a
disp('==================================================================');
disp('                  Running specificCorrelation');
disp('    Target:');
disp(sprintf('        %s',xlsPath));
disp('==================================================================');

%% Parse Inputs
tic;
defaultR = 3;
defaultOut = pwd;

p = inputParser;
p.addRequired('xlsPath',@isstr);
p.addOptional('R',defaultR,@(x) rem(x,1)==0);
p.addOptional('outPath',defaultOut,@isstr);

parse(p,xlsPath,varargin{:});

R = p.Results.R;
outPath = p.Results.outPath;

%% Perform Checks
%  Close any open parallel loops
closeParallel;

%  Check for file existence
if ~exist(xlsPath) | ~exist(xlsTarget); error('Input file does not exist'); end

%  Check whether file can be opened
try
    tmp = xlsfinfo(xlsPath);
    disp('Input is a valid Excel file');
catch
    error('Excel file cannot be read. Ensure compatibility with your MATLAB version');
end
try
    tmp = xlsfinfo(xlsTarget);
    disp('Target is a valid Excel file');
catch
    error('Excel file cannot be read. Ensure compatibility with your MATLAB version');
end

%  Check whether output exists. If not, save current working directory
if ~exist('outPath','var') || isempty(outPath)
    outPath = pwd;
end

%% Load Data File
%   Load triplicate file
[numC,txtC,rawC] = xlsread(xlsPath);
[~,fn,~] = fileparts(xlsPath);
savePath = fullfile(outPath,fn);
mkdir(savePath);

%  Perform check on number of triplicates specified. Try to use number of
%  columns divided by two. If fails, default to triplicates where R = 3.
[nN nR] = size(numC);
if exist('R','var') || ~isempty(R)
    if rem(nR,R) == 0
        disp(sprintf('Specified R = %d matches dataset',R));
    else
        tmp = nR/2;
        if isreal(tmp) && rem(tmp,1)==0
            warning(sprintf('Specified R = %d does not match dataset. Determined R = %d.\n',R,tmp));
            R = tmp;
        else
            R = 3;
            warning(sprintf('Specified R = %d does not match dataset. Using default R = 3.\n',R));
        end
    end
end

geneList = txtC(1:end,1);

%  Check for NaN and Remove them
idxNaN = ~any(isnan(numC),2);
if numel(find(idxNaN == 0)) > 0
    disp(sprintf('Found %d gene(s) with missing observations or replicates < R...discarding\n',numel(find(idxNaN == 0))));
else
    fprintf('All genes have replicates = R...data is excellent.\n');
end
geneList = geneList(idxNaN);
numC = numC(idxNaN,:);

%% Dependent Functions
    function vector = vectorize(corrMat)
        %  Convert correlation matrix to vector by extracting values above
        %  the upper triangle and vectorizing them.
        sz = length(corrMat);
        vector = [];
        parfor vi = 1:sz
            cnt = 0;
            tmp2 = [];
            for vj = 1:sz
                if vi ~= vj & vj > vi
                    cnt = cnt + 1;
                    tmp1 = [vi,vj,corrMat(vi,vj)];
                    tmp2 = vertcat(tmp2,tmp1);
                end
            end
            vector = vertcat(vector,tmp2);
        end
    end  % vectorize end

% Open Parallel Pool
    function openParallel
        %  Queries the maximum number of logical cores avaialble and opens
        %  them for computation
        numworkers = feature('numcores');
        disp(sprintf('Parallel Processing: Found %d logical cores',numworkers));
        parObj = parpool('local',numworkers);
        disp('');
    end  % openParallel end

%  Close Parallel Pool
    function closeParallel
        %  Closes any open parallel pools
        parObj = gcp('nocreate');
        delete(parObj);
    end  % closeParallel end
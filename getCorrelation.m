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
%  Replicate values are transposed because MATLAB computes correlation
%  between columns. Each column is a set of observation, or replicates in
%  this case.
corrMat1 = corr(transpose(num(:,R1_S:R1_E)));
corrMat2 = corr(transpose(num(:,R2_S:R2_E)));

%  Convert correlation into vectors
corrMat1 = vectorize(corrMat1);
corrMat2 = vectorize(corrMat2);

%% Form Correlation Table
tmp = [corrMat1 corrMat2(:,3)];
for i = 1:length(tmp)
    corrMat{i,1} = geneList(tmp(i,1));
    corrMat{i,2} = geneList(tmp(i,2));
    corrMat{i,3} =  tmp(i,3);
    corrMat{i,4} = tmp(i,4);
end

%% Write Array
writetable(cell2table(corrMat),fullfile(outPath,strcat(fn,'.csv')),...
    'WriteVariableNames',false);
disp(sprintf('File saved as %s',fullfile(outPath,strcat(fn,'.csv'))));

    function vector = vectorize(corrMat)
        %  Convert correlation matrix to vector by extracting values above
        %  the upper triangle and vectorizing them.
        tmp = triu(corrMat,1);
        sz = length(tmp);
        cnt = 0;
        for i = 1:sz
            for j = 1:sz
                if i ~= j & j > i
                    cnt = cnt + 1;
                    vector(cnt,1) = i;
                    vector(cnt,2) = j;
                    vector(cnt,3) = tmp(i,j);
                end
            end
        end
    end  % vectorize end
end  % getCorrelation end
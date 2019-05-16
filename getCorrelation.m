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
[nN nR] = size(num);

%  Perform check on number of triplicates specified. Try to use number of
%  columns divided by two. If that fails, default to triplicates where R =
%  3
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

%% Independent Conditions
Cond1 = numStand(:,1:3);
Cond2 = numStand(:,4:6);

%% Find Correlation
corrMat1 = corr(transpose(Cond1));
corrMat2 = corr(transpose(Cond2));

%   Obtain Upper Triangle
corrMat1 = triu(corrMat1,1);
corrMat2 = triu(corrMat2,1);

%   Compress to Sparse Matrix
corrMat1 = sparse(corrMat1);
corrMat2 = sparse(corrMat2);

%% Form Matrix for CytoScape
%   For C
[G1 G2 corr1] = find(corrMat1);
G1 = cell2table(geneList(G1));
G1.Properties.VariableNames{'Var1'} = 'G1';
G2 = cell2table(geneList(G2));
G2.Properties.VariableNames{'Var1'} = 'G2';
corr1 = array2table(corr1);
corr1.Properties.VariableNames{'corr1'} = 'Corr';
data1 = [G1 G2 corr1];
data1 = rmmissing(data1);
data1 = sortrows(data1,3,'descend');

%   For CN
[G1 G2 corr2] = find(corrMat2);
G1 = cell2table(geneList(G1));
G1.Properties.VariableNames{'Var1'} = 'G1';
G2 = cell2table(geneList(G2));
G2.Properties.VariableNames{'Var1'} = 'G2';
corr2 = array2table(corr2);
corr2.Properties.VariableNames{'corr2'} = 'Corr';
data2 = [G1 G2 corr2];
data2 = rmmissing(data2);
data2 = sortrows(data2,3,'descend');

bins = 100;
gca = figure;
hold on;
histogram(table2array(corr1),bins,...
    'EdgeAlpha',0,...
    'FaceColor','r',...
    'FaceAlpha',0.5);
histogram(table2array(corr2),bins,...
    'EdgeAlpha',0,...
    'FaceColor','b',...
    'FaceAlpha',0.5);
hold off;
xlabel('Correlation');
ylabel('Frequency');
legend('C','CTK','Location','northwest');
grid on; box on; axis tight;
title('Absolute Pairwise Gene Correlations in N and NTK');
print('histogram','-djpeg','-r600');

% Take top 20,000 Correlations in each condition
data1 = data1(1:20000,:);
data2 = data2(1:20000,:);

% Save Tables as CSV
writetable(data1,'corrN.xlsx');
writetable(data2,'corrNTK.xlsx');


end
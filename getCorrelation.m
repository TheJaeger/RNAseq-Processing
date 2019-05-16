function corrMat = getCorrelation(xlsPath,varargin)
% getCorrelation computes correlation of all possible gene pairs in an
% Excel sheet containing RNA sequence triplicates
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
% corrMat:
%
% Author: Siddhartha Dhiman
% Email: sdhiman@buffalo.edu & dhiman@musc.edu
% Created with MATLAB 2019a

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
savePath = fullfile(outPath,fn);
mkdir(savePath);

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

geneList = txt(1:end,1);

%  Check for NaN and Remove them
idxNaN = ~any(isnan(num),2);
if numel(find(idxNaN == 0)) > 0
    disp(sprintf('Found %d genes with missing observations or replicates < R...discarding',numel(find(idxNaN == 0))));
else
    fprintf('All genes have replicates = R...data is excellent.');
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
openParallel;
corrMat1 = vectorize(corrMat1);
corrMat2 = vectorize(corrMat2);

%% Form Correlation Table
tmp = [corrMat1 corrMat2(:,3)];
clear num corrMat1 corrMat2
parfor i = 1:length(tmp)
    A = geneList(tmp(i,1));
    B = geneList(tmp(i,2));
    C =  tmp(i,3);
    D = tmp(i,4);
    corrMat(i,:) = [A B C D];
end
closeParallel;
disp(sprintf('Computed %d correlations out of %d possible correlations...discarded %d due to NaNs'...
    ,length(corrMat),nchoosek(nN,2),nchoosek(nN,2)-length(corrMat)));

%% Write Array
fprintf('Writing CSV...');
writetable(cell2table(corrMat),fullfile(savePath,strcat('corr-',fn,'.csv')),...
    'WriteVariableNames',false);
fprintf('done\n');
disp(sprintf('Correlations saved in %s',fullfile(savePath,strcat('corr-',fn,'.csv'))));

%% Plot Histograms
%  Plotting colors
c1 = [86,187,131]/255;   % condition 1 color
c2 = [78,173,241]/255;   % condition 2 color
c3 = [235,235,235]/255;  % background color

%  Legend titles
aLeg = strsplit(fn,'vs');

figure; fig = gcf;

set(fig,'PaperUnits','inches','PaperPosition',[0 0 12 9],...
    'InvertHardcopy','off','Color','white','Visible','off');

[nC1 eC1] = histcounts(cell2mat(corrMat(:,3)),64,...
    'Normalization','Probability');
[nC2 eC2] = histcounts(cell2mat(corrMat(:,4)),64,...
    'Normalization','Probability');
nC1 = smooth(nC1);
nC2 = smooth(nC2);

for k = 1:length(nC1)
    mC1(k) = median([eC1(k) eC1(k+1)]);
end
for k = 1:length(nC2)
    mC2(k) = median([eC2(k) eC2(k+1)]);
end
hold on;
c1Area = area(mC1,nC1,...
    'EdgeColor',c1,'LineWidth',3,...
    'FaceColor',c1,'FaceAlpha',0.55);
c2Area = area(mC2,nC2,...
    'EdgeColor',c2,'LineWidth',3,...
    'FaceColor',c2,'FaceAlpha',0.55);
hold off;
title(sprintf('Histogram of %s Correlations',fn))
xlabel('Correlation');
ylabel('% of Correlations');
grid on; box on; axis tight;
ax = gca;
        set(ax,'Color',c3,...
            'GridColor','white','GridAlpha',1,'MinorGridAlpha',0.15,...
            'fontname','helvetica','FontWeight','bold','fontsize',14);
legend(aLeg,'Location','north');
print(fullfile(savePath,['histo-',fn]),'-dpng','-r800');
disp(sprintf('Histogram saved in %s',fullfile(savePath,['histo-',fn,'.png'])));

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
    end  % openParallel end

%  Close Parallel Pool
    function closeParallel
        %  Closes any open parallel pools
        parObj = gcp('nocreate');
        delete(parObj);
    end  % closeParallel end
toc;
end  % getCorrelation end
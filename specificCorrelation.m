function specificCorrelation(xlsPath,xlsTarget,varargin)
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
disp('    Source:');
disp(sprintf('        %s',xlsPath));
disp('    Target:');
disp(sprintf('        %s',xlsTarget));
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

geneListC = txtC(1:end,1);

%  Check for NaN and Remove them
idxNaN = ~any(isnan(numC),2);
if numel(find(idxNaN == 0)) > 0
    disp(sprintf('Found %d gene(s) with missing observations or replicates < R...discarding\n',numel(find(idxNaN == 0))));
else
    fprintf('All genes have replicates = R...data is excellent.\n');
end
geneListC = geneListC(idxNaN);
numC = numC(idxNaN,:);

%% Standardize Data (Entire Row)
numC = zscore(numC,0,2);

%% Form Independent Conditions
%  R1 --> First condition
%  R2 --> Second condition
%  S --> Starting columnwise index of replicate
%  E --> Ending columnwise index of replicate
R1_S = 1;
R1_E = R;
R2_S = R + 1;
R2_E = R * 2;

%% Load Target Excel and Iterate of Number of Sheets
[~,tSheets] = xlsfinfo(xlsTarget);
nSheets = length(tSheets);
openParallel;
for idxSheet = 2:nSheets    %   Remove first index because it contains information on sheets
    [~,txtT,~] = xlsread(xlsTarget, idxSheet);
    idxTarget = contains(geneListC,txtT);
    disp(sprintf('%s: Found %d genes out of %d',...
        tSheets{idxSheet},nnz(idxTarget),length(txtT)));
    geneList = geneListC(idxTarget);
    num = numC(idxTarget,:);
    
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
    clear num corrMat1 corrMat2
    parfor i = 1:length(tmp)
        A = geneList(tmp(i,1));
        B = geneList(tmp(i,2));
        C =  tmp(i,3);
        D = tmp(i,4);
        corrMat(i,:) = [A B C D];
    end
    disp(sprintf('Computed %d correlations out of %d possible correlations...discarded %d due to NaNs'...
        ,length(corrMat),nchoosek(nN,2),nchoosek(nN,2)-length(corrMat)));
    
    %% Write Array
    fprintf('Writing CSV...');
    writetable(cell2table(corrMat),...
        fullfile(savePath,strcat('cr-',fn,'-',tSheets{idxSheet},'.csv')),...
        'WriteVariableNames',false);
    fprintf('done\n');
    disp(sprintf('Correlations saved in %s',fullfile(savePath,strcat('cr-',fn,'-',tSheets{idxSheet},'.csv'))));
    
    %% Plot Histograms
    %  Plotting colors
    c1 = [86,187,131]/255;   % condition 1 color
    c2 = [78,173,241]/255;   % condition 2 color
    c3 = [235,235,235]/255;  % background color
    
    %  Legend titles
    aLeg = strsplit(fn,'vs');
    
    figure; fig = gcf;
    
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 12 9],...
        'InvertHardcopy','off','Color',c3,'Visible','off');
    
    [nC1 eC1] = histcounts(cell2mat(corrMat(:,3)),64,...
        'Normalization','Probability');
    [nC2 eC2] = histcounts(cell2mat(corrMat(:,4)),64,...
        'Normalization','Probability');
    nC1 = smooth(nC1);
    nC2 = smooth(nC2);
    
    % Correlation zones
    all_y = [max(vertcat(nC1,nC2)) max(vertcat(nC1,nC2))];
    low_x = [-0.3 0.3];
    mid_x_neg = [-0.3 -0.80];
    mid_x_pos = [0.3 0.80];
    hig_x_neg = [-0.80 -1];
    hig_x_pos = [0.80 1];
    
    for k = 1:length(nC1)
        mC1(k) = median([eC1(k) eC1(k+1)]);
    end
    for k = 1:length(nC2)
        mC2(k) = median([eC2(k) eC2(k+1)]);
    end
    hold on;
    c1Area = plot(mC1,nC1,...
        'Color',c1,'LineWidth',3);
    c2Area = plot(mC2,nC2,...
        'Color',c2,'LineWidth',3);
    area(low_x,all_y,...
        'EdgeColor','none',...
        'FaceColor','red','FaceAlpha',0.05);
    area(mid_x_neg,all_y,...
        'EdgeColor','none',...
        'FaceColor','yellow','FaceAlpha',0.05);
    area(mid_x_pos,all_y,...
        'EdgeColor','none',...
        'FaceColor','yellow','FaceAlpha',0.05);
    area(hig_x_neg,all_y,...
        'EdgeColor','none',...
        'FaceColor','green','FaceAlpha',0.05);
    area(hig_x_pos,all_y,...
        'EdgeColor','none',...
        'FaceColor','green','FaceAlpha',0.05);
    hold off;
    title(sprintf('Histogram of %s vs %s Correlations in %s',...
        aLeg{1},aLeg{2},tSheets{idxSheet}));
    xlabel('Correlation');
    ylabel('% of Pairwise Correlations');
    grid on; box on; axis tight;
    ax = gca;
    set(ax,'Color',c3,...
        'GridColor','white','GridAlpha',1,'MinorGridAlpha',0.15,...
        'fontname','helvetica','FontWeight','bold','fontsize',14);
    legend(aLeg,'Location','southeastoutside');
    print(fullfile(savePath,['cr-hist-',fn,'-',tSheets{idxSheet},'.png']),...
        '-dpng','-r800');
    fprintf('Histogram saved in %s\n',...
        fullfile(savePath,['cr-hist-',fn,'-',tSheets{idxSheet},'.png']));
    clear idxTarget geneList corrMat idxTarget;
end
closeParallel;

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
end %specificCorrelation
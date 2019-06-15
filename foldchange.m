function changeMat = foldchange(xlsPath,varargin)
% foldchange computes the fold change of a gene in the transition from one
% one state to another, using RNASeq replicates as input
%--------------------------------------------------------------------------
%
% Usage:
%-------
% change = foldchange(xlsPath,R,outPath)
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
% 3. outPath: output directory to save CSV of fold changes
%             (defaults: current working directory; same name as input file)
%
% Output:
%--------
% changeMat: Table containing gene information and the fold change in
% transition
%
% Author: Siddhartha Dhiman
% Email: sdhiman@buffalo.edu & dhiman@musc.edu
% Created with MATLAB 2019a
disp('==================================================================');
disp('                        Running foldchange                        ');
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

geneList = txt(1:end,1);

%  Check for NaN and Remove them
idxNaN = ~any(isnan(num),2);
if numel(find(idxNaN == 0)) > 0
    disp(sprintf('Found %d gene(s) with missing observations or replicates < R...discarding\n',numel(find(idxNaN == 0))));
else
    fprintf('All genes have replicates = R...data is excellent.\n');
end
geneList = geneList(idxNaN);
num = num(idxNaN,:);

%% Form Independent Conditions
%  R1 --> First condition
%  R2 --> Second condition
%  S --> Starting columnwise index of replicate
%  E --> Ending columnwise index of replicate
R1_S = 1;
R1_E = R;
R2_S = R + 1;
R2_E = R * 2;

%% Compute Means Per Gene, Per Condition
[nN nR] = size(num);
openParallel;   % open parallel pool
parfor i = 1:nN
    tmp1 = mean(num(i,R1_S:R1_E));  % mean of C1
    tmp2 = mean(num(i,R2_S:R2_E));  % mean of C2
    fc{i} = log(tmp2/tmp1);          % log of fold change
end
closeParallel;  % close parallel pool
%% Construct Gene and FC Array
fc = fc';           % transpose fold change vector
changeMat = [geneList fc];      % concatenate gene list and fc

%% Write Array
fprintf('Writing CSV...');
writetable(cell2table(changeMat),fullfile(savePath,strcat('fc-',fn,'.csv')),...
    'WriteVariableNames',false);
fprintf('done\n');
disp(sprintf('Fold changes saved in %s',fullfile(savePath,strcat('fc-',fn,'.csv'))));

%% Plot Histograms
%  Plotting colors
c1 = [86,187,131]/255;   % condition 1 color
c3 = [235,235,235]/255;  % background color

%  Legend titles
aLeg = strsplit(fn,'vs');

figure; fig = gcf;

set(fig,'PaperUnits','inches','PaperPosition',[0 0 12 9],...
    'InvertHardcopy','off','Color',c3,'Visible','off');

[N E] = histcounts(cell2mat(changeMat(:,2)),64,...
    'Normalization','Probability');

for k = 1:length(N)
    mC(k) = median([E(k) E(k+1)]);
end

% Vertical lines (No change)
x_no = [0 0];
y_no = [0 max(N)];

% Vertical lines (+Double)
x_pos = [log(2) log(2)];
y_pos = [0 max(N)];

% Vertical lines (+Double)
x_neg = [-log(2) -log(2)];
y_neg = [0 max(N)];

% Plot
hold on;
plot(x_no,y_no,'-r','LineWidth',2);
plot(x_pos,y_pos,'--','Color',c2,'LineWidth',2);
plot(x_neg,y_pos,'--','Color',c2,'LineWidth',2);
pArea = area(mC,N,...
    'EdgeColor',c1,'LineWidth',3,...
    'FaceColor',c1,'FaceAlpha',0.55);
hold off;
title(sprintf('Histogram of Transition Fold Changes (%s --> %s)',...
    aLeg{1},aLeg{2}));
xlabel('Logarithmic Fold Change');
ylabel('% of Genes');
grid on; box on; axis tight;
ax = gca;
set(ax,'Color',c3,...
    'GridColor','white','GridAlpha',1,'MinorGridAlpha',0.15,...
    'fontname','helvetica','FontWeight','bold','fontsize',14);
print(fullfile(savePath,['fc-hist-',fn]),'-dpng','-r800');
disp(sprintf('Histogram saved in %s',fullfile(savePath,['fc-hist-',fn,'.png'])));

%% Dependent Functions
% Open Parallel Pool
    function openParallel
        %  Queries the maximum number of logical cores avaialble and opens
        %  them for computation
        numworkers = feature('numcores')/2;
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
toc;
end  % foldchange() end
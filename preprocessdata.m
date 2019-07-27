function preprocessdata(xlsPath,varargin)
% preprocessdata removes missing data and cleans rows flooded with zeros.
%--------------------------------------------------------------------------
%
% Usage:
%-------
% corrMat = preprocessdata(dir,outPath)
%
% Required input:
%----------------
% 1. xlsPath: master excel sheet where the excel sheets are stored
%
% Optional input:
%----------------
% 2. FilterList: Path to a list of predefined list for filtering.
%                Triplicates will only be prduced for these genes
%
% 3. Scaling: Logical true or false to specify whether to z-score
%             standardize the data
%
% 4. Threshold: Integer specifying how many non-zeros need to be present
%               before a row is removed. Default: 1
%
% 5. outPath: output directory to save CSV of pairwise correlations
%             (defaults: current working directory; same name as input file)
%
% Output:
%--------
%
% Author: Siddhartha Dhiman
% Email: sdhiman@buffalo.edu & dhiman@musc.edu
% Created with MATLAB 2019a

%% Parse Inputs
tic;
defaultScale = false;
defaultOut = pwd;
defaultTh = 1;

p = inputParser;
p.addRequired('xlsPath',@isstr);
p.addOptional('FilterList',@exists)
p.addOptional('Scaling',false,@islogical);
p.addOptional('Threshold',defaultTh,@(x) rem(x,1)==0);
p.addOptional('OutPath',defaultOut,@isstr);

parse(p,xlsPath,varargin{:});

outPath = p.Results.OutPath;

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

%% Print Configuration
if p.Results.Scaling
    scalingState = 'enabled';
else
    scalingState = 'disabled';
end

if exist('p.Results.FilterList','var') || ~isempty(p.Results.FilterList);
    filterState = 'enabled';
else
    filterState = 'disabled'
end

disp('==================================================================');
disp('                  Running commonGenes');
disp('    Target:');
disp(sprintf('        %s',xlsPath));
disp('    Variables:');
disp(sprintf('         Filter is %s',filterState));
disp(sprintf('         Threshold = %d',p.Results.Threshold));
disp(sprintf('         Scaling is %s',scalingState));
disp('==================================================================');

%% Load Data File
[num,txt,raw] = xlsread(xlsPath);
num(:,1) = [];
[~,fn,~] = fileparts(xlsPath);

hdr = txt(1,:);
hdr(1) = [];
geneList = txt(2:end,2);
nGenesA = numel(geneList);

%% Index Rows to Remove 
[~,keepIdx1] = rmmissing(num); % Rows with missing data
for i = 1:size(num,1)
    if nnz(num(i,:)) <= p.Results.Threshold
        keepIdx2(i,1) = false;
    else
        keepIdx2(i,1) = true;
    end
end
% Add the two indexes and make them logical
keepIdx = logical(keepIdx1 + keepIdx2);

%% Clean Gene List and Numbers and Scale
geneList = geneList(keepIdx); num = num(keepIdx,:);

if p.Results.Scaling
    num = normalize(num,2);
    savePath = fullfile(outPath,[fn,'_clean_scaled.csv']);
else
    savePath = fullfile(outPath,[fn,'_clean.csv']);
end

%% Load Filter List and Filter
if exist('p.Results.FilterList','var') || ~isempty(p.Results.FilterList)
    filterList = readcell(p.Results.FilterList);
    [~,idxFilt,idxTest] = intersect(geneList,filterList);
    disp(sprintf('Filter criteria: %d genes out of %d found',...
        numel(idxFilt), numel(filterList)));
    geneList = geneList(idxFilt);
    num = num(idxFilt,:);
end
nGenesB = numel(geneList);
    
%% Join Data in Table and Write
tab = cell2table(horzcat(geneList, num2cell(num)), 'VariableNames', hdr);
writetable(tab,savePath);
disp(sprintf('%d genes out of %d removed',...
    nGenesA - nGenesB,...
    nGenesA));
toc;
end % cleanData (main)



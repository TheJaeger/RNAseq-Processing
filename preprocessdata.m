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
% 1. Scaling: Logical true or false to specify whether to z-score
%             standardize the data
%
% 2. outPath: output directory to save CSV of pairwise correlations
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

p = inputParser;
p.addRequired('xlsPath',@isstr);
p.addOptional('Scaling',false,@islogical)
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

disp('==================================================================');
disp('                  Running commonGenes');
disp('    Target:');
disp(sprintf('        %s',xlsPath));
disp('    Variables:');
disp(sprintf('         Scaling is %s',scalingState));
disp('==================================================================');

%% Load Data File
[num,txt,raw] = xlsread(xlsPath);
num(:,1) = [];
[~,fn,~] = fileparts(xlsPath);
savePath = fullfile(outPath,[fn,'_clean.csv']);

hdr = txt(1,:);
hdr(1) = [];
geneList = txt(2:end,2);

%% Index Rows to Remove 
[~,keepIdx1] = rmmissing(num); % Rows with missing data
keepIdx2 = any(num,2);         % Rows not completely filled with zeros
% Add the two indexes and make them logical
keepIdx = logical(keepIdx1 + keepIdx2);

%% Clean Gene List and Numbers and Scale
geneList = geneList(keepIdx); num = num(keepIdx,:);

if p.Results.Scaling
    num = normalize(num,2);
end
    
%% Join Data in Table and Write
tab = cell2table(horzcat(geneList, num2cell(num)), 'VariableNames', hdr);
writetable(tab,savePath);
toc;
end % cleanData (main)



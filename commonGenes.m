function commonGenes(xlsPath,R,C,varargin)
% commonGenes computes all genes that commonly occur in all group and spits
% saves the list as CSV in working directory or specified directory
%--------------------------------------------------------------------------
%
% Usage:
%-------
% corrMat = commonGenes(dir,outPath)
%
% Required input:
%----------------
% 1. xlsPath: master excel sheet where the excel sheets are stored 
%
% Optional input:
%----------------
% 2. outPath: output directory to save CSV of pairwise correlations
%             (defaults: current working directory; same name as input file)
%
% Output:
%--------
% geneVec: vector containing list of genes
%
% Author: Siddhartha Dhiman
% Email: sdhiman@buffalo.edu & dhiman@musc.edu
% Created with MATLAB 2019a
disp('==================================================================');
disp('                  Running commonGenes');
disp('    Directory:');
disp(sprintf('        %s',dir));
disp('==================================================================');

%% Parse Inputs
tic;
defaultR = 3;
defaultC = 5;
defaultOut = pwd;

p = inputParser;
p.addRequired('xlsPath',@isstr);
p.addRequired('R',defaultR,@(x) rem(x,1)==0);
p.addRequired('C',defaultC,@(x) rem(x,1)==0);
p.addOptional('OutPath',defaultOut,@isstr);
p.addOptional('FillMissing',false,@islogical)

parse(p,xlsPath,varargin{:});

R = p.Results.R;
outPath = p.Results.outPath;

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

% Check whether R is provided
if ~exist('R','var') || isempty(R)
    error('Please provide variable R for number of replicates');
end

% Check whether C is provided
if ~exist('C','var') || isempty(C)
    error('Please provide variable C for number of conditions');
end

%% Print Configuration
disp('==================================================================');
disp('                  Running commonGenes');
disp('    Target:');
disp(sprintf('        %s',xlsPath));
disp('    Variables:');
disp(sprintf('         R = %d | C = %d',R,C));
disp('==================================================================');

%% Load Data File
[num,txt,~] = xlsread(xlsPath);
[~,fn,~] = fileparts(xlsPath);
savePath = fullfile(outPath,'common_genes.csv');

geneList = txt(2:end,2);

if FillMissing
    
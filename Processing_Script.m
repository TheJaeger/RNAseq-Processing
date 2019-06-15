clc; clear all; close all;

%% Define Paths for Input and Output
inPath = '/Users/sid/Downloads/Data-selected';
outPath = '/Users/sid/Downloads/Correlations';

%% Start
dirIn = dir(fullfile(inPath,'*.xlsx'));

for i = 1:length(dirIn)
    path2File = fullfile(dirIn(i).folder,dirIn(i).name);
    getCorrelation(path2File,3,outPath);
end
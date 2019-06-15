clc; clear all; close all;

%% Define Paths for Input and Output
inPath = 'D:\Datasets\RNAS-Seq\mRNA';
outPath = 'D:\SystemFiles\siddh\Downloads\RNA-Seq'
R = 3;

%% Start
dirIn = dir(fullfile(inPath,'*.xlsx'));

for i = 1:length(dirIn)
    path2File = fullfile(dirIn(i).folder,dirIn(i).name);
    getCorrelation(path2File,R,outPath);
    foldchange(path2File,R,outPath);
end
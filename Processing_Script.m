clc; clear all; close all;

%% Define Paths for Input and Output
inPath = 'D:\Datasets\RNAS-Seq\mRNA\Triplicates';
targetPath = 'D:\Datasets\RNAS-Seq\mRNA\Targets';
outPath = 'D:\SystemFiles\siddh\Downloads\RNA-Seq';
R = 3;

%% Start
dirIn = dir(fullfile(inPath,'*.xlsx'));
dirTarget = dir(fullfile(targetPath,'*xlsx'));

for i = 1:length(dirIn)
    path2File = fullfile(dirIn(i).folder,dirIn(i).name);
    getCorrelation(path2File,R,outPath);
    foldchange(path2File,R,outPath);
end

for i = 1:length(dirTarget)
    tgt = strsplit(dirTarget(i).name,'_');
    source = fullfile(dirIn(1).folder, tgt{2});
    specificCorrelation(source,...
        fullfile(dirTarget(i).folder,dirTarget(i).name),R,outPath);
end
    
clc; clear all; close all;

%% Define Paths for Input and Output
masterPath = 'D:\Datasets\RNAS-Seq\mRNA_master_list.xlsx';
filterPath = 'D:\Datasets\RNAS-Seq\mRNA\Targets_Cleaned\All_Targets.xlsx';
masterOutPath = 'D:\Datasets\RNAS-Seq\mRNA\Preprocessed_Master_List';
inPath = 'D:\Datasets\RNAS-Seq\mRNA\Triplicates_Cleaned';
targetPath = 'D:\Datasets\RNAS-Seq\mRNA\Targets';
outPath = 'D:\Datasets\RNAS-Seq\mRNA\Generated_Data';
R = 3;

%% Preprocess Data
preprocessdata(masterPath,...
    'FilterList', filterPath,...
    'Scaling', true,...
    'Threshold', 8,...
    'OutputPath', masterOutPath);

%% Start
dirIn = dir(fullfile(inPath,'*.xlsx'));
dirTarget = dir(fullfile(targetPath,'*xlsx'));

for i = 1:length(dirIn)
    path2File = fullfile(dirIn(i).folder,dirIn(i).name);
    getCorrelation(path2File,R,outPath);
    %foldchange(path2File,R,outPath);
end

for i = 1:length(dirTarget)
    tgt = strsplit(dirTarget(i).name,'_');
    source = fullfile(dirIn(1).folder, tgt{2});
    specificCorrelation(source,...
        fullfile(dirTarget(i).folder,dirTarget(i).name),R,outPath);
end
    
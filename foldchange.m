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
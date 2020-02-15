# RNAseq-Processing
A pipeline developed to process mRNA or miRNA data derived from RNA sequencing.

## Steps Performed
The following steps were undertaken to generate all data:

1. Run `preprocessdata` on __mRNA_master_list.xlsx__ to generate a cleaned master list named __mRNA_master_list_clean_scaled.xlsx__ using __All_Targets.xlsx__. This produces a master list that is scaled evenly across all conditions (C, N, +N, -TK), and includes at least 8 non-zeros for any given gene.

    '''
    masterPath = 'D:\Datasets\RNAS-Seq\mRNA_master_list.xlsx';
    filterPath = 'D:\Datasets\RNAS-Seq\mRNA\Targets_Cleaned\All_Targets.xlsx';
    masterOutPath = 'D:\Datasets\RNAS-Seq\mRNA\Preprocessed_Master_List';

    preprocessdata(masterPath,...
    'FilterList', filterPath,...
    'Scaling', true,...
    'Threshold', 8,...
    'OutputPath', masterOutPath);
    '''

2. Open __mRNA_master_list_clean_scaled.xlsx__ and create pair-wise observations in the following categories:

    * C vs N
    * C vs CTK
    * N vs NTK
    * N vs NNLS

    Then place them in another folder called __Triplicates_Cleaned__

3. Then run `getCorrelation` to generate pairwise correlation on groupd in folder __Triplicates_Cleaned__

4. Run ``specificCorrelation` to generate pairwise correlations and figures for GO group, using targets in __Targets__

The file `Processing_Script.m` contains all these processing steps highlighted above.
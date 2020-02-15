# RNAseq-Processing
A pipeline developed to process mRNA or miRNA data derived from RNA sequencing.

## Steps Performed
The following steps were undertaken to generate all data:

1. Run `preprocessdata` on _mRNA_master_list.xlsx_ to generate a cleaned master list named _mRNA_master_list_clean_scaled.xlsx_ using _All_Targets.xlsx_. This produces a master list that is scaled evenly across all conditions (C, N, +N, -TK), and includes at least 8 non-zeros for any given gene.

    ```
    masterPath = 'D:\Datasets\RNAS-Seq\mRNA_master_list.xlsx';
    filterPath = 'D:\Datasets\RNAS-Seq\mRNA\Targets_Cleaned\All_Targets.xlsx';
    masterOutPath = 'D:\Datasets\RNAS-Seq\mRNA\Preprocessed_Master_List';

    preprocessdata(masterPath,...
    'FilterList', filterPath,...
    'Scaling', true,...
    'Threshold', 8,...
    'OutputPath', masterOutPath);
    ```

2. Open _mRNA_master_list_clean_scaled.xlsx_ and create pair-wise observations in the following categories:

    * C vs N
    * C vs CTK
    * N vs NTK
    * N vs NNLS

    Then place them in another folder called _Triplicates_Cleaned_

3. Then run `getCorrelation` to generate pairwise correlation on groupd in folder _Triplicates_Cleaned_

4. Run `specificCorrelation` to generate pairwise correlations and figures for GO group, using targets in _Targets_

The file `Processing_Script.m` contains all these processing steps highlighted above.

we have developed a method which is based on read-level to calculate cell type-specific methylation regions, which can be used as features for dissecting biological functions and identifying cell types.The process is divided into four steps:

![ok3j3x3q](https://github.com/user-attachments/assets/a4432e49-5577-4be0-b139-acb9f5692698)

First, segment the whole genome. This step uses a dynamic programming segmentation algorithm to divide the genome into consistent region blocks with the software named wgbstools (https://github.com/nloyfer/wgbs_tools) developed by Netanel et al[1]. 
We use *segment* command to divide the whole genome.The input is beta files of all your samples, block.bed is the output. The region less than 4 CpGs will be discarded.

```
./wgbstools segment --betas ../cell_type_marker/beta/GSM5652183_Saphenous-Vein-Endothel-Z000000SB.hg38.beta \\
../cell_type_marker/beta/GSM5652182_Saphenous-Vein-Endothel-Z000000S7.hg38.beta \\
../cell_type_marker/beta/GSM5652181_Saphenous-Vein-Endothel-Z000000RM.hg38.beta \\
 --min_cpg 4 -o ../cell_type_marker/uxm_marker/plasma_marker/block.bed -@ 50
```
Second, calculate the methylation value at read-level of each sample. The value is defined as alpha value.

```
python get_reads_alpha_value.py path_to_wgbstools path_to_block_file path_to_sample_pat path_to_output_bed_file

example:
python get_reads_alpha_value.py /mnt/data3/methylation/wgbs_tools /mnt/data3/block.bed /mnt/data3/pat_beta/GSM5652303_Blood-Monocytes-Z000000U3.hg38.pat.gz /mnt/data3/2_alpha_bed/GSM5652303_Blood-Monocytes.bed
```

Third, calculate the block region mean methylation level based on alpha value of every reads.
```
python region_mean_alpha.py path_to_blcok_file path_to_output_bed_file threads path_to_output_alpha_csv

example:
python region_mean_alpha.py /mnt/data3/block.bed /mnt/data3/2_alpha_bed/GSM5652237_Liver-Hepatocytes_alpha.bed 50 /mnt/data3/3_region_mean_alpha/GSM5652237_Liver-Hepatocytes_alpha.csv
```

Finally, test the difference of mean methylation levels between two group samples in every block region with Wilcoxon rank-sum test. Besides, the delta mean of each region will be calculated.
```
python test_region.py path_to_alpha_csv path_to_output_marker_file threads all_cell/tissue_types

example:
python test_region.py /mnt/data3/3_region_mean_alpha/ /mnt/data3/output_marker/ 30 ['Vein-Endothel','Hepatocytes','Erythrocyte','Monocytes','Granulocytes']
```
After we get the marker file for each cell type, we use p value < 0.05 and delta mean > 0.5 to filter out the regions with differences of every cell type.

*[1] Loyfer N, Rosenski J, Kaplan T. wgbstools: A computational suite for DNA methylation sequencing data representation, visualization, and analysis, bioRxiv 2024:2024.2005.2008.593132*

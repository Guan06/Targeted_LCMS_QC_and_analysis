Please follow the instruction step-by-step for a bug-free pipeline experience.

### Data reformat  

Before start, make sure that the input file here (i.e. output file from MassHunter) is the same format as the example given here.

Specifically, the number of meta data columns (i.e. number of columns before columns contain measuring results of each compound) should be the same, as well as the order of columns for each compound (in this case: FWHM, Area, RT, S/N, Calc. Conc.).

Changes in these formatting would affect the script `step0_data_reformat.R`.


1) Get the first row from the ../input/*_quant.csv (will be referred as data file below), 
paste to _compound_raw.lst_ and reformat to _compound.lst_ with: 

`cat compound_raw.lst |sed 's/\t/\n/g' |sed '/^$/d' |sed 's/ Results//' | sed '1d' >compound.lst`

Add "#" to the first row of data file, so that in later steps R can ignore reading in this header row.

__ACHTUNG__: If the first row is "Sample", remove the first row (double check before removing).

Make a soft link of note file from __input__ folder.

`ln -s ../inputs/notes.csv`

2) Get the compound list with mapping information and tidy version with script `step0_data_reformat.R`

`Rscript step0_data_reformat.R compound.lst ../input/standard_concentrations_5xdil.csv ../input/2023*_quant.csv`

### Limit of Detection (LOD) 

3) Plot the standard curve and detect the LOD.

`Rscript step1_QC.R`

### Data normalization 
4) Get the reformatted data table with LOD and normalize the concentration by internal standards as well as using pqn.

`Rscript step2_get_sample_lod.R`

Internal standards were added to each sample (except QC samples)

### Results

5) copy the files to the Rproject folder.

```{bash}
cp s1_rt.pdf ../report/data
cp s1_std_curve.pdf ../report/data
cp s1_data_QC.txt ../report/data
cp s2_reformat_data_with_LOD_reNorm.txt ../report/data
cp notes.csv ../report/data
cp s0_reformat_std_concentration.txt ../report/data
cp s3_area_over_time_facet.pdf ../report/data
cp s3_QC_samples_area_over_time_SD_Mean.txt ../report/data
```
# phenotyper
functions to extract and analyze phenotyper data

## procedure to prepare data files ready for processing
- export .txt data files from Ethovision software

- each .txt data file is originally encoded with {unicode}

  open each .txt file and overwrite with {UTF-8} encoding
  
  added advantage is that each .txt file now takes up half of memory size (from approx 1GB to 500MB)
  
- optionally: to further compress the memory load, put all .txt data files in a zipped folder (e.g. 4 zipped {UTF-8} encoded .txt data files will take up approx 150MB memory space instead of 2GB)
 
  the R script is able to unzip the folder content and process the data files from there
  
- place all .txt data files (zipped or unzipped) as well as the .xlsx meta file together in one folder

## dummy .txt data files
- the {UTF-8} encoded .txt files are saved as .RDS files

- use function unload_dummies() to get these data files for practise purposes

## metafile
use function download_meta() to get the template .xlsx file for adding mouse information

## function content:
### helper functions
color_spectrum()

bar_spacing()

### template and example files
import_dummies()

download_meta()

### process functions
import_raw_cw()

cw_entries()

cw_summary()

survival_data()

survival_stat()

survival_summary()

entries_data()

entries_stat()

entries_summary()

### plot functions
accuracy_plot()

survival_plot()

multi_survival_plot()

entries_plot()

time_plot()

### update functions
new_threshold()

new_genotype()

merge_list()

filter_list()

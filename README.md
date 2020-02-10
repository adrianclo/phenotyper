# PHENOTYPER COGNITION WALL (CW)
https://www.noldus.com/phenotyper/add-ons

functions to extract and analyze phenotyper data

## procedure to prepare data files ready for processing
- export .txt data files from Ethovision software as {ANSI} encoded

- In case you have exported each .txt data file as {unicode}, open each .txt file and overwrite with {UTF-8} encoding
  
- optionally: to further compress the memory load, put all .txt data files in a zipped folder 

  e.g. four zipped {UTF-8} encoded .txt data files will take up approx 150MB memory space instead of 2GB
 
  the R script is able to unzip the folder content and process the data files from there
  
- place all .txt data files (zipped or unzipped) as well as the .xlsx meta file together in one folder

## dummy .txt data files
- the {UTF-8} encoded .txt files are saved as .RDS files

- use function cw_dummies() to get these data files for practise purposes

## metafile
use function cw_meta() to get the template .xlsx file for adding mouse information

## cw function content:
### helper functions
color_spectrum()

bar_spacing()

### template and example files
cw_dummies()

cw_meta()

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

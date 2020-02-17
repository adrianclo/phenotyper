# PHENOTYPER COGNITION WALL (CW)
https://www.noldus.com/phenotyper/add-ons

functions to extract and analyze phenotyper data

## procedure to prepare data files ready for processing
- export .txt data files from Ethovision software as {ANSI} encoded

- in case you exported each .txt data file as {unicode}, open each .txt file separately and overwrite with {UTF-8} encoding
  
  **be aware** that opening the files takes time, and some files cannot be opened with the standard text editor (max read capacity exceeded). for this, open the file with editpad lite7 {https://www.editpadlite.com/} and copy chunk by chunk into the standard text editor
  
- **optionally**: to compress the memory load, put all .txt data files in a zipped folder 

  e.g. four zipped {ANSI/UTF-8} encoded .txt data files will take up approx 150MB memory space instead of 2GB
 
  the R script is able to unzip the folder content and process the data files from there (and will delete the unzipped files later)
  
- place all .txt data files (zipped or unzipped) as well as the .xlsx meta file together in one folder

- processed data can be stored as .RDS files to even more compress the memory load
    * option in the import function
    
    * otherwise use saveRDS() manually! 

    * storage less than 1MB per experiment!
  
## dummy .txt data files
- the {ANSI/UTF-8} encoded .txt files are saved as .RDS files

- use function cw_dummies() to get these data files for practise purposes

## metafile
use function cw_meta() to get the template .xlsx file for adding mouse information

## cw function content:
### helper functions
color_spectrum()

gg_color_hue()

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

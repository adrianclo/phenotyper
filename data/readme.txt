# Procedure to obtain data files ready for processing

- Export .txt data files from Ethovision
- Each .txt data file is originally encoded with {unicode}. 
  Open each .txt data file and overwrite with {UTF-8} encoding.
- Added advantage is that each .txt data file takes up half memory size (from approx 1 GB to 500 MB).
- To further compress the memory load: put these .txt data files in a zipped folder.  
  The R script is able to unzip the folder and process the data files from there.

# Repository contains dummy .txt data files

- The {UTF-8} encoded .txt data files are saved as .RDS files
- Use function unload_dummies() to get these dummy .txt data files
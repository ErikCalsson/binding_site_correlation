## Core part: Binding site correlation python tool
The scope of this project is to write a small tool to analyze the overlap between two BED files (BED6, i.e. tab separated CHROM, START, END, NAME, SCORE, STRAND) and return a statistic whether both files are dissimilar and how significant this dissimilarity is.

### Input

Two BED6 formatted files.

### Output

Statisticts about the overlap of the BED files to console -- or to a given file (optional)

## Visualization part: Intergration of interactive elements (widgets) using Dash

### Input

Two BED6 formatted files.

### Output

Statisticts about the overlap of the BED files to an interactive and browseable web-application using plotly Dash. 


# Installation instructions

Download the current version:
```
git clone https://github.com/ErikCalsson/binding_site_correlation.git
```

Install requirements via pip3:
```
pip3 install -r requirements.txt
```

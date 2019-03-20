run-umap.py
===========
The python code to perform dimensionality reduction has the following package dependencies:
* numpy
* pandas
* umap
* sklearn

transcriptomoscope.R
====================
The R code to plot spatial transcriptomics data using Voronoi tessellation has these dependencies:
* deldir
* sp
* rgeos
* alphahull
* RColorBrewer

How to use
==========
Please refer to the command line help available with -h for both scripts, and to the following usage examples.

Usage examples
==============

Simply plot the relative frequencies:
```
~/code/transcriptomoscope/transcriptomoscope.R counts*.tsv.gz
```

Plot the relative frequencies, allowing the plot coordinates to range from 0 to 1000 in both dimensions:
```
~/code/transcriptomoscope/transcriptomoscope.R counts*.tsv.gz -C 10000 -R 10000
```
This is useful when the coordinates of the spots are not in the standard range of 33 columns and 35 rows.


Perform dimensionality reduction with t-SNE and subsequently plot the results:
```
~/code/transcriptomoscope/run-umap.py --dim-red t-SNE  -o t-SNE_ counts*.tsv.gz
~/code/transcriptomoscope/transcriptomoscope.R --one t-SNE_joint_*txt -o t-SNE.pdf
```

Perform dimensionality reduction with UMAP and subsequently plot the results:
```
~/code/transcriptomoscope/run-umap.py --dim-red UMAP  -o umap_ counts*.tsv.gz
~/code/transcriptomoscope/transcriptomoscope.R --one umap_joint_*txt -o umap.pdf
```

Perform dimensionality reduction with UMAP and subsequently plot the results with non-default coordinates:
```
~/code/transcriptomoscope/run-umap.py --dim-red UMAP  -o umap_ counts*.tsv.gz
~/code/transcriptomoscope/transcriptomoscope.R --one umap_joint_*txt -o umap.pdf -C 10000 -R 10000
```

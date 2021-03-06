![logo](https://github.com/duxfrederic/ddmsoft/blob/main/logo/logo.png) 

DDMSoft was developped during an internship at the RWTH Aachen in the summer of 2019. My stay was organized by Jerome J. Crassous and funded by a [BioSoft](http://www.ihrs-biosoft.de/ihrs-biosoft/EN/GSP/GSP_node.html) scholarship. For an introduction and examples see the associated github page: https://duxfrederic.github.io/ddmsoft/.
## What is DDMSoft?
DDMSoft is a Python based GUI that essentially converts a differential dynamic microscopy (DDM) video to a DDM matrix.  Added to that are fitting routines and tools that make the handling of the massive amount of data generated by a DDM analysis easier. DDMSoft is best used on the computer that captures the videos from the microscope because it allows for an instantaneous feedback on the measurement. Here are some of the features implemented in DDMSoft:
 - fast conversion of your videos to DDM matrices
 - easy exploration and selection of the matrices to be further analyzed
 - quick evaluation of the quality of the measurement with a wavenumber dependent plot of correlation function
 - multi-wavenumber fit of the matrix with 
    * single exponential decays (with and without stretch)
    * double exponential decays (with and without stretch)
    * cumulants
    * exponential and flow 
    * CONTIN
- quick with sizing, see <https://youtu.be/YvFRs4CWKe8>.



## Installation
For Windows users, a 'compiled' version is available at <https://nc.fdux.ch/index.php/s/RBpG7TNHM8Eg648>.


For Linux users, the best is to download Anaconda: <https://www.anaconda.com/products/individual>. Anaconda comes with fast algebra libraries and thus will be quick at processing the data. Then:
```bash
$ git clone git@github.com:duxfrederic/dddmsoft.git
$ cd ddmsoft
$ conda env create --file ddmsoft.yml
$ conda activate ddmsoft
$ python DDMsoft.py
```


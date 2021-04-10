CRISPR Specificity Correction (CSC)

SYSTEM REQUIREMENTS

This package currently requires python3 to install and run. This software was developed on a macintosh computer using the Catalina and Big Sur Operating Systems. However, the software is operating system independent and will run on machines of various build and manufacture.

INSTALLATION GUIDE

This package requires scikit-learn and sklearn-contrib-py-earth 0.1.0

CSC will attempt to install pyearth upon running setup.py, however, it is strongly recommended that the user run the following command prior to running the setup script.
    
    system
        conda install sklearn-contrib-py-earth
        
Additional install commands for pyearth are found here:

https://anaconda.org/conda-forge/sklearn-contrib-py-earth

To install system wide, run 

    system
        python setup.py build
        python setup.py install

Alternatively a user can install via pypi with

     system
        pip install CSC-crispr
        
For local installation, run something like

    system
        python setup.py install --user

and then make sure that the local directory with binaries (such as `$HOME/Library/Python/3.8/bin/`) is available in your PATH.

Typical Installation time is less than 5 minutes.

To run CSC from the command line, call the entrypoint function

    system
        csc_process -h

Alternatively CSC can be run from inside the package directory with the following command

    system
        python csc_lite.py -h

CSC can take multiple input file formats. CSC reads in a file with the expected format where the first column contains the 20mer gRNA sequences and the second column is a numerical metric (ie: logFC or other output metric). The CSC is able to correct logFC of gRNAs both for pre-computed human genome-wide libraries (Avana, Brunello, GeckoV1, GeckoV2, TKOV3) and custom libraries even if they do not have a pre-computed file. 

EXAMPLE DEMONSTRATION

If the user utilizes as input a file with gRNA and logFC and the gRNA belongs to a pre-computed library then logFC adjusted for off-target effects will be computed with the following command

    system
        csc_process -i filepath-to-input-file -l name-of-precomputed-library

A concrete example of this functionality is available in the /csc_v2/screen_models/examples/

    system
        csc_process -i example_grna_logfc -l example

If the user utilizes a file containing only gRNA sequences and the gRNA belongs to a pre-computed library then the off-target space up to a Hamming distance of 3 and the GuideScan specificity score are rendered for each gRNA

    system
        csc_process -i example_grna -l example

If a user desires to make their own custom CSC correction for an individual library then they may utilize CSC hg38 (and soon mm10) genome-wide pickle files. Genome wide pickles are located in the following repository.

https://figshare.com/s/942320ff4d38cb93b39a

https://drive.google.com/drive/folders/1U9H3r_CEOpa-MULLCN3awwOdOqG4101T?usp=sharing

    system
        csc_process -i absolute-filepath-to-input-file -g absolute-filepath-to-genome-wide-pickle

If the custom file contains both gRNA sequences and numerical metric (ie: logFC or other output metric) values in the first two columns then CSC adjustments will be modeled and made. If only gRNA sequences are present then the output will be the off-target space up to a Hamming distance of 3 and the GuideScan specificity score for each gRNA

To uninstall CSC execute the following command
   
    system
        pip uninstall csc

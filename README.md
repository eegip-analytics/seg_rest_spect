# seg_rest_spect

Tools for segmentation of resting data and spectral power calculation of Lossless derivative EEGIP data sets.

Copy this 'seg_rest_spect' folder to the 'eegip/[sitename]/derivatives/lossless/derivatives' folder. 

Use EEGLAB with Matlab working directory set to:

eegip/[sitename]/derivatives/lossless

Add to the Matlab path:

addpath('code/dependencies/eeglab_asr_amica/');

Then start EEGLAB:

eeglab

To run the segmentation use the batch_context EEGLAB extension to execute the "seg_rest_[site_ses].htb" script (located in 'derivatives/seg_rest_spect/code/scripts') on the *.qcr.set files in the lossless derivatives directory.

To extract the power values to a csv file run the getSpectPow.m function as follows:

getSpectPow('derivatives/seg_rest_spect/code/misc/fnames.txt','derivatives/seg_rest_spect/output.csv');

... where the fnames.txt file contains the segmented file names to include in the analysis. generated as follows from a terminal:
 
find . -type f -name "*-SEGrest*.set" > derivatives/seg_rest_spect/code/misc/fnames.tx

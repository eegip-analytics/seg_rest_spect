# seg_rest_spect

Tools for segmentation of resting data and spectral power calculation of Lossless derivative EEGIP data sets.

Copy this 'seg_rest_spect' folder to the 'eegip/[sitename]/derivatives/lossless/derivatives' folder. 

Use EEGLAB with Matlab working directory set to:

eegip/[sitename]/derivatives/lossless

Add to the Matlab path:

addpath('code/dependencies/eeglab_asr_amica/');

Then start EEGLAB:

eeglab

To run the segmentation use the batch_context EEGLAB extension to execute the "seg_rest_spect.htb" script (located in 'derivatives/seg_rest_spect/code/scripts') on the *.qcr.set files in the lossless derivatives directory.


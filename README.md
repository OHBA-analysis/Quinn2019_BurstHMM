# Unpacking the transient event dynamics that underlie spontaneous changes and induced responses in electrophysiology
Scripts for the analysis and tutorial in *Quinn et al (2019)*

## Requirements

The analysis requires the following software running on a Unix-Type operating system.
 - MatLab 2018a or greater (may run on earlier versions but not tested)
 - MatLab Signal Processing Toolbox
 - MatLab Wavelet Toolbox
 - distributionPlot
	 -  [https://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m](https://uk.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m)

  ## Getting started
 1. Download or clone this repository to your computer
 2. Ensure that you have a recent MatLab and access to the Signal Processing Toolbox and Wavelet Toolbox
 3. Edit the file paths in the top of ```hmm_0_initialise``` to point to the location of these toolboxes on your computer
 4. run ```hmm_0_initialise``` in MatLab, if this returns without error then you are good to go. If you see warnings then some dependencies may be missing, follow the instructions in the warning message.
 5. Work through the hmm tutorial scripts

## Contents

```hmm_0_initialise``` runs the initial setup and configuration for these analyses

```hmm_1_dynamics_illustration``` creates and plots the power simulations from figure 1

```hmm_2_envelope``` runs an amplitude-envelope HMM on simulated data and creates figure 2

```hmm_3_embedded``` runs an time-delay embedded HMM on simulated data and creates figure 3

```hmm_4_realdata_trialwise``` runs an time-delay embedded HMM and task-evoked analysis on source-space MEG data and creates figures 4,5 and 6
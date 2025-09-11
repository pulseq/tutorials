# Pulseq tutorial "Labels, 2D Mode, ICE integration and GRAPPA"

Welcome to the "Labels, 2D Mode, ICE integration and GRAPPA" tutorial repository! 

This tutorial presents steps that are necessary for converting a 2D GRE sequence, 
which runs but produces unstructured raw data, into a sequence that integrates 
with the scanner's UI and theimage reconstruction environment (ICE). We then integrate 
GRAPPA sampling pattern with embedded auto-calibration lines and mark the data to the 
ICE environment to accept them and produce images.

***s01\_GRE***

***s01*** describes a 2D slice-selective gradient echo (GRE) sequence with a
reasonably optimized timing. It is used as a starting point for further developments. 

***s02\_GRE_preregister***

***s02*** generates an MR pulse sequence equivalent to s01, but does it about 4 to 6 
times faster. This tutorial is designed to demonstrate the technique of pre-registering 
Pulseq objects to accelerate generation of longer MR pulse sequences. 

***s03\_GRE_2D_Mode***

***s03*** describes the steps necessary to switch the sequence from the default 
3D mode to a 2D multi-slice-selective sequence. For that sequence an example data 
set is available in the _data_ directory. The data here are fully unstructured, 
and are totally similar to s01 and s02, which is why no data of s01 were saved.

***s04\_2D_GRE_with_Labels***

***s04*** describes a 2D Cartesian GRE sequence, which uses the LABEL extension 
and produces data that can be seamlessly reconstructed to DICOM images on the scanner. 
The data (available in the data directory) are also compatible with ISMRM-RD 
converters and Gadgetron recon.

***s05\_2D_GRE_with_GRAPPA***

***s05*** is a 2D GRE sequence with the sampling pattern inplementing 2-fold parallel 
imaging accceleration in the way, that allows it to integrate with the iPAT/GRAPPA 
reconstruction in ICE.

***s06\_2D_GRE_with_GRAPPA_and_multiple_runs***

***s06*** shows steps necessary for making the sequence compatible with the **runs** 
feature of the interpreter in a way to remain compatible with the reconstruction in ICE.

## Quick instructions

Source code of the demo sequences and reconstruction scripts is the core of this repository. 
Please download the files to your computer and make them available to Matlab (e.g. by saving 
them in a subdirectory inside your Pulseq-Matlab installation and adding them to the Matlab's 
path). There are two sub-directories:

* seq : contains example pulse sequences specifically prepared for this demo
* data : contains the data acquired with some of the above sequences

## Quick links

Pulseq Matlab repository: 
https://github.com/pulseq/pulseq

## How to follow 

We strongly recommend using a text compare tool like *meld* (see this [Wikipedia page](https://en.wikipedia.org/wiki/Meld_(software)) and compare sequences from subsequent steps to visualise the respective steps.

## Further links

Check out the main *Pulseq* repository at https://github.com/pulseq/pulseq and familarising yourself with the code, example sequences and reconstruction scripts (see 
[pulseq/matlab/demoSeq](https://github.com/pulseq/pulseq/tree/master/matlab/demoSeq) and [pulseq/matlab/demoRecon](https://github.com/pulseq/pulseq/tree/master/matlab/demoRecon)). If you already use Pulseq, consider updating to the current version.


# Pulseq tutorial "From FID to GRE"

Welcome to the "From FID to PRESS" tutorial repository! This was initially developed for the Pulseq software demonstration and hands-on session at the **Italian Chapter of ISMRM 2022** in Pisa.

The tutorial starts with the very basic first steps from an FID and non-slective spin-echo sequence and moves on towards spectroscopic PRESS sequence and a set of imaging examples. The imaging part of the tutorial is based around the GRE sequence, where various degrees of sophistication are introduces step-by-step.

## Quick links

Pulseq Matlab repository: 
https://github.com/pulseq/pulseq

Dropbox link to the measured data and sequences (viewing only): 
https://www.dropbox.com/sh/l04pm0547yygswk/AAAP4mGeT5Ri0rk8uroGmyita?dl=0
 
## Quick instructions

Check out the main *Pulseq* repository at https://github.com/pulseq/pulseq and familarizing yourself with the code, example sequences and reconstructon scripts (see 
[pulseq/matlab/demoSeq](https://github.com/pulseq/pulseq/tree/master/matlab/demoSeq) and [pulseq/matlab/demoRecon](https://github.com/pulseq/pulseq/tree/master/matlab/demoRecon)). If you already use Pulseq, consider updating to the current version.

[Handout materials](Handout.pdf) TODO.

The source code of the demo sequences and reconstruction scripts is the core of this repository. Please download the files to your computer and make them available to Matlab (e.g. by saving them in a subdirectory inside your Pulseq-Matlab installation and adding them to the Matlab's path). There are two sub-directories:

* seq : contains example pulse sequences specifically prepared for this demo
* recon : contains the reconstruction scripts tested with the above sequences
* data : contains raw MR data in the Siemens TWIX format and the corresponding pulse sequences in the Pulseq format

## How to follow 

We strongly recommend using a text compate tool like *meld* (see this [Wikipedia page](https://en.wikipedia.org/wiki/Meld_(software)) and conpare sequences from subsequent steps to visualithe the respective steps.



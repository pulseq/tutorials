# Pulseq tutorial "Radial and Non-Cartesian"

Welcome to the "Radial and Non-Cartesian" tutorial repository! This was initially developed as a part of a live demo at the virtual ISMRM meeting in 2020 and was extended somewhat further thereafter. 

This tutorial presents two simple radial sequences that are derived from the corresponding Cartesian counterparts. Additionally a fast radial GRE sequence is presented and a very basic 2D spiral is itroduced. The slide deck entitled [12_Radial_and_nonCartesian.pdf](./doc/12_Radial_and_nonCartesian.pdf) shows sequence diagrams of all steps and visualises the changes at each step.

***s01\_CartesianSE***

***s01*** describes a 2D slice-selective SE sequence with 7 blocks in
each TR. Note that the polarity of the gy gradient in deltaTR (that is
the polarity of both phase encoding and rephrasing gradients) remains
unchanged due to the 180° refocusing RF pulse.

***s02\_RadialSE***

***s02*** describes a 2D slice-selective radial sequence, built based on
***s01*** using the mr.rotate function. The mr.rotate function rotates
the readout gradient and its pre-phaser around the *z*-axis to a certain
angle, which projects the gradient into the *x*- and *y-*axis and thus
produces two components.

***s03\_CartesianGradientEcho***

***s03*** describes a 2D slice-selective Cartesian GRE sequence, similar
as ***s01*** from the Tutorial 11\_from\_GRE\_to\_EPI.

***s04\_RadialGradientEcho***

***s04*** describes a 2D slice-selective radial GRE sequence built based
on ***s03***.

***s05\_FastRadialGradientEcho***

***s05*** describes a 2D slice-selective radial GRE sequence with the
shortest timing. The mr.addGradients function combines slice-selective
gradient and slice-refocusing gradient into a single "ExtendedTrapezoid"
gradient. The mr.align function is used to position the gxPre gradient
right after the RF pulse and before the end of the slice-refocusing
gradient. The flat time of the readout gradient is extended for
spoiling. The whole sequence has a total of 2 blocks per TR.

***s06\_Spiral***

***s06*** describes a 2D slice-selective spiral sequence with fat
saturation. A Gaussian RF pulse followed by a spoiling gradient in the
*z*-axis is used for fat saturation. A raw Archimedean spiral is
generated and then resampled to stay directly under the slew rate and
maximum gradient limits. The mr.traj2grad function calculates gradient
strength and slew rate based on k-space trajectory. The
mr.makeArbitraryGrad function generates an arbitrary gradient (e.g.
spiral).

## Quick instructions

Source code of the demo sequences and reconstruction scripts is the core of this repository. Please download the files to your computer and make them available to Matlab (e.g. by saving them in a subdirectory inside your Pulseq-Matlab installation and adding them to the Matlab's path). There are two sub-directories:

* seq : contains example pulse sequences specifically prepared for this demo
* recon : contains the reconstruction scripts tested with the above sequences

## Quick links

Pulseq Matlab repository: 
https://github.com/pulseq/pulseq

## How to follow 

We strongly recommend using a text compate tool like *meld* (see this [Wikipedia page](https://en.wikipedia.org/wiki/Meld_(software)) and compare sequences from subsequent steps to visualithe the respective steps.

## Further links

Check out the main *Pulseq* repository at https://github.com/pulseq/pulseq and familarizing yourself with the code, example sequences and reconstructon scripts (see 
[pulseq/matlab/demoSeq](https://github.com/pulseq/pulseq/tree/master/matlab/demoSeq) and [pulseq/matlab/demoRecon](https://github.com/pulseq/pulseq/tree/master/matlab/demoRecon)). If you already use Pulseq, consider updating to the current version.


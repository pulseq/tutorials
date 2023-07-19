![Pulseq live demo at MRI Together banner](./doc/mri_together_esmrmb_banner.png)
# Pulseq tutorial "From GRE to EPI"

Welcome to the "From GRE to EPI" tutorial repository! This was initially developed for the live Pulseq software demonstration and hands-on session at **MRI Together 2021** on-line conference. 

This topic introduces the way to establishm an Echo-Planar Imaging (EPI)
sequence based on a GRE sequence via multi-echo and segmented GRE sequences. 

[Handout materials](./doc/Handout.pdf) that accompanied the original session are available in the PDF format containing the slides from the presentation with some additional material and explanations. In particular, the second part of the document (page 18 and on) describe the code examples. Please keep in mind that these examples were specifically developed for the demo, so you might find some useful information there even if you are familiar with the previous versions of Pulseq. 

Additionally the slide deck entitled [11_from_GRE_to_EPI.pdf](./doc/11_from_GRE_to_EPI.pdf) shows the sequence diagrams of all steps and visualises the changes at each step.

***s01\_GradientEcho***

***s01*** is a single-echo GRE sequence with Ny\*numTE phase encodes.
Its echo time (TE) cycles by the defined TEs (\[4 9 15\] \* 1e-3). numTE
is the number of the defined TEs. Note: the mr.align function is very
useful to set delays of events within a block to achieve desirable
alignment.

***s02\_MultiGradientEcho***

***s02*** is a monopolar multi-echo GRE sequence. The timing is
calculated with the aid of helperT. gxFlyBack gradient is used to
rephase the readout gradient.

***s03\_BipolarMultiGradientEcho***

***s03*** is a bipolar multi-echo GRE sequence. It eliminates gxFlyBack
gradients by reversing the polarity of the even readouts.

***s04a\_SegmentedGradientEcho***

***s04a*** is a readout-segmented single-echo GRE sequence. Instead of
sampling the k-space for one readout during each TR, it samples the
k-space with multiple readouts during each TR (i.e. segmented readouts).
Bipolar readout gradients are used to avoid time loss due to additional
fly-back gradients. The area of phase-encoding (PE) gradients
(pre-dephasing and blip) is determined based on the number of segmented
readouts (nSeg).

***s04b\_SegmentedGradientEcho***

In s04b, the splitGradientAt function is used to split the gyBlip
gradient of ***s04a*** into two parts, each in a separate block, which
reduces the time interval between two adjacent segmented readouts. In
the first segmented readout block, the first part of gyBlip is added
after the readout gradient. In inner segmented readout blocks, the
second and first parts of gyBlip are combined and added before and after
the readout gradient. In the last segmented readout block, the second
part of gyBlip is added before the readout gradient.

In case the gyBlip.riseTime is shorter than the fallTime of the former
block, gyBlip.delay is increased, and thus the peak of gyBlip is moved
to the edge of the former block. In case gyBlip.riseTime is longer than
the fallTime of the former block, a delay is added to the later readout
gradient, such that the last point of gyBlip hits the last point of the
ramp-up of the later readout gradient.

***s04c\_SegmentedGradientEcho***

***s04c*** increases the nSeg of ***s04b*** from 5 to Ny, such that
***s04c*** traverses the whole k-space in one segmented scan (i.e. an
improvised EPI scan).

***s05\_EchoPlanarImaging***

***s05*** is a 2D EPI sequence with the shortest timing. It contains
Ny+3 blocks.

The mr.makeDigitalOutputPulse function defines the output trigger to
play out with every slice excitation.

The shortest timing is achieved by using maximum slew rate and ramp
sampling. gyBlip is split into two parts and distributed to two adjacent
blocks. The dead times of the readout trapezoid gradient are at the
beginning and the end (align with two parts of gyBlip), each equal to
half of gyBlip.Duration. The readout gradient is first calculated based
on the maximum slew rate (with consideration of the dead times). Then,
its amplitude is scaled down to fix the area to be Nx\*deltak.

The ADC dwell time on the readout flat-top is calculated based on
adcDwellNyquist=deltak/gx.amplitude. The dwell time is rounded down to
system.adcRasterTime. Note that the number of ADC samples on Siemens
should be divisible by 4. In addition, the ADC should be aligned with
respect to the readout gradient: *both Pulseq and Siemens define the ADC
samples to happen in the centre of the dwell period*.

***s06\_EPI\_SingleTraj***

***s06*** is a 2D EPI sequence constructed from a single arbitrary
trajectory. It is based on ***s05*** and contains 4 blocks.

A dummy sequence object (seq\_d) creates and exports the single EPI
trajectory in xy-plane with Ny phase encoding steps. A single ADC object
is created to sample from the start to the end of the single trajectory
(excluding the ADC dead times:
adcDur=seq\_d.duration-2\*sys.adcDeadTime). A real sequence object
(seq\_r) is created to combine the slice-selective excitation, Gz
rephasing and Gx and Gy dephasing, the single EPI trajectory and the
single ADC, and Gz spoiler together, for a total of 4 blocks.

## Quick links

Pulseq Matlab repository: 
https://github.com/pulseq/pulseq

Dropbox link to the measured data and sequences (viewing only): 
https://www.dropbox.com/sh/i7f1gpwyigdugps/AACd2jQJg_WjoTY2nqh7O8IHa?dl=0

## Quick instructions

Check out the main *Pulseq* repository at https://github.com/pulseq/pulseq and familarizing yourself with the code, example sequences and reconstructon scripts (see 
[pulseq/matlab/demoSeq](https://github.com/pulseq/pulseq/tree/master/matlab/demoSeq) and [pulseq/matlab/demoRecon](https://github.com/pulseq/pulseq/tree/master/matlab/demoRecon)). If you already use Pulseq, consider updating to the current version.

[Handout materials](Handout.pdf) that accompanied the original session are available in the PDF format containing the slides from the presentation with some additional material and explanations. In particular, the second part of the document (page 18 and on) describe the code examples. Please keep in mind that these examples were specifically developed for the demo, so you might find some useful information there even if you are familiar with the previous versions of Pulseq. 

Source code of the demo sequences and reconstruction scripts is the core of this repository. Please download the files to your computer and make them available to Matlab (e.g. by saving them in a subdirectory inside your Pulseq-Matlab installation and adding them to the Matlab's path). Yhere are two sub-directories:

* seq : contains example pulse sequences specifically prepared for this demo
* recon : contains the reconstruction scripts tested with the above sequences

The raw MR data, sequences in Pulseq format and examples of the reconstructed images can be accessed via the anonymous Dropbox link https://www.dropbox.com/sh/i7f1gpwyigdugps/AACd2jQJg_WjoTY2nqh7O8IHa?dl=0 

The Dropbox repository entails the subdirectory "dataPrerecorded", which contains data and sequence examples generated prior to the session. 

## Quick links

Pulseq Matlab repository: 
https://github.com/pulseq/pulseq

## Quick instructions

The source code of the demo sequences and reconstruction scripts is the core of this repository. Please download the files to your computer and make them available to Matlab (e.g. by saving them in a subdirectory inside your Pulseq-Matlab installation and adding them to the Matlab's path). There are two sub-directories:

* seq : contains example pulse sequences specifically prepared for this demo
* recon : contains the reconstruction scripts tested with the above sequences
* data : contains raw MR data in the Siemens TWIX format and the corresponding pulse sequences in the Pulseq format

## How to follow 

We strongly recommend using a text compate tool like *meld* (see this [Wikipedia page](https://en.wikipedia.org/wiki/Meld_(software)) and compare sequences from subsequent steps to visualithe the respective steps.

## Further links

Check out the main *Pulseq* repository at https://github.com/pulseq/pulseq and familarizing yourself with the code, example sequences and reconstructon scripts (see 
[pulseq/matlab/demoSeq](https://github.com/pulseq/pulseq/tree/master/matlab/demoSeq) and [pulseq/matlab/demoRecon](https://github.com/pulseq/pulseq/tree/master/matlab/demoRecon)). If you already use Pulseq, consider updating to the current version.



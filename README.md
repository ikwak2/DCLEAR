# Distance based Cell LinEAge Reconstruction (DCLEAR)

Il-Youp Kwak (<ikwak2@cau.ac.kr>) and Wuming Gong (<gongx030@umn.edu>)

R/DCLEAR is an R package for Distance based Cell LinEAge Reconstruction(DCLEAR). These codes are created during the participation of [Cell Lineage Reconstruction DREAM challenge](https://www.synapse.org/#!Synapse:syn20692755/wiki/595096).

## DCLEAR Overview 

Figure 1. Overview of DCLEAR modeling architecture. Our model is divided into two parts, 1) estimating distance between cells and 2) constructing tree using distance matrix.
<img src="https://ikwak2.github.io/tmpimages/modeling_overview.png" alt="drawing" width="850"/>

### Estimating distance between cells

Naive approach would be the hamming distance that simply calculate the edit distance.

<img src="https://ikwak2.github.io/tmpimages/hamming.PNG" alt="drawing" width="390"/>

However, the previous approach assume every base difference have same weights. For example, two sequences, '00AB0' and '0-CB0', are different at second and third positions. The second position, we have '0' and '-', and the third position, we have 'A' and 'C'.

For '0' and '-', '-' is point missing and it is possibly '0'. Thus it should have lower weight.
For 'A' and 'C', During the cell propagation, '0' differentiated to 'A' and '0' differentiated to 'C'. Thus it should have larger weight.
We can assign weights as below equation.

<img src="https://ikwak2.github.io/tmpimages/whamming.PNG" alt="drawing" width="500"/>

And we can approximate unknown weights using training data.

### Constructing tree from the distance matrix
With the previously proposed distance matric, we can construct distance matrix among cells. We can apply tree construction algorithms such as Neighbor-Joining(NJ), FastME.


## Usage

   - How to use weighted hamming : [Link](vignette/WHhowto.html)
   - How to use kmer_replacement : [Link](vignette/kmer_replacement_howto.html)
   - Preparation for subchallenge 2 submission : [link](vignette/PrepC2.html)
   - Preparation for subchallenge 3 submission : [link](vignette/PrepC3.html)
   - Another example using subchallenge 2 data : [link](vignette/Example_subchallenge2.html)


## installation

With 'devtools':
```S
devtools::install_github("ikwak2/DCLEAR")
```

## License

The R/DCLEAR package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>


## Presentation

Our talk on the special DREAM session in RECOMB 2020 meeting (https://www.recomb2020.org/) can be found [here](https://www.dropbox.com/s/a93q2lnqni6xf4q/RECOMB_2020_talk_final.pdf?dl=0).  

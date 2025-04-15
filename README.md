**Intro**
This repository contains the code for designing Binders using state-of-the-art RFdiffusion, ProteinMPNN and AF2 'initial guess'. 
Follow the Protocol in 'Binder Design_v1.docx'.

I included the repositories to install all the necessary Python libraries:

https://github.com/RosettaCommons/RFdiffusion?tab=readme-ov-file#the-inpaint_seq-flag
https://github.com/nrbennet/dl_binder_design?tab=readme-ov-file#setup4.

RFDiffusion can be easily run using Docker, the other libraries should be installed using pip or conda. 
ProteinMPNN and AF2 'initial guess' should be installed in Linux or in a WSL environment.

The .LOG file contains an **example output** of a binder designed to target **PD-L1** (PDB: 5O45). This example reproduces the results from  **J. L. Watson et al. Nature (2023)** (doi: 10.1038/s41586-023-06415-8).

Author: Goren Saenz de Pipaon

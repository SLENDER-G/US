# Ultrafast Ultrasound Beamforming as an Inverse Problem using a Low rank and Sparse Model
Code for Reproducing the Results Presented in This Paper.
# Abstract
Inverse problem formulations of ultrasound beamforming are a crucial branch of ultrasound computational imaging and have garnered increasing attention in recent years. This problem includes leveraging computational methodologies and mathematical modeling to unravel complex data embedded within the ultrasound beamforming process. Current sparse-based methods usually pose beamforming as a minimization problem of a fidelity term resulting from the measurement model plus a sparse regularization term that enforces the image contrast and resolution, but it performs poorly in the preservation of speckle texture. To address this issue, this paper proposes to use a low-rank and sparse approach to model the ultrasound plane wave imaging process as an inverse problem. The low-rank property enhances the exploration of correlations between pixels, leading to improved speckle preservation of the reconstructed image. Two optimization algorithms, based on alternating direction method of multipliers (ADMM), are proposed to efficiently solve the resulting optimization problem with either L1 or L0 constraints. The performance of the proposed approaches are evaluated using the data published in Plane Wave Imaging Challenge in Medical Ultrasound (PICMUS) in 2016. Results show that the proposed low-rank and sparse model is better adapt to the ultrasound plane wave beamforming process and has a good performance in speckle preservation. The proposed L1norm solution gives the best image quality in terms of resolution index while the L0-norm solution gives the best image contrast.
# Requirements
MATLAB (codes are tested on MATLAB R2022a)

PICMUS dataset

BFGS solver
# Instruction
1.Download the Required Datasets and Solvers:
Download the PICMUS dataset and the BFGS solver as mentioned in the setup requirements.
Create the Weighting Matrix (Phi):

2.After downloading the required files, you need to generate the weighting matrix Phi. To do this:
Run the following scripts in order:
weighting_matrix.m
summation_1.m
summation_2.m
These scripts will help you create the weighting matrix Phi.

3.Place Phi in the Correct Directory:
After Phi is created, save it in the same folder where the algorithm source codes are located.

We thank the organizers of the ultrasound toolbox and PICMUS challenge for providing publicly available data and codes. The authors sincerely thank Sobhan Goudarzi for making their MATLAB codes available online.
# License
[License](https://github.com/SLENDER-G/US/edit/main/LICENSE.txt) for non-commercial use of the software.Please cite the following [paper](https://doi.org/10.1109/TUFFC.2022.3198874) when using the codes. This code

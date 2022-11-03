# SFWCMF
Code for source localization with the covariance matrix fitting methods and Sliding Frank-Wolfe

Code and data to reproduce the results of XXXX

## Figures

The proposed methods (CMF and COMET with Sliding Frank-Wolfe) are compared with MUSIC, CLEAN-SC, OMP, and differential evolution (see references in the paper)

* `FIG_exp.m`: experimental results (noiseless and noisy versions)
* `FIG_SNR.m`: performance in function of the signal to noise ratio
* `FIG_K.m`: performance in function of the frequency
* `FIG_SNAPS.m`: performance in function of the number of snapshots
* `FIG_RESOL.m`: comparison of the resolution

## Scripts

Most of the scripts are for internal use or for comparison purposes.

The main scripts are:

* `sfw_cmf.m`: SFW algorithm for the CMF problem
* `sfw_comet1.m` and `sfw_comet2.m` : SFW algorithm for the two COMET variants

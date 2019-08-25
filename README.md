# Neighbor Denoising in Short Time Fourier Transform                

This repository contains MATLAB scripts and sample data for applying the denoising method presented in: 

`Mousavi, S. M., and C. A. Langston (2016). 
Adaptive noise estimation and suppression for improving microseismic event detection, 
Journal of Applied Geophysics, 132, 116-124, doi:http://dx.doi.org/10.1016/j.jappgeo.2016.06.008` 

------------------------

BibTeX:

    @article{mousavi2016adaptive,
    title={Adaptive noise estimation and suppression for improving microseismic event detection},
    author={Mousavi, S Mostafa and Langston, Charles A},
    journal={Journal of Applied Geophysics},
    volume={132},
    pages={116--124},
    year={2016},
    publisher={Elsevier}
    }

------------------------

## Paper
(https://www.researchgate.net/publication/305078128_Adaptive_noise_estimation_and_suppression_for_improving_microseismic_event_detection)

## Talk 
(https://earthquake.usgs.gov/contactus/menlo/seminars/1093)

------------------------

`demo.m` includes all info you need to know for running the code. 

you need `MATLAB statistics and signal processing toolboxes` to run this code.

------------------------

## A short description
In this approach for suppresing the noise from seismic data, first the noise level presented in the signal is estimated using 
the minima controlled recursive averaging technique. In this technique, past power values of noisy measurements during a period of signal absence are recursively averaged and the estimate is continued during signal presence. This is done by useing a time-varying frequency-dependent smoothing parameter that is adjusted by the probability of signal presence. The probabilities are obtained using Baye's theorem. 
After the noise estimation, denoising is done by thresholding the Short Time Fourier Transform coefficients based on a risk estimate (usindg Stein's unbiased risk estimate) from neighboring coefficients. 


![Denoising real seismic data. The left column shows presumably induced microseismic events due to wastewater injection in 
central Arkansas in 2010 recorded by a broadband seismometer at the surface. The right column shows the same trace and 
its STFT after denoising.](Fig.png)
Denoising real seismic data. The left column shows presumably induced microseismic events due to wastewater injection in 
central Arkansas in 2010 recorded by a broadband seismometer at the surface. The right column shows the same trace and 
its STFT after denoising.


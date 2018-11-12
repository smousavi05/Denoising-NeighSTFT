# Neighbor Denoising in Short Time Fourier   Transform

------------------------------------------------------

This repository contains MATLAB scripts and sample data for applying the denoising method presented in: 

Mousavi, S. M., and C. A. Langston (2016). Adaptive noise estimation and suppression for improving
microseismic event detection, Journal of Applied Geophysics, 132, 116-124, doi:http://dx.doi.org/10.1016/j.jappgeo.2016.06.008 


`demo.m` includes all info you need to know for running the code. 

you need `MATLAB statistics and signal processing toolboxes` to run this code.

## A short description 
In this approach for suppresing the noise from seismic data, first the noise level presented in the signal is estimated using 
the minima controlled recursive averaging technique. In this technique, past power values of noisy measurements during a period of signal absence are recursively averaged and the estimate is continued during signal presence. This is done by useing a time-varying frequency-dependent smoothing parameter that is adjusted by the probability of signal presence. The probabilities are obtained using Baye's theorem. 
After the noise estimation, denoising is done by thresholding the Short Time Fourier Transform coefficients based on a risk estimate (usindg Stein's unbiased risk estimate) from neighboring coefficients. 

## Paper
(https://www.researchgate.net/publication/305078128_Adaptive_noise_estimation_and_suppression_for_improving_microseismic_event_detection)

## Talk 
(https://earthquake.usgs.gov/contactus/menlo/seminars/1093)

## Abstract 
Microseismic data recorded by surface arrays are often strongly contaminated by unwanted noise. This background noise 
makes the detection of small magnitude events difficult. A noise level estimation and noise reduction algorithm is
presented for microseismic data analysis based upon minimally controlled recursive averaging and neighborhood shrinkage 
estimators. The method might not be compared with more sophisticated and computation- ally expensive denoising algorithm
in terms of preserving detailed features of seismic signal. However, it is fast and data-driven and can be applied in 
real-time processing of continuous data for event detection purposes. Results from application of this algorithm to 
synthetic and real seismic data show that it holds a great promise for improving microseismic event detection.

## A Short Description 
Seismic data recorded by surface arrays are often contaminated by unwanted noise. In many conventional seismic methods, 
the reliability of the seismic data and accuracy of parameter extraction, such as onset time, polarity, and amplitude, 
are directly affected by the background noise level. As a result, the accuracy of event location and other attributes 
derived from seismic traces are also influenced by the noise content. Therefore, there is a great need for developing 
suitable procedures that improve signal-to-noise ratios allowing for robust seismic processing. In this presentation, 
I introduce four different methods for automatic denoising of seismic data. These methods are based on the time-frequency 
thresholding approach. The efficiency and performance of the thresholding-based method for seismic data have been improved 
significantly. Proposed methods are automatic and data driven in the sense that all the filter parameters for denoising are 
dynamically adjusted to the characteristics of the signal and noise. These algorithms are applied to single channel data 
analysis and do not require large arrays of seismometers or coherency of arrivals across an array. Hence, they can be applied
to every type of seismic data and can be combined with other array based methods. Results show these methods can improve 
detection of small magnitude events and accuracy of arrival time picking.

![Denoising real seismic data. The left column shows presumably induced microseismic events due to wastewater injection in 
central Arkansas in 2010 recorded by a broadband seismometer at the surface. The right column shows the same trace and 
its STFT after denoising.](Fig.png)
Denoising real seismic data. The left column shows presumably induced microseismic events due to wastewater injection in 
central Arkansas in 2010 recorded by a broadband seismometer at the surface. The right column shows the same trace and 
its STFT after denoising.


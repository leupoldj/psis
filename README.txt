Here the code and data belonging to the paper from Leupold/Weigel/Bär ("On the choice of the phase difference increment in radiofrequency-spoiled gradient-echo magnetic resonance imaging of liquids with consideration of diffusion") is provided.

i) Simulations
--------------

The simulations are presented in different figures in the paper, a script named fig*_script.m is provided to reproduce any figure. Note all these scripts call the function epg_rfsp.m
Further explanation is given in the header of the individual scripts. (E.g., to reproduce Fig 6, the variable "substance" has to be manually set in the code.)

ii) Data
--------------
The calculation of the SNR data points as displayed in Fig.3 are provided by the following files:
- h2o_images.mat: A MATLAB-matrix containing the images used for ROI analysis of the H2O+CuSO4 phantom
- sil_images.mat: A MATLAB-matrix containing the images used for ROI analysis of the silicone oil phantom
- rois_h2o.m: ROI analysis of the H2O+CuSO4 phantom. Output: SNR values as displayed in Fig. 3
- rois_silicone.m: ROI analysis of the silicone oil phantom. Output: SNR values as displayed in Fig. 3

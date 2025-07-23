# Software-and-computing-
This project is part of my master thesis, in collaboration with ALICE. The Software used is O2Physics and the language used is ROOT c++. The codes inserted here are tasked one to do Particle Identification (readtree.cxx) while the other is generating data using Monte Carlo Method (readtreeMC.cxx) in order to be sure about the Particle Identification. 
Subsequently it is possible to verify the data analysis extrapolating the revelation efficiencies from one code and the raw spectra from the other, the rapport between these two gives us the completa spectra of the particles involved (in this case: pi, ka, pr).
The data used for the analysis comes from ALICE Run 3 with pp collisions at 13 TeV. 
This project is still undergoing development, as of now the objective is to perform Particle Identification with Bayesian statiscs and verifying that it could lead to improved results compared to the standard method.
Unfortunately I was not able to upload the Data tree, even if shrinked, due to the size of the file, but i uploaded some images that the code generated (only for pi+, and only in some cases, like the bin-to-bin analysis, otherwise the amount of images to upload would have been too much), the procedure is the iterated for all the other types of particles (pi-, pr+, pr-, ka+, ka-).
Most of the images are created in both codes, the data and the Monte Carlo one, in order to make sure that the results were compatible with one another.
The uploaded files are: energy loss; the time difference between the revelator, expected and actual of the particle; the study of the beta factor, the fit with a gaussian model divided into pt bins; the analysis bin-to-bin; the creation of the raw spectra; recovering the efficiencies for the revelator; the n sigma for the TPC (with mean and average set to 0 and 1); the n sigma for the TOF (same analysis).




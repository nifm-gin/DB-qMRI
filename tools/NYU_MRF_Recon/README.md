This toolbox reconstructs quantitative maps of arbitrary MR Fingerprinting data with arbitrary k-space trajectories. The algorithm is described in 

[1] J. Asslaender, M.A. Cloos, F. Knoll, D.K. Sodickson, J.Hennig and R. Lattanzi, Low Rank Alternating Direction Method of Multipliers Reconstruction for MR Fingerprinting  Magn. Reson. Med., epub ahead of print, doi:10.1002/mrm.26639, 2017.

The temporal evolution of the magnetization is approximated in a low rank space based on a singular value decomposition (SVD) of the dictionary. The function MRF_dictionary.m provides an example on how to calculate a dictionary and compress it. This function calculates a dictionary based on the pSSFP pattern described in 

[2] J. Asslaender, S. J. Glaser, and J. Hennig, Pseudo Steady-State Free Precession for MR-Fingerprinting, Magn. Reson. Med., p. epub ahead of print, doi:10.1002/mrm.26202, 2016.

The compression matrix resulting from the SVD is then passed on to the constructor of the LR_nuFFT_operator class. This class is based on Jeffrey Fessler's nuFFT implementation with min-max interpolation described in

[3] J. A. Fessler and B. P. Sutton, Nonuniform fast fourier transforms using min-max interpolation, IEEE Trans. Signal Process., vol. 51, no. 2, pp. 560-574, Feb. 2003.

The corresponding code can be found in the folder fessler_nufft, which contains a slightly modified version of Fessler's nuFFT. The original and full version of his famous toolbox can be found on

http://web.eecs.umich.edu/~fessler/code/

and I would like to acknowledge Jeffrey Fessler for his permission to make his code available here. 

The LR_nuFFT_operator exploits the commutative properties of the spatial FFT and the temporal compression. The constructor pre-calculates a sparse matrix that combines the interpolation in the k-space of each time-frame and the compression matrix along time. The mtimes function then performs an FFT of each low rank frame and grids the k-space of each low rank frame onto the k-space trajectory of the entire experiment.

The reconstruction itself is then performed by the function admm_recon.m, which couples the low rank approximation of the signal evolution in each voxel to the dictionary by an ADMM algorithm. It further allows also for spatial regularization with an l21-norm penalty. 

An example script can be found in the 'example' folder, which demonstrates the interface of the toolbox on the basis of a Shepp-Logan phantom. 

If you use this toolbox for the preparation of a manuscript, please cite above mentioned paper [1,3]. For questions, suggestions, bug-reports or if you would like to contribute to this toolbox, please feel free to contact me (see email address below).



License Information

The reconstruction framework is provided free of charge for use in research applications. It must not be used for diagnostic applications. The author takes no responsibility of any kind for the accuracy or integrity of the created data sets. No data created using the framework should be used to make diagnostic decisions. The source code is provided under the GNU General Public License (https://www.gnu.org/copyleft/gpl.html).

The framework, including all of its software components, comes without any warranties. Operation of the software is solely on the user’s own risk. The author takes no responsibility for damage or misfunction of any kind caused by the software.



_______________________________________
(c) Jakob Asslaender, August 2016, February 2017

New York University School of Medicine, Center for Biomedical Imaging

University Medical Center Freiburg, Medical Physics

jakob.asslaender@nyumc.org
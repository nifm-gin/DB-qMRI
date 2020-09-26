# DB-qMRI
Dictionary-based - quantitative magnetic resonance imaging

## Description

The dictionary-based learning (DBL) quantitative MRI method is proposed to bypass inherent magnetic resonance fingerprinting (MRF) limitations in high dimension: reconstruction time and memory requirement and provide both estimates and their confidence levels.

Standard  parameter estimation from magnetic resonance fingerprinting (MRF) data is based on matching the MRF signals to their best counterparts in a grid of coupled simulated signals and parameters, referred to as a dictionary [Ma et al.](https://doi.org/10.1038/nature11971).

To reach a good accuracy, the matching requires an informative dictionary whose cost, in terms of design, storage and exploration, is rapidly prohibitive for even moderate numbers of parameters. In this package, we propose an implementation of two dictionary-based learning (DBL) approaches made of three steps: 1) a quasi-random sampling strategy to produce efficiently an informative  dictionary, 2) a regression model to learn from the dictionary a correspondence between fingerprints and parameters (using either a neural network [Cohen et al.](https://doi.org/10.1002/mrm.27198) or an inverse statistical regression model [Boux et al.](https://hal.univ-brest.fr/INRIA/hal-02314026v2)), and 3) the use of this mapping to provide parameter estimates (and their confidence indices for statistical method).

Details about these methods referred to as dictionary-based matching (DBM), dictionary-based deep learning (DB-DL) and dictionary-based statistical learning (DB-SL) can be find in [Boux et al.](https://hal.univ-brest.fr/INRIA/hal-02314026v2).

## Configuration

The code has been verified with Matlab R2018 and R2019.

```Statistics and Machine Learning Toolbox``` and ```Parallel Computing Toolbox``` toolboxes are required.


## Run
Figures from different experiments can be found in the `./figures` folder. To regenerate paper figures, the best way is to run the `Run.m` script:
```matlab
>> Run
```
Information about figures are described (see comments) in the `Run.m` file. 
Experiments can be launched individually by executing the scripts located in the folder `./Experiments`.

The quantification is achieved, running:
```matlab
>> [Estimation, Parameters] = AnalyzeMRImages(Sequences, Dico, Method)
```
where ```Sequences``` is a 3D or 4D matrix of observed MR signals (the third dimension is the time, others are spatial dimensions), ```Dico``` is a structure that represents the dictionary and ```Method``` is the strings ```'DBM'```, ```'DB-SL'``` or ```'DB-DL'``` to specify the method to use. The fields of ```Dico``` are ```Dico.MRSignals``` that is a 2D matrix of MR signals (second dimension is time) and ```Dico.Parameters.Par``` is a 2D matrix of parameters that match to the corresponding MR signals (second dimension is the parameter dimension). Then, note that the first dimensions of ```Dico.MRSignals``` and ```Dico.Parameters.Par``` must be equals.

```Estimation``` and ```Parameters``` are structures. ```Estimation.Y``` is the matrix of parameter estimates.


## External tools
The `./tools` folder contains external toolboxes located in the subfolder having the same name:

- Antoine Deleforge, the [GLLiM](https://team.inria.fr/perception/gllim_toolbox/) regression.

- Jakob Asslaender, the [NYU_MRF_Recon](https://bitbucket.org/asslaender/nyu_mrf_recon/src/master/) toolbox reconstructs quantitative maps of arbitrary MRF data with arbitrary k-space trajectories. The tool has been modify in order to take into account any sampling during the dictionary generation.


## Contact information
Fabien Boux, <fabien.boux@univ-grenoble-alpes.fr>

# DBL-qMRI
The dictionary-based learning (DBL) method is proposed to bypass inherent magnetic resonance fingerprinting (MRF) limitations in high dimension: reconstruction time and memory requirement and provide both estimates and their confidence levels.

## Configuration

The code has been verified with Matlab R2018 and R2019.

The ```Statistics and Machine Learning Toolbox``` and ```Parallel Computing Toolbox``` toolboxes are required.


## Run
All figures can be found in the `./figures` folder. To reproduce paper figures, the best way is to run the `Main.m` script:
```matlab
>> Main
```
The information about figures is given as comments in the file `Main.m` script.

The `./tools` folder contains the [GLLiM regression](https://team.inria.fr/perception/gllim_toolbox/) of Antoine Deleforge and our adaptation in the MRI context that can be run using:
```matlab
>> [Estimation, Parameters] = AnalyzeMRImages(Sequences,Dico,Method)
```
where ```Sequences``` is a 3D or 4D matrix of observed MR signals (the third dimension is the time, others are spatial dimensions), ```Dico``` is a structure that represents the dictionary and ```Method``` is the strings ```'DBM'``` or ```'DBL'``` that calls the method. The fields of ```Dico``` are ```Dico.MRSignals``` that is a 2D matrix of MR signals (second dimension is time) and ```Dico.Parameters.Par``` is a 2D matrix of parameters that match to the corresponding MR signals (second dimension is the parameter dimension). Then, note that first dimensions of ```Dico.MRSignals``` and ```Dico.Parameters.Par``` must be equals.

```Estimation``` and ```Parameters``` are structures. ```Estimation.DBM.Y``` or ```Estimation.DBL.Y``` are matrices of parameter estimates.

## Contact information
Fabien Boux, <fabien.boux@univ-grenoble-alpes.fr>

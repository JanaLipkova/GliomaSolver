# Glioma solver: Tumor growth with mass effect
This example simulates tumor growth in patient anatomy, accounting for tumor induced brain deformations and distribution of intracranial pressure (ICP).

### Compilation & Running
Set up enviroment, e.g. on the local computer called kraken and compile to get executable ```brain```. Put the executable in the folder ```BrainDeform``` and run the simulation using the ```run_kraken.sh```:
```sh
cd makefile
source setup_kraken.sh
make clean && make helmholtz=hypre -j 4
cp brain ../Examples/BrainDeform 
cd ../Examples/BrainDeform
./run_kraken.sh
```
For running on SLURM cluster, use the scripts with keyword lrz and submit as  ```sbatch run_lrz.sh ```. The  ```run_*.sh ``` contains execution command and simulaiton input parameters.

### Input parameters

| Parameter        | Description     |
| ------------- |:-----------------|
| `N`  | Number of OMP threads to be used|
| `M`  | Number of MPI ranks to be used|
| `program`   | Name of executable, here `brain` |
| `model`   | The name of simulation model to be used, use `deform` for Brain deformation model. |
| `verbose`   | Set to `0/1` to disable/enable simulation printouts |
| `profiler`   | Set to `0/1` to disable/enable code profiling (for performance tests) |
| `bDumpIC`   | Set to `0/1` to disable/enable saving of the initial condition |
| `dumpfreq`   | Frequency of data dumping (w.r.t simulation time, i.e every 50 days) |
| `vtk`   | Set to `0/1` to disable/enable dumping output |
| `CFL`   | CFL condition for advection/convection operator |
| `Tend`  | Final simulation time [day] |
| `rho`   | Proliferation rate [1/day] |
| `Dw`    | Diffusivity in white matter [cm^2/day] |
| `kCSF, kWM, kGM`   | Relaxation in CSF, WM, GM |
| `bMobility`   | Set to `0/1` to disable/enable tissue dependent mobility |
| `ICtype`   | Use `1` for brain anatomy, `0` for simplify model on sphere |
| `PatFileName`   | Path to file with patient anatomy |

### Input anatomy
In the example anatomy of the patient is stored in `/Anatomy/Patient00/`. The input data should be provided in the format `.dat`. See `ProcessingData/source/convert_folder_content_nii2dat.m' for script converting nifty to dat and vice-versa. For each patient provide segmentation of wm, gm, csf and phase-field function (pff), each in separate file with following name extension:
* *_WM.dat for wm
* *_GM.dat for gm
* *_CSF.dat for csf
* _PFF.dat for pff

 
### Output
If `vtk` is set to `1`, the solver will save output in the form `Data_000*.vtu`, where `*` is numbering of output data in chronological order, i.e. `Data_0000.vtu` is the inital conditon and `Data_0001.vtu` is the solution at time specified in 'dumpfreq'. The output files can be open in Paraview.

### References
Please cite:
* Lipkova et al.: *Personalized Radiotherapy Planning for Glioma Using Multimodal Bayesian Model Calibration*, preprint arXiv:1807.00499, (2018)
* Rossinelli D, et al.: *Mrag-i2d: Multi-resolution adapted grids for remeshed vortex methods on multicore architectures.* Journal of Computational Physics 288:1�~@~S18, (2015).








# Glioma solver
Softvare for simulating tumor growth in brain anatomy, including models for:
* tumor growth (Reaction-Diffusion model)
* tumor mass effect and ICP (Brain deformation model)
* others

### Instalations
* Required libaries: `tbb`, `vtk` and for deformation model also `hypre`
* Downloaded the UNIX compatible libraries here: http://tdo.sk/~janka/lib/ 
* Unpack libraries and store to a folder, i.e. /lib/ 
* Install `tbb` libraries by calling ```make clean && make```
* `vtk` library is already precompiled, you just need to uncompressed it, no need to install
* For the deformation solver install also `hypre` (follow instructions in the hypre-2.10.0b)

### Set-up enviroment
In following notes, files with a word *kraken* in the name refer to files for unix local computer, while files with a word *lrz* in the name refer to SLURM cluster files. Pick one, depending on your enviroment.
1) Set path to your libraries either in your ```.bashrc``` or similar file or create set-up script. Folder ```makefile/``` contains two examples of setup scripts: ```setup_kraken.sh``` and ```setup_lrz.sh```. If using the example script, change ```LIB_BASE``` variabel to point to your libraries.
2) Create local makefile for your enviroment.  Folder ```makefile``` contains two examples of local makefiles ```make.kraken``` and ```make.lrz.sh```. 
    * Copy one of the example make file to ```make.your_host_name```, where *your_host_name* is your hostname, i.e. ```make.jana```. 
    * Inside your ```make.jana``` change ```LIB_BASE``` variabel to point to your libraries
3) Main ```Makefile``` calls corresponding local  ```make.* ``` file depending on the hostname. Modify it so it calls your  ```make.jana ```.

### Compilation & Running
Set up enviroment, e.g. on the local computer called kraken. The compile:
```sh
source setup_kraken.sh
make clean && make -j 4
% For deformation model:
make clean && make helmholtz-hypre -j 4
```
Creates executable ```brain```. 

### Examples
Folder `makefile` contains example of solver:
* `TumorGrowth`: tumor growth model in patinet antomy using *reacion-diffusion model*
* `BrainDeform`: tumor growth with in patient anatomy with mass effect ucing *brain-deformaiton model*

### References
Please cite:
* Lipkova et al.: *Personalized Radiotherapy Planning for Glioma Using Multimodal Bayesian Model Calibration*, preprint arXiv:1807.00499, (2018)
* Rossinelli D, et al.: *Mrag-i2d: Multi-resolution adapted grids for remeshed vortex methods on multicore architectures.* Journal of Computational Physics 288:1–18, (2015).
    


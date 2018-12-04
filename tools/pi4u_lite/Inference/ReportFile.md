% Results of the Bayesian Calibration
% GliomaSolver 

# INPUT
The Bayesian calibration was performed for

* the input data provided in the folder 
    * *INPUT_DATA_FOLDER*
* Using *NUMBER_OF_SAMPLES* samples.

# RESULTS
The results include:

1. **The calibrated model parameters** reported in `Table 1` and `Table 2`.
2. **The posterior probability distribution (PDF)** 
    * With manifold shown in `Figure 1.`
    * And samples stored in `Results/Details/posterior.txt`
3. **The most probable tumor cell density**, given by the maximum a posteriori (MAP) estimate is:
    * Saved in `Results/MAP.nii`
    * With a preview shown in `Figure 2.`


### Calibrated parameters

* The calibrated parameters are reported in:
    * `Table 1`: for the tumor growth model 
    * `Table 2`: for image related parameters
* Reported is MAP, mean and standard deviation (std) of the marginal distribution. 
* The units are 
    * *Dw [cm2/day]; rho [1/day]; T [day];* 
    * *(icx, icy, icz)* is reported in dimension-less units.


*`Table 1:` The calibrated parameters of the tumor growth model:*

|   | *Dw*  | *rho*  | *T* |  *icx* | *icy* | *icz* |
|---|---|---|---|---|---|---|
| MAP  | _MAP_p1_  |  _MAP_p2_  | _MAP_p3_  | _MAP_p4_  | _MAP_p5_  | _MAP_p6_ |
| mean | _mean_p1_ |  _mean_p2_ | _mean_p3_ | _mean_p4_ | _mean_p5_ | _mean_p6_ |
| std  | _std_p1_  |  _std_p2_  | _std_p3_  | _std_p4_  | _std_p5_  | _std_p6_ |

*`Table 2:` The calibrated parameters of the image model:*

 |   | *sigma*  | *b*  | *ucT1*  |  *ucFLAIR*  | *sigma-alpha* |
|---|---|---|---|---|---|
| MAP  | _MAP_p7_  | _MAP_p8_  | _MAP_p9_  | _MAP_p10_  | _MAP_p11_  |
| mean | _mean_p7_ | _mean_p8_ | _mean_p9_ | _mean_p10_ | _mean_p11_ |
| MAP  | _std_p7_  | _std_p8_  | _std_p9_  | _std_p10_  | _std_p11_  |


### Posterior probability distribution 

![The results of the Bayesian calibration. **Above the diagonal:** Projection of the TMCMC samples of the posterior distribution in 2D parameter space. The colors indicate likelihood values of the samples. The number in each plot shows the Pearson correlation coefficient between the parameter pairs. **Diagonal:** Marginal distributions obtained with Gaussian kernel estimates. **Below the diagonal:** Projected densities in 2D parameter space constructed by 2D Gaussian kernel estimates.](Details/PosteriorPDF.jpg)

 
### Tumor Cell Density 

![Overview of the input modalities and the inferred MAP tumor cell density. Shown is a middle slice across the MAP tumor cell density.](Details/Overview.jpg)


# References
Please cite

* Lipkova et al.: *Personalized Radiotherapy Planning for Glioma Using Multimodal Bayesian Model Calibration*, preprint arXiv:1807.00499, (2018)
* Rossinelli D, et al.: *Mrag-i2d: Multi-resolution adapted grids for remeshed vortex methods on multicore architectures.* Journal of Computational Physics 288:1–18, (2015).
* P. E. Hadjidoukas et al.: *“Pi4u: A high performance computing framework for bayesian uncertainty quantification of complex models,* Journal of Computational Physics, 2015.
 
#### SEE ALSO
The GliomaSolver source code and all documentation may be downloaded from
<http://COMING-SOON/>.


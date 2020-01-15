## Software for realistic tumor modeling
_Tool for simulating tumor dynamics in patient-specific anatomy reconstructed from the medical scans._

#### Supported tumor models:
 * Tumor growth (Reaction-diffusion model)
 * Tumor mass effect and intracranial pressure (Brain deformation model)
 * Tumor hypoxia and necrosis
 * Bayesian model calibration for personalized radiotherapy planning 

#### Features:
 * **Fast execution** thanks to *highly-parallel architecture* and *adaptive grid refinement*
 * **User-friendly**; just provide patient anatomy and select the tumor model
 * **Developer-friendly**: the models are implemented as *modules*, which allows an *addition of new models* and a *combination* of the existing ones
 * **Transferable**: can be applied to other types of *infiltrative tumors* (e.g. multiple myeloma, liver lesions,...)
 * **Open source**, self-contained *C++ code*
 * Compatible with **Linux/Mac OS**

#### Software home-page
Please visit the [GliomaSolver](http://tdo.sk/~janka/GliomaWebsite//index.html) homepage for installation, tutorials, sample of glioma data, and much more :panda_face:.

#### Data & Resources
For additioanl data and softwares, please visit the solver's homepage, section [References](http://tdo.sk/~janka/GliomaWebsite/references.html). The most used links are also listed below.
 * Data used for the personalized radiotherhpy planning in [1] are avaialble [here](http://tdo.sk/~janka/GliomaSolverHome/GlioblastomaDATA.tgz)
 * A brain atlas and precomputed phase-field function, used for the deformation model, can be found [here](http://tdo.sk/~janka/GliomaSolverHome/TutorialTestData/AtlasPFF.tgz)
 * For automated image-registration tool see [here](https://github.com/JanaLipkova/Registration)
 * For automated skull-stripping and brain tissue segmentation tool see [here](https://github.com/JanaLipkova/s3)
 

#### References
In publications using GliomaSolver or data released with this solver pleace cite:

[1] Lipkova et al., *Personalized Radiotherapy Design for Glioblastoma Using Mathematical Tumor Modelling, Multimodal Scans and Bayesian Inference. IEEE Transactions on Medical Imaging (2019).* 









# README

Thank you for your interest in our work!

This repository hosts the **diskSED** model, as described in [Guolo & Mummery (2025)](https://arxiv.org/abs/2408.17296).

## Requirements

### Basic (Standard XSPEC Minimization)
- Python 3
- Astropy
- NumPy
- pyXSPEC (included as part of [HEASoft](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/))

### Bayesian Analysis (Recommended)
- All of the above
- Bayesian X-ray Analysis ([BXA](https://johannesbuchner.github.io/BXA/index.html))
- [UltraNest](https://johannesbuchner.github.io/UltraNest/index.html) (or alternative Bayesian inference software)
- [Corner](https://corner.readthedocs.io/en/latest/) (for corner plots)

## Installation

### Standard Installation
```bash
pip install git+https://github.com/muryelgp/diskSED.git
```

### Flexible Installation
```bash
git clone https://github.com/muryelgp/diskSED.git
cd diskSED
python3 -m pip install -e .
```


## Tutorials

The **Examples** folder contains several [Jupyter notebooks](https://jupyter.org) demonstrating step-by-step usage of the model:

- **`prep_data.ipynb`**: Demonstrates how to prepare grouped X-ray spectra and UV/optical/NIR SED files.
- **`fitting.ipynb`**: Shows how to load data, configure the model, and perform Bayesian fitting using pyXSPEC/BXA/UltraNest.
- **`fancy_plots.ipynb`**: Provides instructions for generating plots of data, observed and intrinsic model SEDs, corner plots, and parameter histograms.
- **`properties.ipynb`**: Explains how to estimate secondary parameters from the fit results (e.g., Black Hole Mass, outer radius in Rg, Bolometric Luminosity, and Eddington Ratio).
- *Planned Updates:* Advanced techniques like multi-epoch joint fitting, non-standard priors, or non-standard likelihood functions (e.g., including upper limits), aditionon of non-thermal components, etc. See papers in **Citations** for examples.

## Citations

If you use the **diskSED** model in any publication, please cite [Guolo & Mummery (2025)](https://arxiv.org/abs/2408.17296) along with its dependencies.

#### Examples of papers using diskSED
- *Add GSN069*
- *Add eRO-QPE2*

These papers provide examples of how to describe the model and properly cite dependencies. If these works have influenced or motivated your research, you are encouraged to cite them as well.

## Recomended Reading

#### Statistics
* [Statistical Aspects of X-ray Spectral Analysis](https://arxiv.org/abs/2309.05705) by J. Buchner & P. Boorman.
* [Dos and don'ts of reduced chi-squared](https://arxiv.org/pdf/1012.3754) by Andrae et al.

#### Relevant Software Tutorial/Lectures
* [BXA tutorial](https://peterboorman.com/tutorial_bxa.html) by P. Boorman.
* [Practical Inference for Researchers in the Physical Sciences](https://johannesbuchner.github.io/PracticalInferenceForResearchersInThePhysicalSciencesCourse/) by  F. Capel & J. Buchner.

#### Time-dependent Disk Theory
* [The Disk Accretion of a Tidally Disrupted Star onto a Massive Black Hole](https://ui.adsabs.harvard.edu/abs/1990ApJ...351...38C/abstract) by Cannizzo et al.
* [The spectral evolution of disc dominated tidal disruption events](https://arxiv.org/abs/1912.06577) by Mummery & Balbus
* [Fundamental scaling relationships revealed in the optical light curves of tidal disruption events](https://arxiv.org/abs/2308.08255) by Mummery et al.
* [Fitting transients with discs (FitTeD): a public light curve and spectral fitting package based on evolving relativistic discs](https://arxiv.org/abs/2408.15048) by Mummery et al.
  
## Comments & Contact

- If you're looking for the Relativistic version (**kerrSED**), it is unfortunately not publicly available yet.
- The model is implemented **only in Python** and can therefore only be used with **pyXSPEC**. A C/Fortran version does not exist, and therefore the model is not compatible with standard XSPEC. Recommendation: Switch to pyXSPEC for greater flexibility, especially when dealing with high-dimensional models/data-sets and to use Bayesian inference methods.
- If you're interested in implementing **diskSED** in standard XSPEC or have already done so, feel free to contact me so we can include it here.
- Please note that we are not responsible for any misuse of the model, nor do we endorse results derived from it a priori. Users are encouraged to apply caution and sound scientific judgment. The model's assumptions and limitations are described in the original paper, particularly on sections 1 and 5, aditional questions can be send by email.
- For general inquiries, contact me at **mguolop1@jhu.edu**. Please check the tutorials first to see if your question has already been addressed.

Thank you for using **diskSED**!


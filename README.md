# README

Thank you for your interest in our work!

This repository hosts the **diskSED** model, as described in [Guolo & Mummery (2025)](https://arxiv.org/abs/2408.17296).

## Requirements

### Basic (Standard Minimization)
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
- *Planned Updates:* Advanced techniques like multi-epoch fitting, non-standard priors, or non-standard likelihood functions (e.g., including upper limits). See papers in **Citations** for examples.

## Citations

If you use the **diskSED** model in any publication, please cite [Guolo & Mummery (2025)](https://arxiv.org/abs/2408.17296) along with its dependencies.

#### Example of papers using diskSED
- *Add GSN069*
- *Add eRO-QPE2*

These papers provide examples of how to describe the model and properly cite dependencies. If these works have influenced or motivated your research, you are encouraged to cite them as well.

## Recomended Reading

* [Statistical Aspects of X-ray Spectral Analysis](https://arxiv.org/abs/2309.05705) by J. Buchner & P. Boorman.
* [BXA tutorial](https://peterboorman.com/tutorial_bxa.html) by P. Boorman.
* [Practical Inference for Researchers in the Physical Sciences](https://johannesbuchner.github.io/PracticalInferenceForResearchersInThePhysicalSciencesCourse/) by  F. Capel & J. Buchner.
  
## Comments & Contact

- If you're looking for the Relativistic version (**kerrSED**), it is unfortunately not publicly available yet.
- The model is implemented **only in Python** and can therefore only be used with **pyXSPEC**. A C/Fortran version does not exist, and therefore the model is not compatible with standard XSPEC. Recommendation: Switch to pyXSPEC for greater flexibility, especially when dealing with high-dimensional models/data-sets and to use Bayesian inference methods.
- If you're interested in implementing **diskSED** in standard XSPEC or have already done so, feel free to contact me so we can include it here.
- We are **not responsible for any misuse** of the model. Its assumptions and limitations are described in the original paper.
- For general inquiries, contact me at **mguolop1@jhu.edu**. Please check the tutorials first to see if your question has already been addressed.

Thank you for using **diskSED**!


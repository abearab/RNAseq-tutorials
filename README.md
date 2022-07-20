# Setup your system
1. Install `conda`
2. Add relative channels
```bash
# add channels 
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
conda config --add channels python
```
3. Create fresh environment for your task and install packages locally. 

## Exploratory data analysis
__How to use Jupyter with multiple conda environments?__
Install `nb_conda_kernels` in the base environment or separate environment you used for jupyter stuff.
```bash
conda install -c conda-forge nb_conda_kernels
```
Using `nb_conda_kernels`, you can have one Jupyter installed in your system and launch different python or R kernels form any created conda environments even in a single notebook.
> Note: You only need ipykernel, numpy and pandas in each environment in addition to your own packages.
```bash
conda install -n <env-name> -c anaconda ipykernel numpy pandas
```
How to use R (and python) in Jupyter?
Instead, you can include R kernel into an envrinment with R packages. So, install irkernel.
How to tunnel from your local machine to the server?
I could almost always use this: ssh -NfL 1111:localhost:2222 <user>@<server>.ucsf.edu
Replace 2222 with Jupyter port number in the server and 1111 with your desired local port number.

# Tutorials 
## Part 1: Alignment
__[Tutorial notes](part-1-alignment/README.md)__

## Part 2: Differential expression analysis 
__[Tutorial notes](part-2-diff-exp/README.md)__

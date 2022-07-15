# Setup your system
1. Install `conda`
2. Add relative channels
3. Create fresh environment for your task

```bash
# add channels 
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
conda config --add channels python
```
## Exploratory data analysis
Jupyter
- https://anaconda.org/anaconda/jupyter
```bash
conda activate base 
conda install -c anaconda jupyter
# run jupyter
jupyter notebook
```
# Tutorials 
## Part 1: Alignment
__[Tutorial notes](part-1-alignment/README.md)__

## Part 2: Differential expression analysis 
__[Tutorial notes](part-2-diff-exp/README.md)__

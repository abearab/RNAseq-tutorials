# DESeq2 analysis task
- http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
___
0. Use relative environment
```bash
conda create -n deseq -c bioconda bioconductor-deseq2
```

```bash
conda activate deseq
```
1. Run Jupyter
```bash
jupyter lab
```
OR
```bash
nohup jupyter lab &> jupyter.log &
```
You may tunnel to the server using this trick:
```bash
ssh -NfL <port#>:localhost:<port#> user@server.etc
```
2. Open `deseq.ipynb` notebook and follow steps.

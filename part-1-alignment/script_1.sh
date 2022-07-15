```
setup your system:
1. install conda 
2. add relative channels
3. create environment for your task
```
# add channels 
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
conda config --add channels python

# make environment 
conda create -n alignment salmon star kallisto 
# activate it 
conda activate alignment


```
Alignment task: 
1. download reference genome
2. make genome index 
3. align reads (fastq files) to them
```

# Download ref genome to do the alignment
cd ~/genome/
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Salmon 
# https://combine-lab.github.io/salmon/getting_started/

# build index
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i Homo_sapiens.GRCh38.salmon_index

mkdir quants
# alignment
for f in fastq/*; do
	samp=`basename ${f}`
	samp=${samp/.fastq.gz/};
	echo "Processing sample ${samp}"
	salmon quant -i genome/Homo_sapiens.GRCh38.salmon_index/ \
	-l A -r $f -p 8 --validateMappings -o quants/$samp
done

# kallisto
# https://pachterlab.github.io/kallisto/starting

# build index
kallisto index -i Homo_sapiens.GRCh38.kallisto.transcripts.idx \
Homo_sapiens.GRCh38.cdna.all.fa.gz

mkdir kallisto_quant
# alignment
for f in fastq/*; do
	samp=`basename ${f}`
	samp=${samp/.fastq.gz/};
	echo "Processing sample ${samp}"
	echo kallisto quant -i genome/Homo_sapiens.GRCh38.kallisto.transcripts.idx \
	-o kallisto_quant/$samp -b 100 --single -l 180 -s 20 $f
done

```
Install jupyter and start doing exploratory data analysis
# https://anaconda.org/anaconda/jupyter
```
conda activate base 
conda install -c anaconda jupyter
# run jupyter
jupyter notebook
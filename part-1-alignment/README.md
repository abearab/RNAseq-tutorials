# Alignment task
__Salmon__
- https://combine-lab.github.io/salmon/getting_started/

__kallisto__
- https://pachterlab.github.io/kallisto/starting
___
0. Use relative environment
```bash
# make environment 
conda create -n alignment salmon star kallisto 
```
```bash
conda activate alignment
```
1. Provide reference genome
You may download it from ensembl FTP:
```bash
# Download ref genome to do the alignment
cd ~/genome/
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```
2. Make genome index
__Salmon__
```bash
# build index – salmon
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz \
-i Homo_sapiens.GRCh38.salmon_index
```
__kallisto__
```
kallisto index -i Homo_sapiens.GRCh38.kallisto.transcripts.idx \
Homo_sapiens.GRCh38.cdna.all.fa.gz
```
3. Align reads (fastq files) to reference genome
__Salmon__
```bash
mkdir quants
# alignment
for f in fastq/*; do
	samp=`basename ${f}`
	samp=${samp/.fastq.gz/};
	echo "Processing sample ${samp}"
	salmon quant -i genome/Homo_sapiens.GRCh38.salmon_index/ \
	-l A -r $f -p 8 --validateMappings -o quants/$samp
done
```
__kallisto__
```bash
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

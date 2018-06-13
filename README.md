# MNGS Training
See full description at the [wiki](https://github.com/AlfredUg/MNGS/wiki).

### **Introduction**
This is a follow up on the series of bioinformatics trainings that happen at Uganda Virus Research Institute facilitated by H3ABioNet. This particular one focuses on the main procedure of metagenomics NGS data  analysis. Custom scripts are used but the users should take note that, the required software is available and accessible on respective machines. It may also be necessary to change paths to main directories to fit ones file system. This can always be changed at the top of each script.

### **Getting started**
Create a directory named "MNGS2" in your home directory. This is where the analysis will be carried out. As we go on, we shall be creating sub-directories accordingly to ensure that intermediate outputs can be refered to easily whenever needed for downstream analyses. 
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir MNGS2
```
Make a sub-directory `data` in MNGS2 and create soft(symbolic) links to the real data. 
This is helpful since in case one accidentaly deletes them, we do not lose the original data. 
<!--In addition to avoiding loss of data, the links occupy a small memory, as save some space on disk. -->
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir MNGS2/data
ln -s "/path/to/dataset/*" "$HOME/MNGS2/data"
```
To confirm that the links have been successfully created, list the contents of the folder using the `ls` command. If the data is compressed (i.e .gz files), uncompress them using `gunzip *`

### **Quality control check**
Here we use `fastqc` [@fastqc], it is necessary that we store quality control files for easy reference. In MNGS2, create a sub-directory qcresults, this is where the fastqc results will be stored.
After that, do the quality checks;
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/MNGS2/qcresults"
fastqc "$HOME/MNGS2/data/*" -o "$HOME/MNGS2/qcresults"
```
Once this is done, navigate to `qcresults` and download the ".html" files to local machine and open them in any browser. This report can be used to assess quality distribution, length distribution, GC-content, nucleotide distribution. This informs downstream analysis.

### **Data quality trimming**
We use `trim_galore` for quality and adapter trimming. Depending on the qc results, it would be necessary to change some parameters used here accordingly. This script clips off the first 16 bases of the reads from the 5' end. In addition, it removes bases with phred quality less than 25 on the 3' end of the reads. We need to store quality trimmed reads, as such, in MNGS2 directory, make a sub directory `trimmed`, Trim data using the script `trim.sh`.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/MNGS2/trimmed"
bash "$HOME/MNGS2/scripts/trim.sh"
```
If all goes well, trimmed reads will be available in trimmed.
Below is the details of trim.sh script, it iteratively trims all fastq files in the `fastqdir` directory. You may consider looking at the trimmed reads using `fastqc` to check the improvement made by `trim_galore`.

#### **Remove host sequences**
There is need to get rid of host DNA sequences that could have contaminated the data earlier in the sampling and extraction stage. This is done by mapping the reads to host reference genome and picking the unmapped sequences. At this point, we create a directory to store host-free sequence data. `Samtools` [@samtools] and `bedtools` [@bedtools] are the key tools used for this purpose.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/MNGS2/hostfree"
bash "$HOME/MNGS2/scripts/removehost.sh"
```

### **Read mapping**
Mapping to reference sequences is done using `bowtie2` [@bowtie2]. For this training, the references were indexed apriori and are available in the directory (`"$HOME/MNGS2/refs"`). First, make a directory `alignment` which will contain the alignment results. The script `align.sh` aligns reads to each of the reference sequences, converts .sam to .bam files and sorts them, index the sorted files, generates summary of mapped and unmapped reads and ultimately extracts mapped reads.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/MNGS2/alignment"
bash "$HOME/MNGS2/scripts/align.sh"
```

####  **De novo assembly**
Here we use `spades` [@spades] to assemble host free reads to obtain consesus contigs and haplotype assembly.This script used paired end reads stored as fastq files. If the reads are not of this format, one may consider changing acordingly (by editing spades_assemble.sh). Among other outputs, the script produces consensus contigs and conservative regions of diploid genome.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/MNGS2/assembly"
bash "$HOME/MNGS2/scripts/denovo.sh"
```
To evaluate and assess the assembly, we use `quast`. This will provide a summary of the metagenome assembly, including but not limited to N50, N75, L50, L75, GC percentage, number of contigs with size greater than 500bp (Only assesses the consensus, similar procedure can be used to assess other outputs).
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
bash "$HOME/MNGS2/scripts/assess_assemble.sh"
```
<!--Alternative ways of doing this are also available including the used of `idba_ud` metagenome assembler.-->

### **Taxonomic identification**
Here we use `kaiju` metagenome classifier [@kaiju].
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/MNGS2/taxonomy"
bash "$HOME/MNGS2/scripts/taxonomy.sh"
```
<!--### **Molecular evolutionary analysis**{#evol} -->

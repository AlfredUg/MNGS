# MNGS Training

## **Introduction**
This is a follow up on the series of bioinformatics training sessions for the EAST African Network for Bioinformatics Trainin (EANBiT) held at KEMRI wellcome Trust, Kilifi Campus. This session focuses on the main procedure of  shotgun metagenome assembly and gene prediction.

### **Getting started**
Create a directory named "MNGS2" in your home directory. This is where the analysis will be carried out. As we go on, we shall be creating sub-directories accordingly to ensure that intermediate outputs can be refered to easily whenever needed for downstream analyses. 
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir MNGS2
```
Make a sub-directory `data` in MNGS2 and download the data. 
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir MNGS2/data
cd data
wget https://mgexamples.s3.climb.ac.uk/HMP_GUT_SRS052697.25M.1.fastq.gz
wget https://mgexamples.s3.climb.ac.uk/HMP_GUT_SRS052697.25M.2.fastq.gz
```
To confirm that the links have been successfully created, list the contents of the folder using the `ls` command. If the data is compressed (i.e .gz files), uncompress them using `gunzip *` Check to see if the data has been downloaded by using `ls -lh`. What is the size of these files?

### **Quality control check**
Here we use `fastqc` , it is necessary that we store quality control files for easy reference. In MNGS2, create a sub-directory qcresults, this is where the fastqc results will be stored.
After that, do the quality checks;
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir $HOME/MNGS2/qcresults
fastqc $HOME/MNGS2/data/*.fastq -o $HOME/MNGS2/qcresults
```
Once this is done, navigate to `qcresults` and view the ".html" files. This report can be used to assess quality distribution, length distribution, GC-content, nucleotide distribution. This informs downstream quality and adapter trimming.

### **Data quality trimming**
We use `trim_galore` for quality and adapter trimming. Depending on the qc results, it would be necessary to change some parameters used here accordingly. This script clips off the first 16 bases of the reads from the 5' end. In addition, it removes bases with phred quality less than 25 on the 3' end of the reads. We need to store quality trimmed reads, as such, in WGS directory, make a sub directory `trimmed`.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/trimmed"
trim_galore -q 25 -l 75 --dont_gzip --clip_R1 16 --clip_R2 16 --paired ../data/HMP_GUT_SRS052697.25M.1.fastq ../data/HMP_GUT_SRS052697.25M.2.fastq -o .
```

If all goes well, trimmed reads will be available in trimmed. You may consider looking at the trimmed reads using `fastqc` to check the improvement made by `trim_galore`.

#### **Remove host sequences**
There is need to get rid of host DNA sequences that could have contaminated the data earlier in the sampling and extraction stage. This is done by mapping the reads to host reference genome and picking the unmapped sequences. At this point, we create a directory to store host-free sequence data. `Samtools` and `bedtools` are the key tools used for this purpose.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir hostfree
cd hostfree
bowtie2 -x ../hostgenome/host_DB -1 ../trimmed/HMP_GUT_SRS052697.25M.1_trim.fastq -2 ../trimmed/HMP_GUT_SRS052697.25M.2_trim.fastq --threads 20 -S reads_mapped_and_unmapped.sam
```

Extract the unmapped reads, sort them and create a corresponding bam file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools view -bS reads_mapped_and_unmapped.sam | samtools view -b -f 12 -F 256 | samtools sort -n > reads_unmapped.bam
```
Convert bam files to fastq
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
bedtools bamtofastq -i reads_unmapped.bam -fq read_1.fastq -fq2 read_2.fastq
```
Now at this point we very much hope that we dont have any host sequences sequences.

####  **De novo assembly**
Here we use `IDBA-UD` to assemble host free reads to obtain consesus contig.This program takes on a single input. As such we need to merge the forward and backward reads before asssemly.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir $HOME/MNGS2/denovo
cd denovo
fq2fa --merge  ../hostfree/read_1.fastq ../hostfree/read_2.fastq  ../hostfree/reads_12.fa
idba_ud -r ../hostfree/reads_12.fa --num_threads 20 -o .
```

To evaluate and assess the assembly, we use `quast`. This will provide a summary of the metagenome assembly, including but not limited to N50, N75, L50, L75, GC percentage, number of contigs with size greater than 500bp (Only assesses the consensus, similar procedure can be used to assess other outputs).
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
quast.py -t 4 -f --meta contigs.fa -o .
```

### **Gene predicition**
For gene prediction proceeds using prodigal software. In this example, the output includes predicted genes, coded protein sequences and an annotation file. The annotation file has features, corresponding length and location on the genome among others. See full description of  at [GFF format]("http://genome.ucsc.edu/FAQ/FAQformat.html#format3")

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir genes
cd genes
prodigal -a genes.fa -d proteins.fa -i ../denovo/contigs.fa -f gff -p meta > annotation.gff
```

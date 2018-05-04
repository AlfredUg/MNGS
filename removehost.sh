indir='path/to/input/quality/trimmed'
hostgen='path/to/hostgenome' #path to host genome
suffix='_R1_qtrim.1.fastq' #change this accordingly
fmt='fastq'
hostfree='path/to/hostfree'

for sample in $(ls $indir/*suffix); do
        sname=$(basename $sample $suffix)
        bowtie2 -x $hostgen/host_DB -1 $indir/${sname}_R1*.fastq -2 $indir/${sname}_R2*.fastq --threads 20 -S $hostfree/${sname}_mapped_and_unmapped.sam
        samtools view -bS $hostfree/${sname}_mapped_and_unmapped.sam | samtools view -b -f 12 -F 256 | samtools sort -n > $hostfree/${sname}_mapped_unmapped_sorted.bam
        rm -f $hostfree/${sname}_mapped_and_unmapped.sam #remove sam files to save space
        bedtools bamtofastq -i $hostfree/${sname}_mapped_unmapped_sorted.bam -fq $hostfree/${sname}_R1.fastq -fq2 $hostfree/${sname}_R2.fastq
        rm -f $hostfree/${sname}_mapped_unmapped_sorted.bam ; #remove mapped bam files since we have the reads as fastq files
done

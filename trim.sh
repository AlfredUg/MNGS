fastqdir='/home/student9/muiiMNGS2/fastqdir'
trimdir='/home/student9/muiiMNGS2/trimmed'
suffix='_R1_val_1.fq'

for sample in $(ls $fastqdir/*R1*.fq); do
	bname=$(basename $sample $suffix)
	trim_galore -q 25 --length 75 --dont_gzip --clip_R1 16 --clip_R2 16 --paired $fastqdir/${bname}_R1*.fq $fastqdir/${bname}_R2*.fq -o $trimdir;
done

fastqdir='/home/student9/muiiMNGS2/fastqdir'
refdir='/home/student9/muiiMNGS2/refs'
alignment='home/student9/muiiMNGS2/alignment'
suffix='_R1_val_1.fq'

for sample in $(ls $fastqdir/*R1*.fq); do
	p=$(basename $sample $suffix)
	mkdir ${p}_align
	for ref in $(ls $refdir); do
		refname=$(echo $ref)
		bowtie2 -x $refdir/${refname}/${refname} -1 $fastqdir/${p}_R1*.{fmt} -2 $fastqdir/${p}_R2*.{fmt} -S ${p}_align/${p}_${refname}.sam ;
	done
done

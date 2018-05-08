#This script performs de novo assembly of any number of samples using dipspades
#It assumes that reads paired end and are quality trimmed using trim_galore,so feel free to change accordingly if this does not apply to you
indir='$HOME/MNGS2/trimmed' #path to reads to be assembled
outdir='$HOME/MNGS2/assembly' #path we to store the assembly results
suffix='_R1_val_1.fq' #naming extensions, _left.fastq, _right.fq, blah,blah,blah ...
fmt='fq'

for sample in $(ls $indir/*${suffix}); do
	echo $sample
	bname=$(basename $sample $suffix)
	mkdir $outdir/${bname}
	dipspades.py -1 $indir/${bname}_R1*.${fmt} -2 $indir/${bname}_R2*.${fmt} -o $outdir/${bname} ;
done

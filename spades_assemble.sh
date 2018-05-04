inputdir='$HOME/muiiMNGS2/trimmed'
outputdir='$HOME/muiiMNGS2/assembly'
suffix='_R1_val_1.fq' #this depends on the naming, here we assume the quality trimming was done using trim_galore, feel free to change if it is otherwise for you case 
fmt='fq'

for sample in $(ls $inputdir/*$suffix) do;
	bname=$(basename $sample $suffix)
	mkdir $outputdir/${bname}
	dipspades.py -1 $inputdir/${bname}_R1*.${fmt} -2 $inputdir/${bname}_R2*.${fmt} -o $outputdir/${bname} ;	
done


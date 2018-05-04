for ref in $(ls *.fa); do
	refname=$(echo $ref | cut -f 1,2 -d '_')
	mkdir $refname
	cp $ref $refname
	bowtie2-build -f $refname/$ref $refname/$refname;
done

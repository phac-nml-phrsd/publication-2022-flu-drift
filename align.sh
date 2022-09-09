
VIRUS=$1

OPTS='-maxiters 1 -diags -sv -distance1 kbit20_3'


file_in="seqs/"
file_in+=$VIRUS
file_in+="_cleaned.fasta"

file_out="seqs/"
file_out+=$VIRUS
file_out+="_cleaned_align.fasta"


echo $file_in
echo $file_out

muscle $OPTS -in $file_in -out $file_out 

echo -e "\n=== $VIRUS alignment completed ===\n"


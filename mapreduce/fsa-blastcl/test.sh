#for i in `seq 0 $[$1-1]`; do
#echo "$[2+($i*$2)]"
#	head -n "$[2 + ($i * $2)]" /mnt/laptop/blastdb/nr | tail -n2 > test.fasta
./formatdb ../../test/mapreduce/fsa-blastcl/test.fasta
./blast-debug -d /mnt/laptop/blastdb/nr -i ../../test/mapreduce/fsa-blastcl/test.fasta #>> paultestout
#	cat test.fasta;
#	echo "$[2+($i*$2)]"
#done
echo "Done"

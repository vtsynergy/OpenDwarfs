#!/bin/bash
echo "queryr1" >> runtime.txt
for ((i = 32; i <= 256; i = i + i))
do
	for ((j = 72; j <= 112; j = j + 2))
	do
		echo "./blast  -B "$j" -T "$i" -i ../query/nrqueryr1 -d ../database/nr07 > ./log/q1r"$j"_"$i"_1.txt"
		./blast  -B $j -T $i -i ../query/nrqueryr1 -d ../database/nr07 > ./log/q1r"$j"_"$i"_1.txt
#		echo "./blast  -B "$j" -T "$i" -i ../query/nrqueryr1 -d ../database/nr07 > ./log/q1r"$j"_"$i"_2.txt"
#		./blast  -B $j -T $i -i ../query/nrqueryr1 -d ../database/nr07 > ./log/q1r"$j"_"$i"_2.txt
#		echo "./blast  -B "$j" -T "$i" -i ../query/nrqueryr1 -d ../database/nr07 > ./log/q1r"$j"_"$i"_3.txt"
#		./blast  -B $j -T $i -i ../query/nrqueryr1 -d ../database/nr07 > ./log/q1r"$j"_"$i"_3.txt
	done
done


#!/bin/bash
echo "queryr0" >> runtime.txt
for ((i = 32; i <= 256; i = i + i))
do
	for ((j = 56; j <= 112; j = j + 2))
	do
		echo "./blast  -B "$j" -T "$i" -i ../query/nrqueryr0 -d ../database/nr07  > ./log/q0r"$j"_"$i"_1.txt"
		./blast  -B $j -T $i -i ../query/nrqueryr0 -d ../database/nr07 > ./log/q0r"$j"_"$i"_1.txt
#		echo "./blast  -B "$j" -T "$i" -i ../query/nrqueryr0 -d ../database/nr07  > ./log/q0r"$j"_"$i"_2.txt"
#		./blast  -B $j -T $i -i ../query/nrqueryr0 -d ../database/nr07 > ./log/q0r"$j"_"$i"_2.txt
#		echo "./blast  -B "$j" -T "$i" -i ../query/nrqueryr0 -d ../database/nr07  > ./log/q0r"$j"_"$i"_3.txt"
#		./blast  -B $j -T $i -i ../query/nrqueryr0 -d ../database/nr07 > ./log/q0r"$j"_"$i"_3.txt
	done
done



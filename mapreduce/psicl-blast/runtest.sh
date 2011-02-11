#!/bin/bash

#NAME=1K2
#for ((i = 192; i <= 192; i = i + 64))
#do
#	for ((j = 112; j <= 140; j = j + 2))
#	do
#		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_1.txt"
#		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_1.txt
#		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_2.txt"
#		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_2.txt
#		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_3.txt"
#		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_3.txt
#	done
#done


NAME=2K2
for ((i = 192; i <= 192; i = i + 64))
do
	for ((j = 112; j <= 140; j = j + 2))
	do
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_1.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_1.txt
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_2.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_2.txt
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_3.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_3.txt
	done
done

NAME=3K2
for ((i = 192; i <= 192; i = i + 64))
do
	for ((j = 112; j <= 140; j = j + 2))
	do
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_1.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_1.txt
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_2.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_2.txt
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_3.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_3.txt
	done
done

NAME=4K2
for ((i = 192; i <= 192; i = i + 64))
do
	for ((j = 112; j <= 140; j = j + 2))
	do
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_1.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_1.txt
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_2.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_2.txt
		echo "./blast  -B "$j" -T "$i" -i ../query/query"$NAME" -d ../database/nr07  > ./log/q"$NAME"_"$j"_"$i"_3.txt"
		./blast  -B $j -T $i -i ../query/query"$NAME" -d ../database/nr07 > ./log/q"$NAME"_"$j"_"$i"_3.txt
	done
done


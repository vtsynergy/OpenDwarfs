echo "Performance Test of 13 Dwarfs"
#defination
DWAFT_PATH=.
OUTPUT_PATH=data
OUTPUT_FILE=
ARG=
INPUT_FILE=
EXC_FILE=
INPUT_PATH=
PACK_NAME=
TEST_APPS=(gemnoui needle srad lud kmeans samplecl crc cfd bsort oesort bfs csr tdm astar clfft)
PRO_SIZE=(128 256 512 1024 2048 4096)
if [ $# != 1 ]
then 
	echo "Usage: autorun.sh package name"
	exit
else
	PACK_NAME=$1
	echo "package name: $PACK_NAME"
fi

echo "Test Begin..."
#makedir
if [ ! -d "$DWAFT_PATH/$OUTPUT_PATH" ]; then
	mkdir $DWAFT_PATH/$OUTPUT_PATH
else
#clean data
rm $DWAFT_PATH/$OUTPUT_PATH/*
fi


#swat
EXC_FILE=swat
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
INPUT_PATH="$DWAFT_PATH/test/dynamic-programming/swat"
INPUT_FILE1=("$INPUT_PATH/query1K1" \
	"$INPUT_PATH/query2K1" \
	"$INPUT_PATH/query3K1" \
	"$INPUT_PATH/query4K1" \
	"$INPUT_PATH/query5K1")

INPUT_FILE2=("$INPUT_PATH/sampledb1K1" \
	"$INPUT_PATH/sampledb2K1" \
	"$INPUT_PATH/sampledb3K1" \
	"$INPUT_PATH/sampledb4K1" \
	"$INPUT_PATH/sampledb5K1")
for ((  i = 0 ;  i < 5;  i++  ))
do
        echo -e "\n$DWAFT_PATH/$EXC_FILE  ${INPUT_FILE1[$i]} ${INPUT_FILE2[$i]}" >> $OUTPUT_FILE
        $DWAFT_PATH/$EXC_FILE  ${INPUT_FILE1[$i]} ${INPUT_FILE2[$i]}>> $OUTPUT_FILE
done 


#gemnoui
EXC_FILE=gemnoui
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
INPUT_PATH="$DWAFT_PATH/test/n-body-methods/gem"
INPUT_FILE=("$INPUT_PATH/Mb.Hhelix.bondi"\
	"$INPUT_PATH/1uwo_A"\
	"$INPUT_PATH/nucleosome")
ARG="80 1 0"
date > $OUTPUT_FILE
for ((  i = 0 ;  i < 3;  i++  ))
do
	echo -e "\n$DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]} $ARG" >> $OUTPUT_FILE
	$DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]} $ARG >> $OUTPUT_FILE
done

#needle
EXC_FILE=needle
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
for ((  i = 0 ;  i < 6;  i++  ))
do
	echo -e "\n$DWAFT_PATH/$EXC_FILE ${PRO_SIZE[$i]} 10" >> $OUTPUT_FILE
	$DWAFT_PATH/$EXC_FILE ${PRO_SIZE[$i]} 10 >> $OUTPUT_FILE
done

#srad
EXC_FILE=srad
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
for ((  i = 0 ;  i < 6;  i++  ))
do
	ARG=`expr ${PRO_SIZE[$i]} / 2 - 1`
	echo -e "\n$DWAFT_PATH/$EXC_FILE ${PRO_SIZE[$i]} ${PRO_SIZE[$i]} 0 $ARG 0 $ARG 0.5 2" >> $OUTPUT_FILE
	$DWAFT_PATH/$EXC_FILE ${PRO_SIZE[$i]} ${PRO_SIZE[$i]} 0 $ARG 0 $ARG 0.5 2 >> $OUTPUT_FILE 
done

#lud
EXC_FILE=lud
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
INPUT_PATH="$DWAFT_PATH/test/dense-linear-algebra/lud"
INPUT_FILE=("$INPUT_PATH/64.dat"\
	   "$INPUT_PATH/256.dat"\
	   "$INPUT_PATH/512.dat"\ 
	   "$INPUT_PATH/2048.dat")
for ((  i = 0 ;  i < 4;  i++  ))
do
	echo -e "\n$DWAFT_PATH/$EXC_FILE  -i ${INPUT_FILE[$i]}" >> $OUTPUT_FILE
	$DWAFT_PATH/$EXC_FILE  -i ${INPUT_FILE[$i]} >> $OUTPUT_FILE
done

#kmeans
EXC_FILE=kmeans
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
INPUT_PATH="$DWAFT_PATH/test/dense-linear-algebra/kmeans"
INPUT_FILE=("$INPUT_PATH/100" \
	"$INPUT_PATH/204800.txt" \
	"$INPUT_PATH/819200.txt" \
	"$INPUT_PATH/kdd_cup")
for ((  i = 0 ;  i < 4;  i++  ))
do
        echo -e "\n$DWAFT_PATH/$EXC_FILE  -i ${INPUT_FILE[$i]}" >> $OUTPUT_FILE
        $DWAFT_PATH/$EXC_FILE  -i ${INPUT_FILE[$i]} >> $OUTPUT_FILE
done

#crc
EXC_FILE=crc
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
$DWAFT_PATH/$EXC_FILE >> $OUTPUT_FILE

#cfd
EXC_FILE=cfd
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
INPUT_PATH="$DWAFT_PATH/test/unstructured-grids/cfd"
INPUT_FILE=("$INPUT_PATH/fvcorr.domn.097K" \
	"$INPUT_PATH/fvcorr.domn.193K" \
	"$INPUT_PATH/missile.domn.0.2M")
for ((  i = 0 ;  i < 3;  i++  ))
do
        echo -e "\n$DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]}" >> $OUTPUT_FILE
        $DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]} >> $OUTPUT_FILE
done

#bsort
EXC_FILE=bsort
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
for ((  i = 0 ;  i < 6;  i++  ))
do
	echo -e "\n$DWAFT_PATH/$EXC_FILE -n ${PRO_SIZE[$i]}" >> $OUTPUT_FILE
	$DWAFT_PATH/$EXC_FILE -n ${PRO_SIZE[$i]} >> $OUTPUT_FILE
done 

#oesort
EXC_FILE=oesort
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
echo -e "\n$DWAFT_PATH/$EXC_FILE" >> $OUTPUT_FILE
$DWAFT_PATH/$EXC_FILE  >> $OUTPUT_FILE

#bfs
EXC_FILE=bfs
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
INPUT_PATH="$DWAFT_PATH/test/graph-traversal/rodinia-bfs"
INPUT_FILE=("$INPUT_PATH/graph4096.txt" \
	"$INPUT_PATH/graph65536.txt")
for ((  i = 0 ;  i < 2;  i++  ))
do
        echo -e "\n$DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]}" >> $OUTPUT_FILE
        $DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]} >> $OUTPUT_FILE
done

#csr
EXC_FILE=csr
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
$DWAFT_PATH/$EXC_FILE >> $OUTPUT_FILE

#tdm
EXC_FILE=tdm
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
INPUT_PATH="$DWAFT_PATH/test/finite-state-machine/tdm"
INPUT_FILE=("$INPUT_PATH/stream-new-1.csv")
ARG="128 1 a d o"
for ((  i = 0 ;  i < 1;  i++  ))
do
        echo -e "\n$DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]} $INPUT_PATH $ARG" >> $OUTPUT_FILE
        $DWAFT_PATH/$EXC_FILE  ${INPUT_FILE[$i]} $INPUT_PATH $ARG >> $OUTPUT_FILE
done 

#astar
EXC_FILE=astar
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
echo -e "\n$DWAFT_PATH/$EXC_FILE" >> $OUTPUT_FILE
$DWAFT_PATH/$EXC_FILE >> $OUTPUT_FILE

#fft
EXC_FILE=clfft
echo "$EXC_FILE Begin..."
OUTPUT_FILE=$OUTPUT_PATH/$EXC_FILE.txt
date > $OUTPUT_FILE
for ((  i = 1 ;  i < 8;  i++  ))
do
echo -e "\n$DWAFT_PATH/$EXC_FILE --pts $i" >> $OUTPUT_FILE
$DWAFT_PATH/$EXC_FILE --pts $i >> $OUTPUT_FILE
done


#get cpu info
cat /proc/cpuinfo > $OUTPUT_PATH/cpuinfo.txt
echo "All Tests Finish"

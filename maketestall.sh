./autogen.sh
rm -Rf build
mkdir build
cd build
../configure --enable-timing=yes --with-opts="../opts-0.9.10/" CPPFLAGS="-g -I/usr/local/cuda/include"
make

rm -Rf ../testout
mkdir ../testout

#run astar test
./astar > ../testout/astar.out 2>../testout/astar.err

#run crc tests
#Random data set with no verification
./crc > ../testout/crc.out 2>../testout/crc.err
#Random seed (1000), serial version, with verification
./crc -s 1000 -v >> ../testout/crc.out 2>>../testout/crc.err
#Modified polynomial
./crc -p 7 >> ../testout/crc.out 2>> ../testout/crc.err

#run k-means test
./kmeans -o -i ../test/dense-linear-algebra/kmeans/100 >../testout/kmeans.out 2>../testout/kmeans.err

#run lud test
./lud -i ../test/dense-linear-algebra/lud/64.dat >../testout/lud.out 2>../testout/lud.err

#run nw test
./needle 2048 10 >../testout/nw.out 2>../testout/nw.err

#run swat test
#query1K1 vs sampledb1K1
./swat ../test/dynamic-programming/swat/query1K1 ../test/dynamic-programming/swat/sampledb1K1 >../testout/swat.out 2>../testout/swat.err

#query2K1 vs sampledb2K1 with modified penalties
./swat ../test/dynamic-programming/swat/query2K1 ../test/dynamic-programming/swat/sampledb2K1 10.0 0.5 14 >>../testout/swat.out 2>>../testout/swat.err

#run tdm test
./tdm 0 0 -- ../test/finite-state-machine/tdm/stream-new-1.csv ../test/finite-state-machine/tdm/streamIntervals.txt 100 0.5 r s n >../testout/tdm.out 2>../testout/tdm.err

#run bfs test
./bfs ../test/graph-traversal/bfs/medium.txt >../testout/bfs.out 2>../testout/bfs.err

#run gem test
#run small structure
./gemnoui ../test/n-body-methods/gem/Mb.Hhelix.bondi 80 1 0 -dumpVertices ../testout/gem.phi.out >../testout/gem.out 2>../testout/gem.err
#optional large test
./gemnoui ../test/n-body-methods/gem/nucleosome 80 1 0 -dumpVertices ../testout/gemNuc.phi.out >../testout/gemNuc.out 2>../testout/gemNuc.err

#run samplecl test
./scl >../testout/scl.out 2>../testout/scl.err

#run csr test
./csr  >../testout/csr.out 2>../testout/csr.err

#run fft test
./clfft --pts 1 >../testout/fft.out 2>../testout/fft.err

#run srad test
./srad 256 256 0 127 0 127 0.5 2 >../testout/srad.out 2>../testout/srad.err

#run cfd test
./cfd ../test/unstructured-grids/cfd/fvcorr.domn.097K >../testout/cfd.out 2>../testout/cfd.err

cd ..
ls -lh testout

# Benchmark table shared by sycl/run_all.sh and sycl/tools/validate.sh (sourced; needs HERE=<repo>/sycl and REPO).
# Inputs and sizes are those of opencl/run_all.sh (medium = the sizes of scripts/perf_final.sh for most apps).
# crc/spmv inputs are the generated files in data/inputs (override with DATA=).
T="$REPO/opencl/test"; DATA="${DATA:-$REPO/data/inputs}"
declare -A DIR OCLDIR BINNAME ARGS_small ARGS_medium CMP
# name | directory (same layout as opencl/) | binary | small args | medium args | comparison with OpenCL (tools/validate.sh):
#   "--type exact|f32 ..." read-back streams; "text REGEX" matching output text; "block SED" output block
reg() { DIR[$1]=$2; OCLDIR[$1]=$2; BINNAME[$1]=$3; ARGS_small[$1]="$4"; ARGS_medium[$1]="$5"; CMP[$1]="$6"; }
reg bfs     graph-traversal/bfs          bfs     "$T/graph-traversal/bfs/graph4096.txt" "$T/graph-traversal/bfs/graph65536.txt" "--type exact"
reg spmv    sparse-linear-algebra/SPMV   csr     "-i $DATA/csr_small.dat -a" "-i $DATA/csr_medium.dat -a" "--type f32 --rtol 1e-6 --atol 1e-6"
reg crc     combinational-logic/crc      crc     "-i $DATA/crc_small.dat -a" "-i $DATA/crc_medium.dat -a" "--type exact"
reg nw      dynamic-programming/nw       needle  "1024 10" "4096 10" "--type exact"
reg srad    structured-grids/srad        srad    "512 512 0 127 0 127 0.5 2" "2048 2048 0 127 0 127 0.5 2" "--type f32 --rtol 1e-5 --atol 1e-6"
reg lud     dense-linear-algebra/lud     lud     "-i $T/dense-linear-algebra/lud/256.dat -v" "-i $T/dense-linear-algebra/lud/2048.dat" "--type f32 --rtol 1e-4 --atol 1e-4"
reg cfd     unstructured-grids/cfd       cfd     "$T/unstructured-grids/cfd/fvcorr.domn.097K" "$T/unstructured-grids/cfd/fvcorr.domn.193K" "--type f32 --rtol 1e-4 --atol 1e-5"
reg kmeans  dense-linear-algebra/kmeans  kmeans  "-i $T/dense-linear-algebra/kmeans/204800.txt -r -o" "-i $T/dense-linear-algebra/kmeans/819200.txt -r -o" "--type exact"
reg nqueens branch-and-bound/nqueens     nqueens "12" "14" "text [0-9]+-queen has [0-9]+ solutions \([0-9]+ unique\)"
reg hmm     graphical-models/hmm         bwa_hmm "-n 60 -v n" "-n 120 -v n" "--type f32 --rtol 1e-4 --atol 1e-6"
reg swat    dynamic-programming/swat     swat    "$T/dynamic-programming/swat/query1K1 $T/dynamic-programming/swat/sampledb1K1" "$T/dynamic-programming/swat/query3K1 $T/dynamic-programming/swat/sampledb3K1" "block s/^Done\\. *//;/^=====/,/^Max position/p"
reg tdm     finite-state-machine/tdm     tdm     "$T/finite-state-machine/tdm/sim-64-size50.csv $T/finite-state-machine/tdm/ivl.txt $T/finite-state-machine/tdm/30-episodes.txt 128" "$T/finite-state-machine/tdm/sim-64-size500.csv $T/finite-state-machine/tdm/ivl.txt $T/finite-state-machine/tdm/30-episodes.txt 128" "--type exact"
reg gem     n-body-methods/gem           gemnoui "$T/n-body-methods/gem/1uwo_A 80 1 0" "$T/n-body-methods/gem/Mb.Hhelix.bondi 80 1 0" "--type f32 --rtol 1e-4 --atol 1e-5"
reg fft     spectral-methods/fft         clfft   "--pts 1 --passes 5" "--pts 3 --passes 5" "--type f32 --rtol 1e-4 --atol 1e-3"
APPS="nqueens crc lud kmeans nw swat tdm bfs hmm gem spmv srad cfd fft"

#ifndef _constants_
#define _constants_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

#define constants_minSignedInt -2147483647
#define constants_maxSignedInt 2147483647
#define constants_max_int2 32767
#define constants_gappedExtensionDummyValue -666
#define constants_sentinalScore -9999
#define constants_initialAllocUngappedExtensions 10000
#define constants_initialAllocAlignments 5000
#define constants_initialAllocGoodAlignments 1000
#define constants_initialAllocFinalAlignments 500
#define constants_initialAllocCodewordQueryPositions 5
#define constants_initialTracebackAlloc 10000
#define constants_volumeMaxSize 2000000000
#define constants_maxNumVolumes 100
#define constants_maximumTracebackSize 10000000
#define constants_unpackRegionExtend 1000
#define constants_initialAllocUnpackRegions 10000
#define constants_databaseVersion 3
#define constants_initialSequenceData 100000

extern float Robinson_prob[];
extern float Nucleotide_prob[];

#ifdef __MY_EXTERN_C__
}
#endif


#endif

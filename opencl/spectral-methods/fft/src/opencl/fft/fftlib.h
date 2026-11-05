#ifndef FFTLIB_H
#define FFTLIB_H

#include "OptionParser.h"

extern cl_device_id fftDev;
extern cl_context fftCtx;
extern cl_command_queue fftQueue;
extern Event fftEvent, ifftEvent, chkEvent;

struct cplxflt {
    float x;
    float y;
};

struct cplxdbl {
    double x;
    double y;
};
void createKernelWithSource();
void init(OptionParser& op, bool dp, int fftn);
void init2(OptionParser& op, bool dp, int fftn1, int fftn2);
void forward(void* work, void* temp, int n_ffts, int fftn);
void forward2(void* work, void* temp,  int n_ffts, int fftn1, int fftn2);
void inverse(void* work, int n_ffts);
int check(void* work, void* check, int half_n_ffts, int half_n_cmplx);
void allocDeviceBuffer(void** bufferp, unsigned long bytes);
void freeDeviceBuffer(void* buffer);
void allocHostBuffer(void** bufp, unsigned long bytes); 
void freeHostBuffer(void* buf);
void copyToDevice(void* to_device, void* from_host, unsigned long bytes);
void copyFromDevice(void* to_host, void* from_device, unsigned long bytes);
void getLocalDimension(size_t &localsz, size_t &globalsz, int fftn1, int fftn2);
void getGlobalDimension(size_t  &localsz, size_t &globalsz, int BS, int n_g);
void getLocalRadix(unsigned int n, unsigned int *radixArray, unsigned int *numRadices, unsigned int maxRadix);
void getGlobalRadix(int n, int *radix, int *R1, int *R2, int *numRadices);
void setGlobalOption(string &arg, int fftn1, int fftn2);
#endif


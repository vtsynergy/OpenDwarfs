extern int printf(constant char *format, ...);
__kernel void sort(
    __global unsigned int *input,
    unsigned int range,
    unsigned int i
){
    unsigned int tid= get_global_id(0);
    unsigned int ixj= tid^range;
    unsigned int left = input[tid];
    unsigned int right = input[ixj];
    unsigned int tmp;

    if(ixj>tid){
       if((tid&i)==0){
          if(left>right){
             tmp=left;
             left=right;
             right=tmp;
          }
       }else
          if(left<right){
             tmp=left;
             left=right;
             right=tmp;
          }
        input[tid] =left;
        input[ixj] =right;
    }
}
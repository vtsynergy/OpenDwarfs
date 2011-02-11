#define SHARED_SIZE_LIMIT 512


__kernel void bitonicSortShared(
    __global unsigned int *d_DstKey,
    __global unsigned int *d_DstVal,
    __global unsigned int *d_SrcKey,
    __global unsigned int *d_SrcVal,
    unsigned int arrayLength,
    unsigned int dir
){


    __local unsigned int s_key[SHARED_SIZE_LIMIT];
    __local unsigned int s_val[SHARED_SIZE_LIMIT];
    unsigned int size,t;
    unsigned int stride,pos,ddd;


    d_SrcKey += get_group_id(0)  * SHARED_SIZE_LIMIT + get_local_id(0);
    d_SrcVal += get_group_id(0)  * SHARED_SIZE_LIMIT + get_local_id(0);
    d_DstKey += get_group_id(0)  * SHARED_SIZE_LIMIT + get_local_id(0);
    d_DstVal += get_group_id(0)  * SHARED_SIZE_LIMIT + get_local_id(0);
    s_key[get_local_id(0) +                       0] = d_SrcKey[                      0];
    s_val[get_local_id(0) +                       0] = d_SrcVal[                      0];
    s_key[get_local_id(0) + (SHARED_SIZE_LIMIT / 2)] = d_SrcKey[(SHARED_SIZE_LIMIT / 2)];
    s_val[get_local_id(0) + (SHARED_SIZE_LIMIT / 2)] = d_SrcVal[(SHARED_SIZE_LIMIT / 2)];

    for(size = 2; size <arrayLength; size <<= 1){
        ddd = dir ^ ( (get_local_id(0) & (size / 2)) != 0 );
        for(stride = size / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));

            if( (s_key[pos +      0] > s_key[pos + stride]) == ddd ){
               t = s_key[pos +      0]; s_key[pos +      0] = s_key[pos + stride]; s_key[pos + stride] = t;
               t = s_val[pos +      0]; s_val[pos +      0] = s_val[pos + stride]; s_val[pos + stride] = t;
            }

        }
    }


    {
        for(stride = arrayLength / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            if( (s_key[pos +      0] > s_key[pos + stride]) == dir ){
               t = s_key[pos +      0]; s_key[pos +      0] = s_key[pos + stride]; s_key[pos + stride] = t;
               t = s_val[pos +      0]; s_val[pos +      0] = s_val[pos + stride]; s_val[pos + stride] = t;
            }
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    d_DstKey[                      0] = s_key[get_local_id(0) +                       0];
    d_DstVal[                      0] = s_val[get_local_id(0) +                       0];
    d_DstKey[(SHARED_SIZE_LIMIT / 2)] = s_key[get_local_id(0) + (SHARED_SIZE_LIMIT / 2)];
    d_DstVal[(SHARED_SIZE_LIMIT / 2)] = s_val[get_local_id(0) + (SHARED_SIZE_LIMIT / 2)];
}
__kernel void precompute(__global unsigned char* g_table,
                    const unsigned char crc)
{
    unsigned int tid = get_global_id(0);
    unsigned char num = tid;
    unsigned char crcCalc = 0x0;

    for(unsigned int k = 0; k < 8; k++)
    {
        //If the k-th bit is 1
        if((num >> (7-k)) % 2 == 1)
        {
            num ^= crc >> (k + 1);
            crcCalc ^= crc << (7-k);
        }
    }

    g_table[tid] = crcCalc;
}

__kernel void compute(__global unsigned char* g_num,
                    __global unsigned char* g_table,
                    __global unsigned char* g_answer,
                    const unsigned int num_size)
{
    unsigned int tid = get_global_id(0);
    if(tid < num_size)
    {
        unsigned char loc = g_num[tid];
        for(int i = 0; i < num_size - tid; i++)
        {
            loc = g_table[loc];
        }
        g_answer[tid] = loc;
    }
}

__kernel void reduce(__global unsigned char* g_answer,
                    __global unsigned char* g_reducedSet,
                  const unsigned int num_size)
{
    unsigned int gid = get_group_id(0);
    unsigned int lid = get_local_id(0);
    unsigned int tid = get_global_id(0);
    if(tid < num_size)
    {
        int mod = 2;
        int remainder = 1;
        for(mod = 2, remainder = 1; remainder <= get_local_size(0); mod *= 2, remainder *= 2)
        {
            if(lid % mod == remainder)
            {
                g_answer[tid - remainder] ^= g_answer[tid];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        if(lid == 0)
            g_reducedSet[gid] = g_answer[tid];
    }
}

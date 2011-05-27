__kernel void compute(__global unsigned char* g_num,
                    __global unsigned char* g_tables,
			const unsigned int numTables,
                    __global unsigned char* g_answer,
                    const unsigned int num_size)
{
    unsigned int tid = get_global_id(0);
    if(tid < num_size)
    {
		int val = num_size - tid;
		int i = 0;
		unsigned char temp = g_num[tid];
		for(i = 0; i < numTables; i++)
		{
			if((val >> i) % 2 == 1)
				temp = g_tables[i * 256 + temp];
		}
		g_answer[tid] = temp;
    }
}

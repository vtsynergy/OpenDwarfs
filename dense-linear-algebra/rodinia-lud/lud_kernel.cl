#define BLOCK_SIZE 16

__kernel void 
lud_diagonal(__global float *m, int matrix_dim, int offset)
{
  int i,j;
  __local float shadow[BLOCK_SIZE][BLOCK_SIZE];

  int array_offset = offset*matrix_dim+offset;
  for(i=0; i < BLOCK_SIZE; i++){
    shadow[i][get_local_id(0)]=m[array_offset+get_local_id(0)];
    array_offset += matrix_dim;
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  for(i=0; i < BLOCK_SIZE-1; i++) {

    if (get_local_id(0)>i){
      for(j=0; j < i; j++)
        shadow[get_local_id(0)][i] -= shadow[get_local_id(0)][j]*shadow[j][i];
      shadow[get_local_id(0)][i] /= shadow[i][i];

      barrier(CLK_LOCAL_MEM_FENCE);

      for(j=0; j < i+1; j++)
        shadow[i+1][get_local_id(0)] -= shadow[i+1][j]*shadow[j][get_local_id(0)];

      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }

  /* 
     The first row is not modified, it
     is no need to write it back to the
     global memory

   */
  array_offset = (offset+1)*matrix_dim+offset;
  for(i=1; i < BLOCK_SIZE; i++){
    m[array_offset+get_local_id(0)]=shadow[i][get_local_id(0)];
    array_offset += matrix_dim;
  }
}

__kernel void
lud_perimeter(__global float *m, int matrix_dim, int offset)
{
  __local float dia[BLOCK_SIZE][BLOCK_SIZE];
  __local float peri_row[BLOCK_SIZE][BLOCK_SIZE];
  __local float peri_col[BLOCK_SIZE][BLOCK_SIZE];

  int i,j, array_offset;
  int idx;

  if (get_local_id(0) < BLOCK_SIZE) {
    idx = get_local_id(0);
    
    array_offset = offset*matrix_dim+offset;
    for (i=0; i < BLOCK_SIZE/2; i++){
      dia[i][idx]=m[array_offset+idx];
      array_offset += matrix_dim;
    }
    
    array_offset = offset*matrix_dim+offset;
    for (i=0; i < BLOCK_SIZE; i++) {
      peri_row[i][idx]=m[array_offset+(get_group_id(0)+1)*BLOCK_SIZE+idx];
      array_offset += matrix_dim;
    }

  } else {
    idx = get_local_id(0)-BLOCK_SIZE;
    
    array_offset = (offset+BLOCK_SIZE/2)*matrix_dim+offset;
    for (i=BLOCK_SIZE/2; i < BLOCK_SIZE; i++){
      dia[i][idx]=m[array_offset+idx];
      array_offset += matrix_dim;
    }
    
    array_offset = (offset+(get_group_id(0)+1)*BLOCK_SIZE)*matrix_dim+offset;
    for (i=0; i < BLOCK_SIZE; i++) {
      peri_col[i][idx] = m[array_offset+idx];
      array_offset += matrix_dim;
    }
  
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  if (get_local_id(0) < BLOCK_SIZE) { //peri-row
    idx=get_local_id(0);
    for(i=1; i < BLOCK_SIZE; i++){
      for (j=0; j < i; j++)
        peri_row[i][idx]-=dia[i][j]*peri_row[j][idx];
    }
    
    array_offset = (offset+1)*matrix_dim+offset;
    for(i=1; i < BLOCK_SIZE; i++){
      m[array_offset+(get_group_id(0)+1)*BLOCK_SIZE+idx] = peri_row[i][idx];
      array_offset += matrix_dim;
    }
  } else { //peri-col
    idx=get_local_id(0) - BLOCK_SIZE;
    for(i=0; i < BLOCK_SIZE; i++){
      for(j=0; j < i; j++)
        peri_col[idx][i]-=peri_col[idx][j]*dia[j][i];
      peri_col[idx][i] /= dia[i][i];
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    
    array_offset = (offset+(get_group_id(0)+1)*BLOCK_SIZE)*matrix_dim+offset;
    for(i=0; i < BLOCK_SIZE; i++){
      m[array_offset+idx] =  peri_col[i][idx];
      array_offset += matrix_dim;
    }
  }

}

__kernel void
lud_internal(__global float *m, int matrix_dim, int offset)
{
  __local float peri_row[BLOCK_SIZE][BLOCK_SIZE];
  __local float peri_col[BLOCK_SIZE][BLOCK_SIZE];

  int i;
  float sum;

  int global_row_id = offset + (get_group_id(1)+1)*BLOCK_SIZE;
  int global_col_id = offset + (get_group_id(0)+1)*BLOCK_SIZE;

  peri_row[get_local_id(1)][get_local_id(0)] = m[(offset+get_local_id(1))*matrix_dim+global_col_id+get_local_id(0)];
  peri_col[get_local_id(1)][get_local_id(0)] = m[(global_row_id+get_local_id(1))*matrix_dim+offset+get_local_id(0)];

  barrier(CLK_LOCAL_MEM_FENCE);

  sum = 0;
  for (i=0; i < BLOCK_SIZE; i++)
    sum += peri_col[get_local_id(1)][i] * peri_row[i][get_local_id(0)];
  m[(global_row_id+get_local_id(1))*matrix_dim+global_col_id+get_local_id(0)] -= sum;


}


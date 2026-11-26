/* CUDA translation of srad_kernel.cl (mechanical: 16x16 tile per block; get_group_id->blockIdx, get_local_id->threadIdx, get_num_groups->gridDim). */
#define SRAD_TILE 16

__global__ void
srad_cuda_1(
		float* __restrict__ E_C, 
		float* __restrict__ W_C, 
		float* __restrict__ N_C, 
		float* __restrict__ S_C,
		float* __restrict__ J_cuda, 
		float* __restrict__ C_cuda, 
		int cols, 
		int rows, 
		float q0sqr
	   ) 
{

	//block id
	int bx = (int)blockIdx.x;
	int by = (int)blockIdx.y;

	//thread id
	int tx = (int)threadIdx.x;
	int ty = (int)threadIdx.y;

	//indices
	int index   = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + tx;
	int index_n = cols * SRAD_TILE * by + SRAD_TILE * bx + tx - cols;   //up
	int index_s = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * SRAD_TILE + tx; //down
	int index_w = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty - 1;  //left
	int index_e = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + SRAD_TILE;  //right

	float n, w, e, s, jc, g2, l, num, den, qsqr, c;


	if ( by == 0 ){
		index_n = SRAD_TILE * bx + tx;
	}
	else if ( by == (int)gridDim.y - 1 ){
		index_s = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ( SRAD_TILE - 1 ) + tx;
	}

	if ( bx == 0 ){
		index_w = cols * SRAD_TILE * by + cols * ty;
	}
	else if ( bx == (int)gridDim.x - 1 ){
		index_e = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + SRAD_TILE-1;
	}

	//Corresponding to [ty+1][x], [ty-1][x], etc.
	int index_yp1_x = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * (ty + 1) + tx;
	int index_ym1_x = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * (ty - 1) + tx;
	int index_y_xp1 = cols * SRAD_TILE * by + SRAD_TILE * bx + cols  * ty + (tx + 1);
	int index_y_xm1 = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + (tx - 1);


	jc = J_cuda[index];

	if ( ty == 0 && tx == 0 ){ //nw
		n  = J_cuda[index_n] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_w]  - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}	    
	else if ( ty == 0 && tx == SRAD_TILE-1 ){ //ne
		n  = J_cuda[index_n] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_y_xm1] - jc;
		e  = J_cuda[index_e] - jc;
	}
	else if ( ty == SRAD_TILE -1 && tx == SRAD_TILE - 1){ //se
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_s] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_e]  - jc;
	}
	else if ( ty == SRAD_TILE -1 && tx == 0 ){//sw
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_s] - jc;
		w  = J_cuda[index_w]  - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}
	else if ( ty == 0 ){ //n
		n  = J_cuda[index_n] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}
	else if ( tx == SRAD_TILE -1 ){ //e
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_e] - jc;
	}
	else if ( ty == SRAD_TILE -1){ //s
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_s] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}
	else if ( tx == 0 ){ //w
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_w] - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}
	else{  //the data elements which are not on the borders 
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}


	g2 = ( n * n + s * s + w * w + e * e ) / (jc * jc);

	l = ( n + s + w + e ) / jc;

	num  = (0.5f*g2) - ((1.0f/16.0f)*(l*l)) ;
	den  = 1 + (.25f*l);
	qsqr = num/(den*den);

	// diffusion coefficent (equ 33)
	den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;
	c = 1.0f / (1.0f+den) ;

	// saturate diffusion coefficent
	if (c < 0){C_cuda[index] = 0;}
	else if (c > 1) {C_cuda[index] = 1;}
	else {C_cuda[index] = c;}

	//barrier(CLK_LOCAL_MEM_FENCE);

	E_C[index] = e;
	W_C[index] = w;
	S_C[index] = s;
	N_C[index] = n;

}


__global__ void
srad_cuda_2(
		float* __restrict__ E_C, 
		float* __restrict__ W_C, 
		float* __restrict__ N_C, 
		float* __restrict__ S_C,	
		float* __restrict__ J_cuda, 
		float* __restrict__ C_cuda, 
		int cols, 
		int rows, 
		float lambda,
		float q0sqr
	   ) 
{
	//block id
	int bx = (int)blockIdx.x;
	int by = (int)blockIdx.y;

	//thread id
	int tx = (int)threadIdx.x;
	int ty = (int)threadIdx.y;


	//indices
	int index   = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + tx;
	int index_s = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * SRAD_TILE + tx;
	int index_e = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + SRAD_TILE;
	float cc, cn, cs, ce, cw, d_sum;

	//The below correspond to the index that was in localmem [ty+1][tx] and [ty][tx+1], so we load from global memory for the corresponding conditional branch
	int index_yp1_x = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * (ty+1) + tx;
	int index_y_xp1 = cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + (tx+1);

	cc = C_cuda[index];
	cn = cc;
	cw = cc;


//This kind of branching degrades fpga performance.
	if ( ty == SRAD_TILE -1 && tx == SRAD_TILE - 1){ //se
		if ( by == (int)gridDim.y - 1 )
			cs = C_cuda[cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + tx];
		else
			cs = C_cuda[index_s];

		if ( bx == (int)gridDim.x - 1 )
			ce = C_cuda[cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + tx];
		else
			ce = C_cuda[index_e];
	} 
	else if ( tx == SRAD_TILE -1 ){ //e
		cs  = C_cuda[index_yp1_x];

		if ( bx == (int)gridDim.x - 1 )
			ce = C_cuda[cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + tx];
		else
			ce = C_cuda[index_e];
	}
	else if ( ty == SRAD_TILE -1){ //s
		if ( by == (int)gridDim.y - 1 )
			cs = C_cuda[cols * SRAD_TILE * by + SRAD_TILE * bx + cols * ty + tx];
		else
			cs = C_cuda[index_s];
		ce  = C_cuda[index_y_xp1];
	}
	else{ //the data elements which are not on the borders 
		cs  = C_cuda[index_yp1_x];
		ce  = C_cuda[index_y_xp1];
	}

	// divergence (equ 58)
	d_sum = cn * N_C[index] + cs * S_C[index] + cw * W_C[index] + ce * E_C[index];

	// image update (equ 61)
	J_cuda[index] = J_cuda[index] + 0.25f * lambda * d_sum;

	//barrier(CLK_GLOBAL_MEM_FENCE);

}

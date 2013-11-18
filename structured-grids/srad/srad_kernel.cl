
	__kernel void
srad_cuda_1(
		__global float *E_C, 
		__global float *W_C, 
		__global float *N_C, 
		__global float *S_C,
		__global float * J_cuda, 
		__global float * C_cuda, 
		int cols, 
		int rows, 
		float q0sqr
	   ) 
{

	//block id
	int bx = get_group_id(0);
	int by = get_group_id(1);

	//thread id
	int tx = get_local_id(0);
	int ty = get_local_id(1);

	//indices
	int index   = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx;
	int index_n = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + tx - cols;
	int index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * BLOCK_SIZE + tx;
	int index_w = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty - 1;
	int index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE;

	float n, w, e, s, jc, g2, l, num, den, qsqr, c;


	if ( by == 0 ){
		index_n = BLOCK_SIZE * bx + tx;
	}
	else if ( by == get_num_groups(1) - 1 ){
		index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ( BLOCK_SIZE - 1 ) + tx;
	}

	if ( bx == 0 ){
		index_w = cols * BLOCK_SIZE * by + cols * ty;
	}
	else if ( bx == get_num_groups(0) - 1 ){
		index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE-1;
	}

	//Corresponding to [ty+1][x], [ty-1][x], etc.
	int index_yp1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty + 1) + tx;
	int index_ym1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty - 1) + tx;
	int index_y_xp1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx + 1);
	int index_y_xm1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx - 1);


	jc = J_cuda[index];

	if ( ty == 0 && tx == 0 ){ //nw
		n  = J_cuda[index_n] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_w]  - jc; 
		e  = J_cuda[index_y_xp1] - jc;
	}	    
	else if ( ty == 0 && tx == BLOCK_SIZE-1 ){ //ne
		n  = J_cuda[index_n] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_e] - jc;
	}
	else if ( ty == BLOCK_SIZE -1 && tx == BLOCK_SIZE - 1){ //se
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_s] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_e]  - jc;
	}
	else if ( ty == BLOCK_SIZE -1 && tx == 0 ){//sw
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
	else if ( tx == BLOCK_SIZE -1 ){ //e
		n  = J_cuda[index_ym1_x] - jc;
		s  = J_cuda[index_yp1_x] - jc;
		w  = J_cuda[index_y_xm1] - jc; 
		e  = J_cuda[index_e] - jc;
	}
	else if ( ty == BLOCK_SIZE -1){ //s
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

	__kernel void
srad_cuda_2(
		__global float *E_C, 
		__global float *W_C, 
		__global float *N_C, 
		__global float *S_C,	
		__global float * J_cuda, 
		__global float * C_cuda, 
		int cols, 
		int rows, 
		float lambda,
		float q0sqr
	   ) 
{
	//block id
	int bx = get_group_id(0);
	int by = get_group_id(1);

	//thread id
	int tx = get_local_id(0);
	int ty = get_local_id(1);

	//indices
	int index   = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx;
	int index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * BLOCK_SIZE + tx;
	int index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE;
	float cc, cn, cs, ce, cw, d_sum;

	//The below correspond to the index that was in localmem [ty+1][tx] and [ty][tx+1], so we load from global memory for the corresponding conditional branch
	int index_yp1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty+1) + tx;
	int index_y_xp1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx+1);

	cc = C_cuda[index];
	cn = cc;
	cw = cc;

	if ( ty == BLOCK_SIZE -1 && tx == BLOCK_SIZE - 1){ //se
		if ( by == get_num_groups(1) - 1 )
			cs = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			cs = C_cuda[index_s];

		if ( bx == get_num_groups(0) - 1 )
			ce = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			ce = C_cuda[index_e];
	} 
	else if ( tx == BLOCK_SIZE -1 ){ //e
		cs  = C_cuda[index_yp1_x];

		if ( bx == get_num_groups(0) - 1 )
			ce = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			ce = C_cuda[index_e];
	}
	else if ( ty == BLOCK_SIZE -1){ //s
		if ( by == get_num_groups(1) - 1 )
			cs = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
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

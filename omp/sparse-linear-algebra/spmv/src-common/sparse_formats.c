#include "../inc/sparse_formats.h"

triplet* triplet_new_array(const size_t N) {
	//dispatch on location
	return (triplet*) malloc(N * sizeof(triplet));
}

int triplet_comparator(const void *v1, const void *v2)
{
	const triplet* t1 = (triplet*) v1;
	const triplet* t2 = (triplet*) v2;

	if(t1->i < t2->i)
		return -1;
	else if(t1->i > t2->i)
		return +1;
	else if(t1->j < t2->j)
		return -1;
	else if(t1->j > t2->j)
		return +1;
	else
		return 0;
}

int unsigned_int_comparator(const void* v1, const void* v2)
{
	const unsigned int int1 = *((unsigned int*) v1);
	const unsigned int int2 = *((unsigned int*) v2);

	if(int1 < int2)
		return -1;
	else if(int1 > int2)
		return +1;
	else
		return 0;
}

void write_csr(const csr_matrix* csr,const unsigned int num_csr,const char* file_path)
{
	FILE* fp;
	int i,j;
	fp = fopen(file_path,"w");
	check(fp != NULL,"sparse_formats.write_csr() - Cannot Open File");
	fprintf(fp,"%u\n\n",num_csr);

	for(j=0; j<num_csr; j++)
	{
		fprintf(fp,"%u\n%u\n%u\n%u\n%lf\n%lf\n%lf\n",csr[j].num_rows,csr[j].num_cols,csr[j].num_nonzeros,csr[j].density_ppm,csr[j].density_perc,csr[j].nz_per_row,csr[j].stddev);

		for(i=0; i<=csr[j].num_rows; i++)
			fprintf(fp,"%u ",csr[j].Ap[i]);
		fprintf(fp,"\n");

		for(i=0; i<csr[j].num_nonzeros; i++)
			fprintf(fp,"%u ",csr[j].Aj[i]);
		fprintf(fp,"\n");

		for(i=0; i<csr[j].num_nonzeros; i++)
			fprintf(fp,"%f ",csr[j].Ax[i]);
		fprintf(fp,"\n\n");
	}

	fclose(fp);
}

csr_matrix* read_csr(unsigned int* num_csr,const char* file_path)
{
	FILE* fp;
	int i,j,read_count;
	csr_matrix* csr;

	check(num_csr != NULL,"sparse_formats.read_csr() - ptr to num_csr is NULL!");

	fp = fopen(file_path,"r");
	check(fp != NULL,"sparse_formats.read_csr() - Cannot Open Input File");

	read_count = fscanf(fp,"%u\n\n",num_csr);
	check(read_count == 1,"sparse_formats.read_csr() - Input File Corrupted! Read count for num_csr differs from 1");
	csr = malloc(sizeof(struct csr_matrix)*(*num_csr));

	for(j=0; j<*num_csr; j++)
	{
		read_count = fscanf(fp,"%u\n%u\n%u\n%u\n%lf\n%lf\n%lf\n",&(csr[j].num_rows),&(csr[j].num_cols),&(csr[j].num_nonzeros),&(csr[j].density_ppm),&(csr[j].density_perc),&(csr[j].nz_per_row),&(csr[j].stddev));
		check(read_count == 7,"sparse_formats.read_csr() - Input File Corrupted! Read count for header info differs from 7");

		read_count = 0;
		csr[j].Ap = int_new_array(csr[j].num_rows+1,"sparse_formats.read_csr() - Heap Overflow! Cannot allocate space for csr.Ap");
		for(i=0; i<=csr[j].num_rows; i++)
			read_count += fscanf(fp,"%u ",csr[j].Ap+i);
		check(read_count == (csr[j].num_rows+1),"sparse_formats.read_csr() - Input File Corrupted! Read count for Ap differs from csr[j].num_rows+1");

		read_count = 0;
		csr[j].Aj = int_new_array(csr[j].num_nonzeros,"sparse_formats.read_csr() - Heap Overflow! Cannot allocate space for csr.Aj");
		for(i=0; i<csr[j].num_nonzeros; i++)
			read_count += fscanf(fp,"%u ",csr[j].Aj+i);
		check(read_count == (csr[j].num_nonzeros),"sparse_formats.read_csr() - Input File Corrupted! Read count for Aj differs from csr[j].num_nonzeros");

		read_count = 0;
		csr[j].Ax = float_new_array(csr[j].num_nonzeros,"sparse_formats.read_csr() - Heap Overflow! Cannot allocate space for csr.Ax");
		for(i=0; i<csr[j].num_nonzeros; i++)
			read_count += fscanf(fp,"%f ",csr[j].Ax+i);
		check(read_count == (csr[j].num_nonzeros),"sparse_formats.read_csr() - Input File Corrupted! Read count for Ax differs from csr[j].num_nonzeros");
	}

	fclose(fp);
	return csr;
}

void print_timestamp(FILE* stream)
{
	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	fprintf(stream,"Current time: %s",asctime(timeinfo));
}

unsigned long gen_rand(const long LB, const long HB)
{
	int range = HB - LB + 1;
	check((HB >= 0 && LB >= 0 && range > 0),"sparse_formats.gen_rand() - Invalid Bound(s). Exiting...");
	return (rand() % range) + LB;
}

csr_matrix laplacian_5pt(const unsigned int N)
{

	csr_matrix csr;
	csr.num_rows = N*N;
	csr.num_cols = N*N;
	csr.num_nonzeros = 5*N*N - 4*N;

	csr.Ap = int_new_array(csr.num_rows+4,"sparse_formats.laplacian_5pt() - Heap Overflow! Cannot allocate space for csr.Ap");
	csr.Aj = int_new_array(csr.num_nonzeros,"sparse_formats.laplacian_5pt() - Heap Overflow! Cannot allocate space for csr.Aj");
	csr.Ax = float_new_array(csr.num_nonzeros,"sparse_formats.laplacian_5pt() - Heap Overflow! Cannot allocate space for csr.Ax");

	unsigned int nz = 0;
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int indx = 0;

	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			indx = N*i + j;

			if (i > 0){
				csr.Aj[nz] = indx - N;
				csr.Ax[nz] = -1;
				nz++;
			}

			if (j > 0){
				csr.Aj[nz] = indx - 1;
				csr.Ax[nz] = -1;
				nz++;
			}

			csr.Aj[nz] = indx;
			csr.Ax[nz] = 4;
			nz++;

			if (j < N - 1){
				csr.Aj[nz] = indx + 1;
				csr.Ax[nz] = -1;
				nz++;
			}

			if (i < N - 1){
				csr.Aj[nz] = indx + N;
				csr.Ax[nz] = -1;
				nz++;
			}

			csr.Ap[indx + 1] = nz;
		}
	}
	return csr;
}


int bin_search(const triplet* data, int size, const triplet* key)
{
	triplet* mid_triplet;
	int lo,hi,m;
	lo = 0;
	hi = size-1;
	while(lo <= hi) //binary search to determine if element exists and, if not, what is the proper index for insertion
	{
		m = lo + ((hi - lo)/2);
		if(triplet_comparator(key,&(data[m])) > 0)
			lo = m + 1;
		else if (triplet_comparator(key,&(data[m])) < 0)
			hi = m - 1;
		else
			return m;
	}
	return (-1*lo - 1);
}

coo_matrix rand_coo(const unsigned int N,const unsigned long density, FILE* log)
{
	coo_matrix coo;
	triplet tmp, *current_triplet, *mid_triplet;

	unsigned int ind;
	int m;

	coo.num_rows = N;
	coo.num_cols = N;
	coo.density_ppm = density;
	coo.num_nonzeros = (((double)(N*density))/1000000.0)*N;
	printf("NUM_nonzeros: %d\n",coo.num_nonzeros);

	coo.non_zero = triplet_new_array(coo.num_nonzeros);
	check(coo.non_zero != NULL,"sparse_formats.rand_coo_bin_insertion(): Heap Overflow - Cannot allocate memory for coo.non_zero\n");
	print_timestamp(log);
	fprintf(log,"Memory Allocated. Generating Data...\n");

	current_triplet = &(coo.non_zero[0]); //Generate random first element
	(current_triplet->i) = gen_rand(0,N-1);
	(current_triplet->j) = gen_rand(0,N-1);

	for(ind=1; ind<coo.num_nonzeros; ind++)
	{
		current_triplet = &(coo.non_zero[ind]); //element to be inserted
		(current_triplet->i) = gen_rand(0,N-1);
		(current_triplet->j) = gen_rand(0,N-1);

		m = bin_search(coo.non_zero,ind,current_triplet);
		if(m < 0)
		{
			m = -1*m - 1;
		}
		else
		{
			ind--;
			continue;
		}

		if(m < ind)
		{
			tmp = *current_triplet;
			memmove(coo.non_zero + m + 1,coo.non_zero+m,sizeof(triplet)*(ind-m));
			coo.non_zero[m] = tmp;
		}
	}

	for(ind=0; ind<coo.num_nonzeros; ind++)
	{
		current_triplet = &(coo.non_zero[ind]);
		(current_triplet->v) = 1.0 - 2.0 * (rand() / (2147483647 + 1.0));
		while((current_triplet->v) == 0.0)
			(current_triplet->v) = 1.0 - 2.0 * (rand() / (2147483647 + 1.0));
	}

	print_timestamp(log);
	fprintf(log,"Matrix Completed. Returning...\n");

	return coo;
}

void print_coo_metadata(const coo_matrix* coo, FILE* stream) {
	fprintf(stream,"\nCOO Matrix Metadata:\n\nNRows=%d\tNCols=%d\tNNZ=%d\tDensity (ppm)=%d\tDensity (fract)=%g\n\n",coo->num_rows,coo->num_cols,coo->num_nonzeros,coo->density_ppm,(double)(coo->density_ppm/1000000.0));
}

void print_csr_metadata(const csr_matrix* csr, FILE* stream) {
	fprintf(stream,"\nCSR Matrix Metadata:\n\nNRows=%lu\tNCols=%lu\tNNZ=%lu\tDensity=%lu ppm = %g%%\tAverage NZ/Row=%g\tStdDev NZ/Row=%g\n\n",csr->num_rows,csr->num_cols,csr->num_nonzeros,csr->density_ppm,csr->density_perc,csr->nz_per_row,csr->stddev);
}

void print_coo(const coo_matrix* coo, FILE* stream)
{
	unsigned int ind;
	fprintf(stream,"\nPrinting COO Matrix in COO Form:\n\nNRows=%d\nNCols=%d\nNNZ=%d\nDensity (ppm)=%d\nDensity (fract)=%g\n",coo->num_rows,coo->num_cols,coo->num_nonzeros,coo->density_ppm,(double)(coo->density_ppm/1000000.0));
	for(ind=0; ind<coo->num_nonzeros; ind++)
		fprintf(stream,"(%2d,%2d,%5.2f)\n",coo->non_zero[ind].i,coo->non_zero[ind].j,coo->non_zero[ind].v);
}

void print_coo_std(const coo_matrix* coo,FILE* stream)
{
	int ind,ind2,nz_count=0;
	float val;

	fprintf(stream,"\nPrinting COO Matrix in Standard Form:\n\nNRows=%d\nNCols=%d\nNNZ=%d\nDensity (ppm)=%d\nDensity (fract)=%g\n",coo->num_rows,coo->num_cols,coo->num_nonzeros,coo->density_ppm,(double)(coo->density_ppm/1000000.0));

	for(ind=0; ind<coo->num_rows; ind++)
	{
		fprintf(stream,"[");
		for(ind2=0; ind2<coo->num_cols; ind2++)
		{
			if(ind == coo->non_zero[nz_count].i && ind2 == coo->non_zero[nz_count].j)
				val = coo->non_zero[nz_count++].v;
			else
				val = 0.0;
			fprintf(stream,"%6.2f",val);
		}
		fprintf(stream,"]\n");
	}
}

void print_csr_arr_std(const csr_matrix* csr, const unsigned int num_csr, FILE* stream)
{
	unsigned int k;
	for(k=0; k<num_csr; k++)
		print_csr_std(&csr[k],stream);
}

void print_csr_std(const csr_matrix* csr,FILE* stream)
{
	int ind,ind2,nz_count=0,row_count=0,next_nz_row;
	float val,density;
	density = ((float)(csr->num_nonzeros))/(((float)(csr->num_rows))*((float)(csr->num_cols)));

	print_csr_metadata(csr,stream);

	while(csr->Ap[row_count+1] == nz_count)
		row_count++;

	for(ind=0; ind<csr->num_rows; ind++)
	{
		fprintf(stream,"[");
		for(ind2=0; ind2<csr->num_cols; ind2++)
		{
			if(ind == row_count && ind2 == csr->Aj[nz_count])
			{
				val = csr->Ax[nz_count++];
				while(csr->Ap[row_count+1] == nz_count)
					row_count++;
			}
			else
				val = 0.0;
			fprintf(stream,"%6.2f",val);
		}
		fprintf(stream,"]\n");
	}
	fprintf(stream,"\n");
}

csr_matrix coo_to_csr(const coo_matrix* coo,FILE* log)
{
	int ind,row_count,newline_count;

	csr_matrix csr;
	csr.num_rows = coo->num_rows;
	csr.num_cols = coo->num_cols;
	csr.num_nonzeros = coo->num_nonzeros;

	csr.Ap = int_new_array(csr.num_rows+1,"sparse_formats.coo_to_csr() - Heap Overflow! Cannot allocate space for csr.Ap");
	csr.Aj = int_new_array(csr.num_nonzeros,"sparse_formats.coo_to_csr() - Heap Overflow! Cannot allocate space for csr.Aj");
	csr.Ax = float_new_array(csr.num_nonzeros,"sparse_formats.coo_to_csr() - Heap Overflow! Cannot allocate space for csr.Ax");

	print_timestamp(log);
	fprintf(log,"Memory Allocated. Copying column indices & values...\n");

	for(ind=0; ind<coo->num_nonzeros; ind++)
	{
		csr.Ax[ind] = coo->non_zero[ind].v;
		csr.Aj[ind] = coo->non_zero[ind].j;
	}

	print_timestamp(log);
	fprintf(log,"Calculating Row Pointers...\n");

	row_count = 0;
	ind = 0;
	while(row_count <= coo->non_zero[ind].i)
		csr.Ap[row_count++] = 0;

	for(ind=1; ind<coo->num_nonzeros; ind++)
	{
		newline_count = coo->non_zero[ind].i - coo->non_zero[ind-1].i;
		while(newline_count > 0)
		{
			csr.Ap[row_count++] = ind;
			newline_count--;
		}
	}
	csr.Ap[row_count] = csr.num_nonzeros;

	print_timestamp(log);
	fprintf(log,"Conversion Complete. Returning...\n");

	return csr;
}

csr_matrix rand_csr(const unsigned int N,const unsigned int density, const double normal_stddev,unsigned long* seed,FILE* log)
{
	unsigned int i,j,nnz_ith_row,nnz,update_interval,rand_col;
	double nnz_ith_row_double,nz_error,nz_per_row_doubled,high_bound;
	int kn[128];
	float fn[128],wn[128];
	char* used_cols;
	csr_matrix csr;

	csr.num_rows = N;
	csr.num_cols = N;
	csr.density_perc = (((double)(density))/10000.0);
	csr.nz_per_row = (((double)N)*((double)density))/1000000.0;
	csr.num_nonzeros = round(csr.nz_per_row*N);
	csr.stddev = normal_stddev * csr.nz_per_row; //scale normalized standard deviation by average NZ/row

	fprintf(log,"Average NZ/Row: %-8.3f\n",csr.nz_per_row);
	fprintf(log,"Standard Deviation: %-8.3f\n",csr.stddev);
	fprintf(log,"Target Density: %u ppm = %g%%\n",density,csr.density_perc);
	fprintf(log,"Approximate NUM_nonzeros: %d\n",csr.num_nonzeros);

	csr.Ap = int_new_array(csr.num_rows+1,"rand_csr() - Heap Overflow! Cannot Allocate Space for csr.Ap");
	csr.Aj = int_new_array(csr.num_nonzeros,"rand_csr() - Heap Overflow! Cannot Allocate Space for csr.Aj");

	csr.Ap[0] = 0;
	nnz = 0;
	nz_per_row_doubled = 2*csr.nz_per_row; //limit nnz_ith_row to double the average because negative values are rounded up to 0. This
	high_bound = MINIMUM(csr.num_cols,nz_per_row_doubled); //limitation ensures the distribution will be symmetric about the mean, albeit not truly normal.
	used_cols = malloc(csr.num_cols*sizeof(char));
	check(used_cols != NULL,"rand_csr() - Heap Overflow! Cannot allocate space for used_cols");

	r4_nor_setup(kn,fn,wn);
	srand(*seed);

	update_interval = round(csr.num_rows / 10.0);
	if(!update_interval) update_interval = csr.num_rows;

	for(i=0; i<csr.num_rows; i++)
	{
		if(i % update_interval == 0) fprintf(log,"\t%d of %d (%5.1f%%) Rows Generated. Continuing...\n",i,csr.num_rows,((double)(i))/csr.num_rows*100);

		nnz_ith_row_double = r4_nor(seed,kn,fn,wn); //random, normally-distributed value for # of nz elements in ith row, NORMALIZED
		nnz_ith_row_double *= csr.stddev; //scale by standard deviation
		nnz_ith_row_double += csr.nz_per_row; //add average nz/row
		if(nnz_ith_row_double < 0)
			nnz_ith_row = 0;
		else if(nnz_ith_row_double > high_bound)
			nnz_ith_row = high_bound;
		else
			nnz_ith_row = (unsigned int) round(nnz_ith_row_double);

		csr.Ap[i+1] = csr.Ap[i] + nnz_ith_row;
		if(csr.Ap[i+1] > csr.num_nonzeros)
			csr.Aj = realloc(csr.Aj,sizeof(unsigned int)*csr.Ap[i+1]);

		for(j=0; j<csr.num_cols; j++)
			used_cols[j] = 0;

		for(j=0; j<nnz_ith_row; j++)
		{
			rand_col = gen_rand(0,csr.num_cols - 1);
			if(used_cols[rand_col])
			{
				j--;
			}
			else
			{
				csr.Aj[csr.Ap[i]+j] = rand_col;
				used_cols[rand_col] = 1;
			}
		}
		qsort((&(csr.Aj[csr.Ap[i]])),nnz_ith_row,sizeof(unsigned int),unsigned_int_comparator);
	}

	nz_error = ((double)abs((signed int)(csr.num_nonzeros - csr.Ap[csr.num_rows]))) / ((double)csr.num_nonzeros);
	if(nz_error >= .05)
		fprintf(stderr,"WARNING: Actual NNZ differs from Theoretical NNZ by %5.2f%%!\n",nz_error*100);
	csr.num_nonzeros = csr.Ap[csr.num_rows];
	fprintf(log,"Actual NUM_nonzeros: %d\n",csr.num_nonzeros);
	csr.density_perc = (((double)csr.num_nonzeros)*100.0)/((double)csr.num_cols)/((double)csr.num_rows);
	csr.density_ppm = (unsigned int)round(csr.density_perc * 10000.0);
	fprintf(log,"Actual Density: %u ppm = %g%%\n",csr.density_ppm,csr.density_perc);

	free(used_cols);
	csr.Ax = float_new_array(csr.num_nonzeros,"rand_csr() - Heap Overflow! Cannot Allocate Space for csr.Ax");
	for(i=0; i<csr.num_nonzeros; i++)
	{
		csr.Ax[i] = 1.0 - 2.0 * (rand() / (2147483647 + 1.0));
		while(csr.Ax[i] == 0.0)
			csr.Ax[i] = 1.0 - 2.0 * (rand() / (2147483647 + 1.0));
	}

	return csr;
}

void free_csr(csr_matrix* csr,const unsigned int num_csr)
{
	int k;
	for(k=0; k<num_csr; k++)
	{
		free(csr[k].Ap);
		free(csr[k].Aj);
		free(csr[k].Ax);
	}
	free(csr);
}

/***********************************************************************
 * BFS.c
 *
 * A port of the Rodinia BFS test from CUDA to OpenCL.
 * 
 * Written by: Kenneth Lee
 * Last Modified: 2/4/2010
 **********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <utility>
#define __NO_STD_VECTOR // Use cl::vector and cl::string and 
#define __NO_STD_STRING  // not STL versions, more on this later
#include <CL/cl.h>
//#include <oclUtils.h>
#include <ctime>

const char* kernelSource1 = "bfs_kernel.cl";
const char* kernelSource2 = "kernel2.cl";

unsigned int no_of_nodes;
unsigned int edge_list_size;
FILE *fp;

//Structure for Nodes in the graph
struct Node
{
	int starting;     //Index where the edges of the node start
	int no_of_edges;  //The degree of the node
};

static cl_context createContext(void)
{
     cl_device_type dType = CL_DEVICE_TYPE_GPU;

     /*
     * Have a look at the available platforms and pick either
     * the AMD one if available or a reasonable default.
     */

    cl_uint numPlatforms;
    cl_platform_id platform = NULL;
    clGetPlatformIDs(0, NULL, &numPlatforms);
    if (0 < numPlatforms) 
    {
        cl_platform_id* platforms = new cl_platform_id[numPlatforms];
        clGetPlatformIDs(numPlatforms, platforms, NULL);
        
	for (unsigned i = 0; i < numPlatforms; ++i) 
        {
            char pbuf[100];
            clGetPlatformInfo(platforms[i],
                                       CL_PLATFORM_VENDOR,
                                       sizeof(pbuf),
                                       pbuf,
                                       NULL);

            platform = platforms[i];
            if (!strcmp(pbuf, "Advanced Micro Devices, Inc.") && !strcmp(pbuf, "NVIDIA Corporation")) 
            {
                break;
            }
        }
        delete[] platforms;
    }

    /*
     * If we could find our platform, use it. Otherwise pass a NULL and get whatever the
     * implementation thinks we should be using.
     */

    cl_context_properties cps[3] = 
    {
        CL_CONTEXT_PLATFORM, 
        (cl_context_properties)platform, 
        0
    };
    /* Use NULL for backward compatibility */
    cl_context_properties* cprops = (NULL == platform) ? NULL : cps;

    cl_context context = clCreateContextFromType(cprops,
                                      dType,
                                      NULL,
                                      NULL,
                                      NULL);
    return context;
}

static const char* oclLoadProgSource(const char* filename, size_t* out_length)
{
	FILE* file;
	file = fopen(filename, "r");
	if(!file)
	{
		printf("Error when opening file: %s", filename);
		exit(1);
	}
	fseek(file, 0, SEEK_END);
	size_t length = ftell(file);
	rewind(file);
	char* source = (char*)malloc(length);
	fread(source, 1, length, file);
	fclose(file);
	*out_length = length;
	//printf("%s\n", source);
	return source;
}

void BFSGraph(int argc, char** argv);

/******************************************************************************
 * MAIN METHOD
 *****************************************************************************/
int main(int argc, char** argv)
{
    no_of_nodes = 0;
    edge_list_size = 0;
    BFSGraph(argc, argv);
    return 0;
}

/******************************************************************************
 * Apply BFS on a Graph using OpenCL
 *****************************************************************************/
void BFSGraph(int argc, char ** argv)
{
    if(argc != 2)
    {
        printf("No input file given.\n");
        exit(1);
    }
    printf("Reading File\n");
    //Read in Graph from a file
    fp = fopen(argv[1], "r");
    if(!fp)
    {
	printf("Error Reading graph file\n");
    	return;
    }

    int source = 0;
    fscanf(fp, "%d", &no_of_nodes);

    //allocate host memory
    Node* h_graph_nodes = (Node*) malloc(sizeof(Node) * no_of_nodes);
    int* h_graph_mask = (int*) malloc(sizeof(int) * no_of_nodes);
    int* h_updating_graph_mask = (int*) malloc(sizeof(int) * no_of_nodes);
    int* h_graph_visited = (int*) malloc(sizeof(int) * no_of_nodes);

    int start, edgeno;
    //initalize the memory
    for(unsigned int i = 0; i < no_of_nodes; i++)
    {
	fscanf(fp, "%d %d", &start, &edgeno);
	h_graph_nodes[i].starting = start;           //Hmmmm... Seems to me that start and edgeno are undefined...
	h_graph_nodes[i].no_of_edges = edgeno;
	h_graph_mask[i] = 0;
	h_updating_graph_mask[i] = 0;
	h_graph_visited[i] = 0;
    }

    //read the source node from the file
    fscanf(fp, "%d", &source);
    source = 0;					     //Hmmmm.... Seems that the fscanf is not really used....

    //set the source node as true in the masks
    h_graph_mask[source] = 1;
    h_graph_visited[source] = 1;

    fscanf(fp, "%d", &edge_list_size);

    int id, cost;
    int* h_graph_edges = (int*) malloc(sizeof(int) * edge_list_size);
    for(unsigned int i = 0; i < edge_list_size; i++)
    {
	fscanf(fp, "%d", &id);
	fscanf(fp, "%d", &cost);
	h_graph_edges[i] = id;
    }

    if(fp)
	fclose(fp);

    printf("Read File\n");


    //Setup OpenCL
    //Create a context to OpenCL
    cl_context GPUContext = createContext();//clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);
    
    //Get the list of GPU devices assocated with the context
    size_t parmDataBytes;
    clGetContextInfo(GPUContext, CL_CONTEXT_DEVICES, 0, NULL, &parmDataBytes);
    //printf("Size: %ld\n", parmDataBytes);
    cl_device_id* GPUDevices = (cl_device_id*)malloc(parmDataBytes);
    clGetContextInfo(GPUContext, CL_CONTEXT_DEVICES, parmDataBytes, GPUDevices, NULL);

    //Create a command-queue on the first GPU device
    cl_int err = 0;
    cl_command_queue commandQueue = clCreateCommandQueue(GPUContext, GPUDevices[0], 0, &err);//GPUDevices[0], 0, NULL);
    if(err != CL_SUCCESS)
    {
	    printf("Error creating Command Queue\n");
    }


    //Copy the Node list to device memory
    cl_mem d_graph_nodes = clCreateBuffer(GPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	    sizeof(Node) * no_of_nodes, h_graph_nodes, &err);
    //Copy the Edge List to device memory
    cl_mem d_graph_edges = clCreateBuffer(GPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	    sizeof(int) * edge_list_size, h_graph_edges, &err);
    //Copy the Mask to device memory
    cl_mem d_graph_mask = clCreateBuffer(GPUContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    sizeof(int) * no_of_nodes, h_graph_mask, &err);
    //Copy the updating graph mask to device memory
    cl_mem d_updating_graph_mask = clCreateBuffer(GPUContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    sizeof(int) * no_of_nodes, h_updating_graph_mask, &err);
    //Copy the Visited nodes to device memory
    cl_mem d_graph_visited = clCreateBuffer(GPUContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    sizeof(int) * no_of_nodes, h_graph_visited, &err);
    //Allocate memory for the result on host side
    int* h_cost = (int*) malloc(sizeof(int) * no_of_nodes);
    for(unsigned int i = 0; i < no_of_nodes; i++)
    	h_cost[i] = -1;
    h_cost[source] = 0;
    //Allocate device memory for result
    cl_mem d_cost = clCreateBuffer(GPUContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
	    sizeof(int) * no_of_nodes, h_cost, &err);
    //Make a bool to check if the execution is over
    cl_mem d_over = clCreateBuffer(GPUContext, CL_MEM_WRITE_ONLY,
	    sizeof(int), NULL, &err);

    printf("Copied Everything to GPU memory\n");

    //setup execution parameters (compile code)
    size_t szKernelLength;
    //const char* cPathAndName = shrFindFilePath(kernelSource1, argv[0]);
    const char* cSourceCL = oclLoadProgSource(kernelSource1, &szKernelLength);
   
    cl_program kernel1Program = clCreateProgramWithSource(GPUContext, 1, (const char**) &cSourceCL, &szKernelLength, &err);
    if(err != CL_SUCCESS)
	printf("ERROR COMPILING KERNEL 1 (%d)\n", err);

    // Build the program (OpenCL JIT compilation)

    err = clBuildProgram(kernel1Program, 1, GPUDevices, NULL, NULL, NULL);	
    if(err != CL_SUCCESS)
    {
	printf("Error building kernel 1(%d).\n", err);
	char programLog[1024];
   	clGetProgramBuildInfo(kernel1Program, GPUDevices[0], CL_PROGRAM_BUILD_LOG, 1024, programLog, 0);
	printf("%s\n",programLog);
    }
  
   cl_kernel kernel1 = clCreateKernel(kernel1Program, "kernel1", &err);
    if(err != CL_SUCCESS)
	printf("Error creating Kernel 1(%d).\n", err);

    /*
    //////////////////////
    //Kernel 2 Compilation
    //////////////////////
    cPathAndName = shrFindFilePath(kernelSource2, argv[0]);
    cSourceCL = oclLoadProgSource(cPathAndName, "", &szKernelLength);
    
    cl_program kernel2Program = clCreateProgramWithSource(GPUContext, 1, (const char**) &cSourceCL, &szKernelLength, &err);
    
    if(err != CL_SUCCESS)
	printf("ERROR COMPILING KERNEL 2\n");

    // Build the program (OpenCL JIT compilation)
    err = clBuildProgram(kernel2Program, 0, NULL, NULL, NULL, NULL);
    if(err != CL_SUCCESS)
	printf("Error building kernel 2.\n");
	*/
    cl_kernel kernel2 = clCreateKernel(kernel1Program, "kernel2", &err);
    if(err != CL_SUCCESS)
	printf("Error creating kernel 2.\n");

    //Set Arguments for Kernel1 and 2
    clSetKernelArg(kernel1, 0, sizeof(cl_mem), (void*)&d_graph_nodes);
    clSetKernelArg(kernel1, 1, sizeof(cl_mem), (void*)&d_graph_edges);
    clSetKernelArg(kernel1, 2, sizeof(cl_mem), (void*)&d_graph_mask);
    clSetKernelArg(kernel1, 3, sizeof(cl_mem), (void*)&d_updating_graph_mask);
    clSetKernelArg(kernel1, 4, sizeof(cl_mem), (void*)&d_graph_visited);
    clSetKernelArg(kernel1, 5, sizeof(cl_mem), (void*)&d_cost);
    clSetKernelArg(kernel1, 6, sizeof(unsigned int), (void*)&no_of_nodes);
   
    clSetKernelArg(kernel2, 0, sizeof(cl_mem), (void*)&d_graph_mask);
    clSetKernelArg(kernel2, 1, sizeof(cl_mem), (void*)&d_updating_graph_mask);
    clSetKernelArg(kernel2, 2, sizeof(cl_mem), (void*)&d_graph_visited);
    clSetKernelArg(kernel2, 3, sizeof(cl_mem), (void*)&d_over);
    clSetKernelArg(kernel2, 4, sizeof(unsigned int), (void*)&no_of_nodes);

    int k = 0;
    int stop;
    size_t WorkSize[1] = {no_of_nodes}; // one dimensional Range
    size_t localWorkSize[1] = {64};
    cl_event syncEvent;
    do
    {
	stop = 0;
	//Copy stop to device
	clEnqueueWriteBuffer(commandQueue, d_over, CL_TRUE, 0, sizeof(int), (void*)&stop, 0, NULL, NULL);

	//Run Kernel1 and Kernel2
	cl_int err = clEnqueueNDRangeKernel(commandQueue, kernel1, 1, NULL,
		WorkSize, localWorkSize, 0, NULL, NULL);

	if(err != CL_SUCCESS)
	    printf("Error occurred running kernel1.(%d)\n", err);

	err = clEnqueueNDRangeKernel(commandQueue, kernel2, 1, NULL,
		WorkSize, localWorkSize, 0, NULL, NULL);
	
	if(err != CL_SUCCESS)
	    printf("Error occurred running kernel2.\n");
	
	//Copy stop from device
	clEnqueueReadBuffer(commandQueue, d_over, CL_TRUE, 0, sizeof(int), (void*)&stop, 0, NULL, NULL);
	
	k++;
    }while(stop == 1);

    printf("Kernel Executed %d times\n", k);

    //copy result form device to host
    clEnqueueReadBuffer(commandQueue, d_cost, CL_TRUE, 0, sizeof(int)*no_of_nodes, (void*)h_cost, 0, NULL, NULL);

    //Store the result into a file
    FILE* fpo = fopen("result.txt", "w");
    for(unsigned int i = 0; i < no_of_nodes; i++)
    	fprintf(fpo, "%d) cost:%d\n", i, h_cost[i]);
    fclose(fpo);
    printf("Result stored in result.txt\n");

    //cleanup memory
    //Free Host memory
    free(h_graph_nodes);
    free(h_graph_edges);
    free(h_graph_mask);
    free(h_updating_graph_mask);
    free(h_graph_visited);
    free(h_cost);
    //Free memory memory
    //free(GPUDevices);
    clReleaseKernel(kernel1);
    clReleaseKernel(kernel2);
    clReleaseProgram(kernel1Program);
    //clReleaseProgram(kernel2Program);
    clReleaseCommandQueue(commandQueue);
    clReleaseContext(GPUContext);
    clReleaseMemObject(d_graph_nodes);
    clReleaseMemObject(d_graph_edges);
    clReleaseMemObject(d_graph_mask);
    clReleaseMemObject(d_updating_graph_mask);
    clReleaseMemObject(d_graph_visited);
    clReleaseMemObject(d_cost);
    clReleaseMemObject(d_over);
}

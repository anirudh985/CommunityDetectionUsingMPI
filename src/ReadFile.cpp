/*
* @Author : AADITYA CHAUHAN
*
* The code is tailored to work for west1505 Dataset
* The dataset contains 3010 Vertices and 10890 Edges,
* Original Dataset has the first 2 lines removes, to ease offset calculation.
* c R-MAT graph saved from XMT	p sp 3010 10890. Total Datasize = 85364048
*/

#include "ReadFile.h"

#define NUM_THREADS 8

using namespace std;

void buildGraphFromFile(Graph* G){

	FILE *filePointer;
	long fileSize;
	char *buffer;
	size_t result;

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Run code formattor
//	 runFormattor();

	filePointer = fopen("inputToLouvain_copapers", "rb");
	if (filePointer == NULL) {
		fputs("File error", stderr);
		exit(1);
	}

	//get file size
	fseek(filePointer, 0, SEEK_END);
	fileSize = ftell(filePointer);
	rewind(filePointer);

	//Allocate memory for allocating file
	buffer = (char *)malloc(sizeof(char) * fileSize);
	if (buffer == NULL) {
		cout << "Unable to Allocate Memory" << "\n";
	}

	// copy the file into the buffer:
	result = fread(buffer, 1, fileSize, filePointer);
	if (result != fileSize) {
		cout << "Unable to Read from File" << "\n";
		free(buffer);
		fclose(filePointer);
		exit(1);
	}

	readFileInParallel(buffer, fileSize, getLineSize(), size, rank, G);

	fclose(filePointer);
	free(buffer);
}

void readFileInParallel(char* buffer, size_t size, int lineSize, int numberOfProcess, int rank, Graph* G) {
	unsigned long numOfEdges = (size / lineSize);

	// todo : set graph number of edges to the number of lines here
	// number of Edges will be the total number of lines / number of processes

	unsigned long* numOfEdgesOnEachProc = (unsigned long*)malloc(numberOfProcess * sizeof(unsigned long));
	unsigned long* numOfEdgesOnEachProcPrefixSum = (unsigned long*)malloc((numberOfProcess + 1) * sizeof(unsigned long));

	calculateRanges(numOfEdgesOnEachProc, numOfEdgesOnEachProcPrefixSum, rank, numberOfProcess, numOfEdges);

	char *x = (char *)malloc(MAX_WORD_SIZE);
	char *y = (char *)malloc(MAX_WORD_SIZE);

	strncpy(x, buffer + numOfEdgesOnEachProcPrefixSum[rank] * lineSize, MAX_WORD_SIZE);
	strncpy(y, buffer + (numOfEdgesOnEachProcPrefixSum[rank + 1] - 1) * lineSize, MAX_WORD_SIZE);
	unsigned long numOfVerticesOnEachProc = strtoul(y, NULL, 10) - strtoul(x, NULL, 10) + 1;
//	unsigned long numOfVertices = 3010;

	strncpy(x, buffer, MAX_WORD_SIZE);
	strncpy(y, buffer + (numOfEdges - 1) * lineSize, MAX_WORD_SIZE);
	unsigned long numOfVertices = strtoul(y, NULL, 10) - strtoul(x, NULL, 10) + 1;

//	cout << "Number of vertices for rank " << rank << " is " << numOfVertices << "\n";

	unsigned long* vertexStartPointers = (unsigned long*)malloc(sizeof(unsigned long) * numOfVerticesOnEachProc + 1);
	assert(vertexStartPointers != 0);
	unsigned long* startVertices = (unsigned long*)malloc(sizeof(unsigned long) * numOfEdgesOnEachProc[rank]);
	assert(startVertices != 0);
	unsigned long* destinationVertices = (unsigned long*)malloc(sizeof(unsigned long) * numOfEdgesOnEachProc[rank]);
	assert(destinationVertices != 0);
	long* weights = (long*)malloc(sizeof(long) * numOfEdgesOnEachProc[rank]);
	assert(weights != 0);

	unsigned long prevLong = -1l;
	int vCounter = 0;
	int eCounter = 0;

	char *s1 = (char *)malloc(MAX_WORD_SIZE);
	char *s2 = (char *)malloc(MAX_WORD_SIZE);
	char *s3 = (char *)malloc(MAX_WORD_SIZE);
	//	#pragma omp parallel for
	//	for (int i = omp_get_thread_num(); i < numberOfLines; i += omp_get_num_threads()) {
	unsigned long j = 0;
	for (unsigned long i = numOfEdgesOnEachProcPrefixSum[rank]; i < numOfEdgesOnEachProcPrefixSum[rank + 1]; i++, j++) {
		strncpy(s1, buffer + i * lineSize, MAX_WORD_SIZE);
		strncpy(s2, buffer + i * lineSize + MAX_WORD_SIZE + 1, MAX_WORD_SIZE);
		strncpy(s3, buffer + i * lineSize + 2 * MAX_WORD_SIZE + 2, MAX_WORD_SIZE);

		unsigned long source = strtoul(s1, NULL, 10);
		unsigned long destination = strtoul(s2, NULL, 10);
		long weight = strtol(s3, NULL, 10);
		assert(weight == 1);

		if (prevLong != source) {
			vertexStartPointers[vCounter++] = j;
		}

		startVertices[j] = source;
		destinationVertices[j] = destination;
//		weights[j] = (rand() % 1000) + 1;
		weights[j] = weight;

		prevLong = source;
	}
	vertexStartPointers[vCounter] = j;

	G->numOfVertices = numOfVertices;
	G->numOfEdges = numOfEdges;
	G->vertexStartPointers = vertexStartPointers;
	G->startVertices = startVertices;
	G->destinationVertices = destinationVertices;
	G->weights = weights;

	free(s1);
	free(s2);
	free(s3);

}

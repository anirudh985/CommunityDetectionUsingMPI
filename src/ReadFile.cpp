/*
* @Author : AADITYA CHAUHAN
*
* The code is tailored to work for west1505 Dataset
* The dataset contains 3010 Vertices and 10890 Edges,
* Original Dataset has the first 2 lines removes, to ease offset calculation.
* c R-MAT graph saved from XMT	p sp 3010 10890. Total Datasize = 85364048
*/

#include <sstream>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <omp.h>
#include "ReformatFile.h"
#include <mpi.h>

#define NUM_THREADS 8

using namespace std;
void readFileInParallel(char* buffer, size_t size, int lineSize, int numberOfProcess, int processID);

int main(int argc, char *argv[]) {

	FILE *filePointer;
	long fileSize;
	char *buffer;
	size_t result;

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Run code formattor
	// runFormattor();

	filePointer = fopen("out.txt", "rb");
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
	}

	readFileInParallel(buffer, fileSize, getLineSize(), size, rank);

	fclose(filePointer);
	free(buffer);
	return 0;
}

void readFileInParallel(char* buffer, size_t size, int lineSize, int numberOfProcess, int rank) {
	int numberOfLines = ((size / lineSize) / numberOfProcess);

	// todo : set graph number of edges to the number of lines here
	// number of Edges will be the total number of lines / number of processes

	int lineNum1 = ((rank + 1) * numberOfLines) - 1;
	int lineNum2 = ((rank) * numberOfLines) - 1;

	unsigned long numOfVertices = strtoul(buffer + lineNum1 * lineSize, NULL, 10) - strtoul(buffer + lineNum2 * lineSize, NULL, 10);

	cout << "Number of vertices for rank " << rank << " is " << numOfVertices << "\n";

	unsigned long numOfEdges = numberOfLines;
	unsigned long* vertexStartPointers = (unsigned long*)malloc(sizeof(unsigned long) * numOfVertices + 1);
	unsigned long* startVertices = (unsigned long*)malloc(sizeof(unsigned long) * numOfEdges);
	unsigned long* destinationVertices = (unsigned long*)malloc(sizeof(unsigned long) * numOfEdges);
	double* weights = (double*)malloc(sizeof(double) * numOfEdges);

	unsigned long prevLong = -1l;
	int vCounter = 0;
	int eCounter = 0;

	//	#pragma omp parallel for
	//	for (int i = omp_get_thread_num(); i < numberOfLines; i += omp_get_num_threads()) {
	for (int i = 0; i < numberOfLines; i++) {
		char *s1 = (char *)malloc(MAX_WORD_SIZE);
		char *s2 = (char *)malloc(MAX_WORD_SIZE);
		char *s3 = (char *)malloc(2 * MAX_WORD_SIZE);

		strncpy(s1, buffer + i * lineSize, MAX_WORD_SIZE);
		strncpy(s2, buffer + i * lineSize + MAX_WORD_SIZE + 1, MAX_WORD_SIZE);
		strncpy(s3, buffer + i * lineSize + 2 * MAX_WORD_SIZE + 2, 2 * MAX_WORD_SIZE);

		unsigned long source = strtoul(s1, NULL, 10);
		unsigned long destination = strtoul(s2, NULL, 10);
		double weight = strtod(s3, NULL);

		if (prevLong != source) {
			vertexStartPointers[++vCounter] = i;
		}

		startVertices[i] = source;
		destinationVertices[i] = destination;
		weights[i] = weight;

		prevLong = source;
	}
	getchar();
}
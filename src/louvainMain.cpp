/*
 * louvainMain.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: osu8229
 */

#include "louvainMain.h"

using namespace std;

int main(int argc, char* argv[]){

	MPI::Init(argc, argv);
//	MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);

	int rank = MPI::COMM_WORLD.Get_rank();

	int numberOfThreads;

	#pragma omp parallel
	{
		numberOfThreads = omp_get_num_threads();
	}

	if(rank == 0)
		cout << "Number of Threads " << numberOfThreads << endl;

	Graph *G = new Graph();

// TODO: Write a function to create Graph from input File
	buildGraphFromFile(G);
//	buildGraph(G, rank);
	unsigned long numOfVertices = G->numOfVertices;
	unsigned long *finalCommunities = (unsigned long *) malloc(numOfVertices * sizeof(unsigned long));

	runLouvain(G, finalCommunities, 1E-12);


	if(rank == 0){
		char outFile[256];
		sprintf(outFile,"%s_clustInfo", "inputGraph");
		printf("Cluster information will be stored in file: %s\n", outFile);
		FILE* out = fopen(outFile,"w");
		for(long i = 0; i < numOfVertices; i++) {
		  fprintf(out,"%ld\n", finalCommunities[i]);
		}
		fclose(out);
	}

	MPI::Finalize();

	return 0;
}

void buildGraph(Graph *G, int rank){
	G->numOfVertices = 8;
	G->numOfEdges = 19;

	if(rank == 0){
		G->vertexStartPointers = (unsigned long *)malloc(2 * sizeof(unsigned long));
		G->startVertices = (unsigned long *)malloc(5 * sizeof(unsigned long));
		G->destinationVertices = (unsigned long *)malloc(5 * sizeof(unsigned long));
		G->weights = (long *)malloc(5 * sizeof(long));

		G->startVertices[0] = 0;	G->destinationVertices[0] = 4;	G->weights[0] = 10;
		G->startVertices[1] = 0;	G->destinationVertices[1] = 2;	G->weights[1] = 7;
		G->startVertices[2] = 0;	G->destinationVertices[2] = 5;	G->weights[2] = 6;
		G->startVertices[3] = 1;	G->destinationVertices[3] = 3;	G->weights[3] = 5;
		G->startVertices[4] = 1;	G->destinationVertices[4] = 4;	G->weights[4] = 9;

		G->vertexStartPointers[0] = 0;
		G->vertexStartPointers[1] = 3;
	}
	else if(rank == 1){
		G->vertexStartPointers = (unsigned long *)malloc(2 * sizeof(unsigned long));
		G->startVertices = (unsigned long *)malloc(5 * sizeof(unsigned long));
		G->destinationVertices = (unsigned long *)malloc(5 * sizeof(unsigned long));
		G->weights = (long *)malloc(5 * sizeof(long));

		G->startVertices[0] = 1;	G->destinationVertices[0] = 2;	G->weights[0] = 2;
		G->startVertices[1] = 2;	G->destinationVertices[1] = 0;	G->weights[1] = 7;
		G->startVertices[2] = 2;	G->destinationVertices[2] = 1;	G->weights[2] = 2;
		G->startVertices[3] = 2;	G->destinationVertices[3] = 5;	G->weights[3] = 7;
		G->startVertices[4] = 2;	G->destinationVertices[4] = 2;	G->weights[4] = 15;

		G->vertexStartPointers[0] = 0;
		G->vertexStartPointers[1] = 1;
	}
	else if(rank == 2){
		G->vertexStartPointers = (unsigned long *)malloc(3 * sizeof(unsigned long));
		G->startVertices = (unsigned long *)malloc(5 * sizeof(unsigned long));
		G->destinationVertices = (unsigned long *)malloc(5 * sizeof(unsigned long));
		G->weights = (long *)malloc(5 * sizeof(long));

		G->startVertices[0] = 3;	G->destinationVertices[0] = 1;	G->weights[0] = 5;
		G->startVertices[1] = 3;	G->destinationVertices[1] = 5;	G->weights[1] = 6;
		G->startVertices[2] = 4;	G->destinationVertices[2] = 0;	G->weights[2] = 10;
		G->startVertices[3] = 4;	G->destinationVertices[3] = 1;	G->weights[3] = 9;
		G->startVertices[4] = 5;	G->destinationVertices[4] = 0;	G->weights[4] = 6;

		G->vertexStartPointers[0] = 0;
		G->vertexStartPointers[1] = 2;
		G->vertexStartPointers[2] = 4;
	}
	else if(rank == 3){
		G->vertexStartPointers = (unsigned long *)malloc(3 * sizeof(unsigned long));
		G->startVertices = (unsigned long *)malloc(4 * sizeof(unsigned long));
		G->destinationVertices = (unsigned long *)malloc(4 * sizeof(unsigned long));
		G->weights = (long *)malloc(4 * sizeof(long));

		G->startVertices[0] = 5;	G->destinationVertices[0] = 3;	G->weights[0] = 6;
		G->startVertices[1] = 5;	G->destinationVertices[1] = 2;	G->weights[1] = 7;
		G->startVertices[2] = 6;	G->destinationVertices[2] = 7;	G->weights[2] = 4;
		G->startVertices[3] = 7;	G->destinationVertices[3] = 6;	G->weights[3] = 4;

		G->vertexStartPointers[0] = 0;
		G->vertexStartPointers[1] = 2;
		G->vertexStartPointers[2] = 3;
	}







}

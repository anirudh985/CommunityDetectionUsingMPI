/*
 * louvainMultiPhase.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: osu8229
 */

#include "louvainMultiPhase.h"

using namespace std;

void runLouvain(Graph *G, unsigned long *finalCommunities, double threshold){

	int numProcs = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

	double totTimeClustering = 0, totTimeBuildingPhase = 0, tmpTime = 0;
	int tmpItr = 0, totItr = 0;

	unsigned long numOfVertices = G->numOfVertices;
	unsigned long numOfEdges = G->numOfEdges;

	double prevModularity = -1;
	double currModularity = -1;
	unsigned long numClusters = 0;

	Graph *newGraph;


	unsigned long *C = (unsigned long *) malloc (numOfVertices * sizeof(unsigned long));
	assert(C != 0);

	//TODO: Need to see if this is necessary
	#pragma omp parallel for
	for (unsigned long i = 0; i < numOfVertices; i++) {
		C[i] = -1;
	}

	int phaseNumber = 1;

	while(true){
		prevModularity = currModularity;

		currModularity = louvain(G, C, currModularity, threshold, &tmpTime, &tmpItr, rank, numProcs);
		totTimeClustering += tmpTime;
		totItr += tmpItr;

//		if(rank == 0)
//			cout << "Phase: " << phaseNumber << " --> " << "Modularity: " << currModularity << endl;

		numClusters = renumberClustersContiguously(C, G->numOfVertices);
//		printf("Number of unique clusters: %ld\n", numClusters);

		if(phaseNumber == 1){
			#pragma omp parallel for
			for (unsigned long i = 0; i < numOfVertices; i++) {
				finalCommunities[i] = C[i];
			}
		}
		else{
			#pragma omp parallel for
			for (unsigned long i = 0; i < numOfVertices; i++) {
				assert(finalCommunities[i] < G->numOfVertices);
				if(finalCommunities[i] >= 0)
					finalCommunities[i] = C[finalCommunities[i]];
			}
		}

		if((phaseNumber > 10) || (totItr > 10000)){
			break;
		}

		if(currModularity > prevModularity){
			newGraph = new Graph();
			tmpTime = consolidateGraphForNextPhase(G, newGraph, C, numClusters);
			totTimeBuildingPhase += tmpTime;

			delete G;

			G = newGraph;
			G->vertexStartPointers = newGraph->vertexStartPointers;
			G->startVertices = newGraph->startVertices;
			G->destinationVertices = newGraph->destinationVertices;
			G->weights = newGraph->weights;
			G->numOfVertices = newGraph->numOfVertices;
			G->numOfEdges = newGraph->numOfEdges;

			//Free up the previous cluster & create new one of a different size
			free(C);
			C = (unsigned long *) malloc (numClusters * sizeof(unsigned long)); assert(C != 0);

			#pragma omp parallel for
			for (unsigned long i = 0; i < numClusters; i++) {
				C[i] = -1;
			}
			phaseNumber++;

			}
		else{
			break;
		}

	}


	if(rank == 0){
		printf("********************************************\n");
		printf("*********    Compact Summary   *************\n");
		printf("********************************************\n");
		printf("Total number of phases         : %d\n", phaseNumber);
		printf("Total number of iterations     : %d\n", totItr);
		printf("Final number of clusters       : %ld\n", numClusters);
		printf("Final modularity               : %lf\n", prevModularity);
		printf("Total time for clustering      : %lf\n", totTimeClustering);
		printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
		printf("********************************************\n");
		printf("TOTAL TIME                     : %lf\n", (totTimeClustering + totTimeBuildingPhase));
		printf("********************************************\n");
	}

	free(C);
	delete G;
}

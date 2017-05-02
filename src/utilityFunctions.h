/*
 * utilityFunctions.h
 *
 *  Created on: Apr 16, 2017
 *      Author: osu8229
 */

#ifndef UTILITYFUNCTIONS_H_
#define UTILITYFUNCTIONS_H_

#include <iostream>
#include <assert.h>
#include <mpi.h>
#include <omp.h>
#include <map>
#include "Graph.h"
#include <cstdlib>

inline unsigned long getSizeOfArray(unsigned long* arr){
	return (unsigned long) sizeof(arr)/sizeof(unsigned long);
}

inline int ceil(double number){
	return (number - (int)number) > 0 ? (int)number + 1 : (int)number;
}

void calculateRanges(unsigned long* valuesOnEachProc, unsigned long* valuesOnEachProcPrefixSum, int rank, int numProcs, unsigned long valueToDivideAmongProcs);

int lookup(unsigned long* arr, unsigned long target, int sizeOfArr);

unsigned long vertexFollowing(Graph *G, unsigned long *communities);

unsigned long getNumOfVerticesOnProc(unsigned long numOfEdges, unsigned long *startVertices);

unsigned long getNumOfEdgesOnProc(unsigned long numOfEdges, int rank, int numProcs);

unsigned long renumberClustersContiguously(unsigned long *C, unsigned long size);

double consolidateGraphForNextPhase(Graph *G, Graph *newG, unsigned long *C, unsigned long numClusters);

void calculateStartAndEndIndices(unsigned long* startIndices, unsigned long* endIndices, int numProcs,
						unsigned long valueToDivideAmongProcs);

#endif /* UTILITYFUNCTIONS_H_ */

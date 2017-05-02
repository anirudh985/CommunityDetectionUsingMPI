/*
 * parallelLouvain.h
 *
 *  Created on: Apr 26, 2017
 *      Author: osu8229
 */

#ifndef PARALLELLOUVAIN_H_
#define PARALLELLOUVAIN_H_

#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <map>
#include "utilityFunctions.h"
#include "Graph.h"

using namespace std;

double louvain(Graph *G, unsigned long *communityInfo, double Lower,
				double threshold, double *totTime, int *numItr,
				int rank, int numProcs);

void sumVertexDegrees(unsigned long *vertexStartPointers,
					  unsigned long *startVertices,
					  long *weights,
					  long *vertexDegrees,
					  unsigned long numOfVertices,
					  unsigned long numOfVerticesOnProc,
					  unsigned long numOfEdgesOnThisProc,
					  unsigned long *sizeOfCommunities,
					  long *degreesOfCommunities,
					  int *offsetToSend);

double calculateConstantForSecondTerm(unsigned long numOfVertices,
									  long *degreesOfCommunities);

void initialCommunityAssignment(unsigned long *pastCommunityAssignment,
								unsigned long *currCommunityAssignment,
								unsigned long numOfVertices);

void sendRecvUpdateFromPrevProc(unsigned long *startVertices,
								long *vertexDegrees,
								int rank,
								int *offsetToSend);

void sendRecvUpdateFromNextProc(unsigned long *startVertices,
								long *vertexDegrees,
								unsigned long numOfVerticesOnProc,
								unsigned long numOfEdgesOnThisProc,
								int rank);

double buildCommunityDegreeMap(unsigned long startPointer, unsigned long endPointer,
		   	   	   	   	   	   map<unsigned long, long> &communityDegreeMap,
		   	   	   	   	   	   unsigned long *destinationVertices,
		   	   	   	   	   	   long *weights,
							   unsigned long *currCommunityAssignment, unsigned long currVertex);

unsigned long findTargetCommunityOfCurrVertex(map<unsigned long, long> &communityDegreeMap,
								long selfLoop,
								unsigned long *sizeOfCommunities,
								long *degreesOfCommunities,
								long vertexDegree,
								unsigned long currVertex,
								double constantForSecondTerm,
								unsigned long *currCommunityAssignment);

double calculateModularity(double *eii, double constantForSecondTerm,
						   long *degreesOfCommunities, unsigned long numOfVertices,
						   unsigned long numOfVerticesOnProc);

#endif /* PARALLELLOUVAIN_H_ */

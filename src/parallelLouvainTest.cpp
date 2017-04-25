/*
 * parallelLouvainTest.cpp
 *
 *  Created on: Apr 14, 2017
 *      Author: osu8229
 */

#include <mpi.h>
#include <omp.h>
#include <map>
#include "utilityFunctions.h"

using namespace std;


// Assumes MPI::Init is already called in main
double louvain(Graph *G, unsigned long *communityInfo, double Lower,
				double threshold, double *totTime, int *numItr){

	int numberOfThreads = omp_get_num_threads();

	unsigned long numOfVertices = G->numOfVertices;
	unsigned long numOfEdges = G->numOfEdges;
	unsigned long* vertexStartPointers = G->vertexStartPointers; // size of this is numOfVerticesOnProc + 1
	unsigned long* startVertices = G->startVertices;			   // size of this is numOfEdgesOnProc
	unsigned long* destinationVertices = G->destinationVertices; // size of this is numOfEdgesOnProc
	double* weights = G->weights;


	unsigned long numOfVerticesOnProc = getNumOfVerticesOnProc(numOfEdges, startVertices);
	double* vertexDegrees = (double *) malloc(numOfVerticesOnProc * sizeof(double));
	unsigned long totalEdgeWeightTwice;
	double constantForSecondTerm;
	double prevModularity = -1;
	double currModularity = -1;
	double thresMod = threshold; //Input parameter
	int numItrs = 0;

	unsigned long* sizeOfCommunities = (unsigned long *) malloc(numOfVertices * sizeof(unsigned long));
	double* degreesOfCommunities = (double *) malloc(numOfVertices * sizeof(double));

	long* updateSizeOfCommunities = (long *) malloc(numOfVertices * sizeof(long));
	double* updateDegreesOfCommunities = (double *) malloc(numOfVertices * sizeof(double));

	double* eii = (double *) malloc(numOfVerticesOnProc * sizeof(double));

	int offsetToSend = 0;

	sumVertexDegrees(vertexStartPointers,
					 weights,
					 vertexDegrees,
					 numOfVertices,
					 numOfVerticesOnProc,
					 sizeOfCommunities,
					 degreesOfCommunities,
					 &offsetToSend);

	constantForSecondTerm = calculateConstantForSecondTerm();


	unsigned long* pastCommunityAssignment = (unsigned long *) malloc(numOfVertices * sizeof(unsigned long));
	unsigned long* currCommunityAssignment = (unsigned long *) malloc(numOfVertices * sizeof(unsigned long));
	unsigned long* targetCommunityAssignment = (unsigned long *) malloc(numOfVerticesOnProc * sizeof(unsigned long));

	initialCommunityAssignment(pastCommunityAssignment,
							   currCommunityAssignment,
							   numOfVertices);



	// START MAXIMIZING MODULARITY
	while(true){
		numItrs++;

		#pragma omp parallel for
		for(unsigned long i = 0; i < numOfVertices; i++){
			eii[i] = 0;
			updateSizeOfCommunities[i] = 0;
			updateDegreesOfCommunities[i] = 0;
		}

		// FIND THE TARGET COMMUNITIES REGARDLESS OF WHETHER ALL DEST ARE ON THIS PROC
		#pragma omp parallel for
		for(unsigned long i = 0; i < numOfVerticesOnProc; i++){
			unsigned long currVertex = startVertices[vertexStartPointers[i]];

			unsigned long startPointer = vertexStartPointers[i];
			unsigned long endPointer = vertexStartPointers[i+1];

			double selfLoop = 0;

			map<unsigned long, double> communityDegreeMap;
			map<unsigned long, double>::iterator commDegreeMapIterator;

			if(startPointer != endPointer){
				communityDegreeMap[currCommunityAssignment[currVertex]] = 0;

				selfLoop = buildCommunityDegreeMap(startPointer, endPointer,
												   communityDegreeMap,
												   destinationVertices,
												   weights,
												   currCommunityAssignment, currVertex);

				eii[i] = communityDegreeMap[currCommunityAssignment[currVertex]];

				targetCommunityAssignment[i] = findTargetCommunityOfCurrVertex(communityDegreeMap,
																			selfLoop,
																			sizeOfCommunities,
																			degreesOfCommunities,
																			vertexDegrees[i],
																			currVertex,
																			constantForSecondTerm);
			}
			else{
				targetCommunityAssignment[i] = -1;
			}

	         //Update
	        if(targetCommunityAssignment[i] != currCommunityAssignment[currVertex]  && targetCommunityAssignment[i] != -1) {
	        	__sync_fetch_and_add(&updateDegreesOfCommunities[targetCommunityAssignment[i]], vertexDegrees[i]);
	        	__sync_fetch_and_add(&updateSizeOfCommunities[targetCommunityAssignment[i]], 1);
	        	__sync_fetch_and_sub(&updateDegreesOfCommunities[currCommunityAssignment[currVertex]], vertexDegrees[i]);
	        	__sync_fetch_and_sub(&updateSizeOfCommunities[currCommunityAssignment[currVertex]], 1);
	        }

	        communityDegreeMap.clear();
		}// End of for

		// CALCULATE modularity
		currModularity = calculateModularity(eii, constantForSecondTerm,
											 degreesOfCommunities, numOfVertices);

		if(currModularity - prevModularity < thresMod){
			break;
		}

		prevModularity = currModularity;
		if(prevMod < Lower)
			prevMod = Lower;

		MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, updateDegreesOfCommunities,
					  numOfVertices, MPI::DOUBLE, MPI::SUM); //TODO
		MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, updateSizeOfCommunities,
							  numOfVertices, MPI::LONG, MPI::SUM);

		#pragma omp parallel for
		for(unsigned long i = 0; i < numOfVertices; i++){
			degreesOfCommunities[i] += updateDegreesOfCommunities[i];
			sizeOfCommunities[i] += updateSizeOfCommunities[i];
		}

		unsigned long* tmp;
		tmp = pastCommunityAssignment;
		pastCommunityAssignment = currCommunityAssignment; //Previous holds the current
		currCommunityAssignment = tmp; //Current holds the chosen assignment


		int numOfElementsToSend = numOfVerticesOnProc - offsetToSend;
		int recvCounts[size];

		MPI::COMM_WORLD.Alltoall(&numOfElementsToSend, 1, MPI::INT, &recvCounts, 1, MPI::INT);

		//Calculate Prefix Sum
		int displacements[size];
		for(int i = 0; i < size; i++){
			if(i == 0) displacements[0] = 0;
			else{
				displacements[i] = displacements[i-1] + recvCounts[i-1];
			}
		}

		MPI_Allgatherv(&targetCommunityAssignment[offsetToSend],
					   numOfElementsToSend,
					   MPI::UNSIGNED_LONG,
					   currCommunityAssignment,
					   recvCounts,
					   displacements,
					   MPI::UNSIGNED_LONG);//TODO

	}// end of while(true)
//	*totTime = something;//TODO
	*numItr = numItrs;

	#pragma omp parallel for
	for (long i=0; i<NV; i++) {
		communityInfo[i] = pastCommunityAssignment[i];
	}

	free(pastCommunityAssignment);
	free(currCommunityAssignment);
	free(targetCommunityAssignment);
	free(vertexDegrees);
	free(sizeOfCommunities);
	free(degreesOfCommunities);
	free(updateSizeOfCommunities);
	free(updateDegreesOfCommunities);

	return prevModularity;

}

void sumVertexDegrees(unsigned long &vertexStartPointers,
					  unsigned long &startVertices,
					  double &weights,
					  double &vertexDegrees,
					  unsigned long numOfVertices,
					  unsigned long numOfVerticesOnProc,
					  unsigned long &sizeOfCommunities,
					  double &degreesOfCommunities,
					  int *offsetToSend){

	#pragma omp parallel for
	for(unsigned long i = 0; i < numOfVerticesOnProc; i++){
		unsigned long start = vertexStartPointers[i];
		unsigned long end = vertexStartPointers[i+1];
		double totalWeight = 0;
		for(unsigned long j = start; j < end; j++){
			totalWeight += weights[j];
		}
		vertexDegrees[i] = totalWeight;
	}

	int size = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();


	if(rank == 0){
		// send to and receive from 1
		sendRecvUpdateFromNextProc(startVertices, vertexDegrees, numOfVerticesOnProc, rank);
	}

	if(rank == size - 1){
		// send to and receive from rank - 1
		sendRecvUpdateFromPrevProc(startVertices, vertexDegrees, rank, offsetToSend);
	}

	if(rank > 0 && rank < size - 1){
		sendRecvUpdateFromNextProc(startVertices, vertexDegrees, numOfVerticesOnProc, rank);
		sendRecvUpdateFromPrevProc(startVertices, vertexDegrees, rank, offsetToSend);
	}

	int numOfElementsToSend = numOfVerticesOnProc - (*offsetToSend);
	int recvCounts[size];

	MPI::COMM_WORLD.Alltoall(&numOfElementsToSend, 1, MPI::INT, &recvCounts, 1, MPI::INT);

	//Calculate Prefix Sum
	int displacements[size];
	for(int i = 0; i < size; i++){
		if(i == 0) displacements[0] = 0;
		else{
			displacements[i] = displacements[i-1] + recvCounts[i-1];
		}
	}

	MPI_Allgatherv(&vertexDegrees[offsetToSend],
				   numOfElementsToSend,
				   MPI::DOUBLE,
				   degreesOfCommunities,
				   recvCounts,
				   displacements,
				   MPI::DOUBLE);//TODO

	#pragma omp parallel for
	for(unsigned long i = 0; i < numOfVertices; i++){
		sizeOfCommunities[i] = 1;
	}

}

double calculateConstantForSecondTerm(unsigned long numOfVerticesOnProc,
									  double &vertexDegrees){

	double partialTotalWeight = 0;
	double globalTotalWeight = 0;
	#pragma omp parallel for private(i) reduction(+:partialTotalWeight)
	for(unsigned long i = 0; i < numOfVerticesOnProc; i++)
		partialTotalWeight = partialTotalWeight + vertexDegrees[i];

	MPI::COMM_WORLD.Allreduce(&partialTotalWeight,
							  &globalTotalWeight,
							  1,
							  MPI::DOUBLE,
							  MPI::SUM,
							  MPI::COMM_WORLD);

	return globalTotalWeight;
}

void initialCommunityAssignment(unsigned long &pastCommunityAssignment,
								unsigned long &currCommunityAssignment,
								unsigned long numOfVertices){
	#pragma omp parallel for
	for(unsigned long i = 0; i < numOfVertices; i++){
		pastCommunityAssignment[i] = i;
		currCommunityAssignment[i] = i;
	}
}

void sendRecvUpdateFromPrevProc(unsigned long &startVertices,
								double &vertexDegrees,
								int rank,
								int *offsetToSend){
	unsigned long neighborVertex = 0;

	double received = 0;

	MPI::COMM_WORLD.SendRecv(&startVertices[0],
							 1,
							 MPI::UNSIGNED_LONG,
							 rank - 1,
							 0,
							 &neighborVertex,
							 1,
							 MPI::UNSIGNED_LONG,
							 rank - 1,
							 0);
	if(startVertices[0] == neighborVertex){
		(*offsetToSend)++;	//TODO
		MPI::COMM_WORLD.SendRecv(&vertexDegrees[0],
								 1,
								 MPI::DOUBLE,
								 rank - 1,
								 0,
								 &received,
								 1,
								 MPI::DOUBLE,
								 rank - 1,
								 0);

		vertexDegrees[0] += received;
	}

}

void sendRecvUpdateFromNextProc(unsigned long &startVertices,
								double &vertexDegrees,
								unsigned long numOfVerticesOnProc,
								int rank){
	unsigned long neighborVertex = 0;

	double received = 0;

	MPI::COMM_WORLD.SendRecv(&startVertices[numOfVerticesOnProc - 1],
							 1,
							 MPI::UNSIGNED_LONG,
							 rank + 1,
							 0,
							 &neighborVertex,
							 1,
							 MPI::UNSIGNED_LONG,
							 rank + 1,
							 0);

	if(startVertices[numOfVerticesOnProc - 1] == neighborVertex){
		MPI::COMM_WORLD.SendRecv(&vertexDegrees[numOfVerticesOnProc - 1],
								 1,
								 MPI::DOUBLE,
								 rank + 1,
								 0,
								 &received,
								 1,
								 MPI::DOUBLE,
								 rank + 1,
								 0);

		vertexDegrees[numOfVerticesOnProc - 1] += received;
	}

}

double buildCommunityDegreeMap(unsigned long startPointer, unsigned long endPointer,
		   	   	   	   	   	   map<unsigned long, double> &communityDegreeMap,
		   	   	   	   	   	   unsigned long &destinationVertices,
		   	   	   	   	   	   double &weights,
							   unsigned long &currCommunityAssignment, unsigned long currVertex){
	double selfLoop = 0;
	map<unsigned long, double>::iterator commDegreeMapIterator;
	for(unsigned long i = startPointer; i < endPointer; i++){
		if(destinationVertices[i] == currVertex){
			selfLoop += weights[i];
		}

		commDegreeMapIterator = communityDegreeMap.find(currCommunityAssignment[destinationVertices[i]]);
		if(commDegreeMapIterator != communityDegreeMap.end()){
			commDegreeMapIterator->second += weights[i];
		}
		else{
			communityDegreeMap[currCommunityAssignment[destinationVertices[i]]] = weights[i];
		}
	}

	return selfLoop;
}

unsigned long findTargetCommunityOfCurrVertex(map<unsigned long, double> &communityDegreeMap,
								double selfLoop,
								unsigned long &sizeOfCommunities,
								double &degreesOfCommunities,
								double vertexDegree,
								unsigned long currVertex,
								double constantForSecondTerm,
								unsigned long &currCommunityAssignment){

	map<unsigned long, double>::iterator commDegreeMapIterator;
	long maxIndex = currVertex;	//Assign the initial value as self community
	double curGain = 0;
	double maxGain = 0;
	unsigned long currVertexCommunity = currCommunityAssignment[currVertex];
	double eix = communityDegreeMap[currVertexCommunity]- selfLoop;
	double ax = degreesOfCommunities[currVertex] - vertexDegree;
	double eiy = 0;
	double ay = 0;

	commDegreeMapIterator = communityDegreeMap.begin();

	do{
		if(currVertexCommunity != commDegreeMapIterator->first) {
			ay = degreesOfCommunities[commDegreeMapIterator->first];
			eiy = commDegreeMapIterator->second;
		    curGain = 2*(eiy - eix) - 2*vertexDegree*(ay - ax)*constantForSecondTerm;

		    if( (curGain > maxGain) ||
		          ((curGain == maxGain) && (curGain != 0) && (storedAlready->first < maxIndex)) ) {
		    	maxGain = curGain;
		    	maxIndex = storedAlready->first;
		      }
		    }
		    commDegreeMapIterator++;
	}while(commDegreeMapIterator != communityDegreeMap.end());

	if(sizeOfCommunities[maxIndex] == 1 && sizeOfCommunities[currVertexCommunity] == 1 && maxIndex > currVertexCommunity){
	    maxIndex = currVertexCommunity;
	  }

	return maxIndex;

}

double calculateModularity(double &eii, double constantForSecondTerm,
						   double &degreesOfCommunities, unsigned long numOfVertices,
						   unsigned long numOfVerticesOnProc){

	int size = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

	double partialE2_XX = 0;
	double partialA2_X = 0;

	double e_xx__a2_x[2] = {0, 0};

	#pragma omp parallel for reduction(+:partialE2_XX)
	for(unsigned long i = 0; i < numOfVerticesOnProc; i++){
		partialE2_XX += eii[i];
	}

	unsigned long block = numOfVertices / size;
	#pragma omp parallel for reduction(+:partialA2_X)
	for(unsigned long i = rank * block; i < (rank + 1) * block && i < numOfVertices; i++){
		partialA2_X += degreesOfCommunities[i];
	}

	double partial_e2_a2[2] = {partialE2_XX, partialA2_X};

	MPI::COMM_WORLD.Allreduce(&partial_e2_a2, &e_xx__a2_x, 2, MPI::DOUBLE, MPI::SUM);
	return (e_xx__a2_x[0]*constantForSecondTerm) - (e_xx__a2_x[1]*constantForSecondTerm*constantForSecondTerm);
}

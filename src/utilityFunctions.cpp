/*
 * utilityFunctions.cpp
 *
 *  Created on: Apr 17, 2017
 *      Author: osu8229
 */


#include "utilityFunctions.h"

using namespace std;

unsigned long vertexFollowing(Graph *G, unsigned long * communities){

	unsigned long numOfVertices = G->numOfVertices;
	unsigned long numOfEdges = G->numOfEdges;
	unsigned long* vertexStartPointers = G->vertexStartPointers; // size of this is numOfVerticesOnProc + 1
	unsigned long* startVertices = G->startVertices;			   // size of this is numOfEdgesOnProc
	unsigned long* destinationVertices = G->destinationVertices; // size of this is numOfEdgesOnProc
	long* weights = G->weights;

	unsigned long numVerticesToMerge = 0;

	#pragma omp parallel for
	for(unsigned long i = 0; i < numOfVertices; i++){
		communities[i] = i;
	}

	unsigned long numOfVerticesOnProc = getNumOfVerticesOnProc(numOfEdges, startVertices);






	return 0;
}


unsigned long getNumOfVerticesOnProc(unsigned long numOfEdges,
							unsigned long *startVertices){
	int numProcs = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

	unsigned long numOfEdgesOnProc = getNumOfEdgesOnProc(numOfEdges, rank, numProcs);

	unsigned long startVertexOnProc = startVertices[0];
	unsigned long endVertexOnProc = startVertices[numOfEdgesOnProc - 1];
	return endVertexOnProc - startVertexOnProc + 1;
}

unsigned long getNumOfEdgesOnProc(unsigned long numOfEdges, int rank, int numProcs){
	int remainder = numOfEdges % numProcs;
	unsigned long numOfEdgesOnProc = numOfEdges / numProcs;
	if(rank < remainder) numOfEdgesOnProc++;
	return numOfEdgesOnProc;
}

// TODO: Refactor the below function

//WARNING: Will overwrite the old cluster vector
//Returns the number of unique clusters
unsigned long renumberClustersContiguously(unsigned long *C, unsigned long size) {
#ifdef PRINT_DETAILED_STATS_
  printf("Within renumberClustersContiguously()\n");
#endif
  double time1 = omp_get_wtime();
  //Count the number of unique communities and internal edges
  map<unsigned long, unsigned long> clusterLocalMap; //Map each neighbor's cluster to a local number
  map<unsigned long, unsigned long>::iterator storedAlready;
  unsigned long numUniqueClusters = 0;

  //Do this loop in serial
  //Will overwrite the old cluster id with the new cluster id
  for(unsigned long i=0; i<size; i++) {
    assert(C[i]<size);
    if (C[i] >= 0) { //Only if it is a valid number
      storedAlready = clusterLocalMap.find(C[i]); //Check if it already exists
      if( storedAlready != clusterLocalMap.end() ) {	//Already exists
	C[i] = storedAlready->second; //Renumber the cluster id
      } else {
	clusterLocalMap[C[i]] = numUniqueClusters; //Does not exist, add to the map
	C[i] = numUniqueClusters; //Renumber the cluster id
	numUniqueClusters++; //Increment the number
      }
    }//End of if()
  }//End of for(i)
  time1 = omp_get_wtime() - time1;
//  printf("Time to renumber clusters: %lf\n", time1);

  return numUniqueClusters; //Return the number of unique cluster ids
}//End of renumberClustersContiguously()


void calculateRanges(unsigned long* valuesOnEachProc, unsigned long* valuesOnEachProcPrefixSum, int rank, int numProcs, unsigned long valueToDivideAmongProcs){
	unsigned long quotient = valueToDivideAmongProcs / numProcs;
	int remainder = valueToDivideAmongProcs % numProcs;
	valuesOnEachProcPrefixSum[0] = 0;
//	May be faster without openmp
//	#pragma omp parallel for
	for(int i = 0; i < numProcs; i++){
		if(i < remainder){
			valuesOnEachProc[i] = quotient + 1;
		}
		else{
			valuesOnEachProc[i] = quotient;
		}
		valuesOnEachProcPrefixSum[i + 1] = valuesOnEachProcPrefixSum[i] + valuesOnEachProc[i];
	}
	return;
}

// Currently supports the numOfEdges in the integer range
// Expects both G and newG to be initialized (only G and newG should be pointing to some memory, not the values inside them)
double consolidateGraphForNextPhase(Graph *G, Graph *newG, unsigned long *C, unsigned long numClusters){

	int numProcs = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

	double total = 0, totItr = 0;

	double time1, time2, time;

	time1 = omp_get_wtime();

	unsigned long numOfVertices = G->numOfVertices;
	unsigned long numOfEdges = G->numOfEdges;
	unsigned long* vertexStartPointers = G->vertexStartPointers; // size of this is numOfVerticesOnProc + 1
	unsigned long* startVertices = G->startVertices;			   // size of this is numOfEdgesOnProc
	unsigned long* destinationVertices = G->destinationVertices; // size of this is numOfEdgesOnProc
	long* weights = G->weights;
//	unsigned long numOfEdgesOnProc = getSizeOfArray(G->startVertices);
	// OR
	unsigned long numOfEdgesOnProc = getNumOfEdgesOnProc(numOfEdges, rank, numProcs);



	map<unsigned long,long>** neighborClustersMap = (map<unsigned long,long>**) malloc(numClusters*sizeof(map<unsigned long,long>*));
	assert(neighborClustersMap != 0);

	  //Add for a self loop with zero weight
	#pragma omp parallel for
	for (unsigned long i = 0; i < numClusters; i++) {
		neighborClustersMap[i] = new map<unsigned long,long>();
		(*(neighborClustersMap[i]))[i] = 0;
	}

	  //Create an array of locks for each map corresponding to each cluster
	omp_lock_t *neighClustLocks = (omp_lock_t *) malloc (numClusters * sizeof(omp_lock_t));
	assert(neighClustLocks != 0);

	  //Initialize locks
	#pragma omp parallel for
	for(unsigned long i = 0; i < numClusters; i++) {
		omp_init_lock(&neighClustLocks[i]);
	}

	#pragma omp parallel for
	for(unsigned long i = 0; i < numOfEdgesOnProc; i++){
		unsigned long currVertex = startVertices[i];
		unsigned long currVertexCluster = C[currVertex];
        assert(C[currVertexCluster] < numClusters);
        unsigned long currVertexNeighborCluster = C[destinationVertices[i]];
        assert(C[currVertexNeighborCluster] < numClusters);

        map<unsigned long, long>::iterator it;

        omp_set_lock(&neighClustLocks[currVertexCluster]);  // Locking the currVertexCluster Map

        it = neighborClustersMap[currVertexCluster]->find(currVertexNeighborCluster);
        if(it != neighborClustersMap[currVertexCluster]->end()){
        	it->second += weights[i];
        }
        else{
        	(*(neighborClustersMap[currVertexCluster]))[currVertexNeighborCluster] = weights[i];
        }
        omp_unset_lock(&neighClustLocks[currVertexCluster]); // Unlocking the currVertexCluster Map
	}

	#pragma omp parallel for
	for (unsigned long i = 0; i < numClusters; i++) {
		omp_destroy_lock(&neighClustLocks[i]);
	}
	free(neighClustLocks);

	// size of sizesOfMaps = (numClusters + 1) -- to accomodate for prefix sum calculation
	unsigned long* sizesOfMaps = (unsigned long*)malloc((numClusters + 1) * sizeof(unsigned long));
	unsigned long totalSize = 0;

	sizesOfMaps[0] = 0;
	#pragma omp parallel for reduction(+:totalSize)
	for(unsigned long i = 0; i < numClusters; i++){
		sizesOfMaps[i + 1] = neighborClustersMap[i]->size();
		totalSize += neighborClustersMap[i]->size();
	}

	//TODO: Parallel In-Place Prefix Sum --- currently doing it sequentially
	// Prefix Sum of sizesOfMaps (not changing the variable name to save space)
	for(unsigned long i = 0; i < numClusters; i++){
		sizesOfMaps[i + 1] += sizesOfMaps[i];
	}

	unsigned long* sendStartVerticesWithAllBinRanges = (unsigned long*)malloc(totalSize * sizeof(unsigned long));
	unsigned long* sendDestinationVerticesWithAllBinRanges = (unsigned long*)malloc(totalSize * sizeof(unsigned long));
	long* sendWeightsWithAllBinRanges = (long *)malloc(totalSize * sizeof(long));

	#pragma omp parallel for
	for(unsigned long i = 0; i < numClusters; i++){
		unsigned long startWriteIndex = sizesOfMaps[i];
		map<unsigned long, long>::iterator it;
		unsigned long j = 0;
		for(it = neighborClustersMap[i]->begin(); it != neighborClustersMap[i]->end(); ++it, j++){	// TODO: Proof read
			sendStartVerticesWithAllBinRanges[startWriteIndex + j] = i;
			sendDestinationVerticesWithAllBinRanges[startWriteIndex + j] = it->first;
			sendWeightsWithAllBinRanges[startWriteIndex + j] = it->second;
		}
	}


	// Calculate SendCounts
//	unsigned long binRange = numClusters / numProcs;

	unsigned long* binRanges = (unsigned long*)malloc(numProcs * sizeof(unsigned long));
	unsigned long* binRangePrefixSum = (unsigned long*)malloc((numProcs + 1) * sizeof(unsigned long));

	calculateRanges(binRanges, binRangePrefixSum, rank, numProcs, numClusters);

	int* sendCounts = (int*)malloc(numProcs * sizeof(int));

	// May be faster if used without openmp
//	#pragma omp parallel for
	for(int i = 0; i < numProcs; i++){
		sendCounts[i] = sizesOfMaps[binRangePrefixSum[i+1]] - sizesOfMaps[binRangePrefixSum[i]];
	}

	// Calculate SendDispls
	int* sendDispls = (int*)malloc(numProcs * sizeof(int));

	// May be faster if used without openmp
//	#pragma omp parallel for
	for(int i = 0; i < numProcs; i++){
		sendDispls[i] = sizesOfMaps[binRangePrefixSum[i]];
	}

	// Communicate Counts & Calculate RecvCounts
	int* recvCounts = (int*)malloc(numProcs * sizeof(int));
	assert(recvCounts != 0);

	MPI::COMM_WORLD.Alltoall(sendCounts, 1, MPI::INT,
	    						   recvCounts, 1, MPI::INT);

	// Calculate RecvDispls
	int* recvDispls = (int*)malloc(numProcs * sizeof(int));
	recvDispls[0] = 0;
	for(int i = 0; i < numProcs - 1; i++){
		recvDispls[i + 1] = recvDispls[i] + recvCounts[i];
	}

	unsigned long numEdgesPerProcAfterBinWiseShuffling = recvDispls[numProcs - 1] + recvCounts[numProcs - 1];
	unsigned long* perProcStartVerticesRecv = (unsigned long*)malloc(numEdgesPerProcAfterBinWiseShuffling * sizeof(unsigned long));
	unsigned long* perProcDestinationVerticesRecv = (unsigned long*)malloc(numEdgesPerProcAfterBinWiseShuffling * sizeof(unsigned long));
	long* perProcWeightsRecv = (long *)malloc(numEdgesPerProcAfterBinWiseShuffling * sizeof(long));

	assert(perProcDestinationVerticesRecv != 0);
	assert(perProcStartVerticesRecv != 0);
	assert(perProcWeightsRecv != 0);

	// Communicate Data
	MPI::COMM_WORLD.Alltoallv(sendStartVerticesWithAllBinRanges, sendCounts, sendDispls, MPI::UNSIGNED_LONG,
							   perProcStartVerticesRecv, recvCounts, recvDispls, MPI::UNSIGNED_LONG);

	MPI::COMM_WORLD.Alltoallv(sendDestinationVerticesWithAllBinRanges, sendCounts, sendDispls, MPI::UNSIGNED_LONG,
							   perProcDestinationVerticesRecv, recvCounts, recvDispls, MPI::UNSIGNED_LONG);

	MPI::COMM_WORLD.Alltoallv(sendWeightsWithAllBinRanges, sendCounts, sendDispls, MPI::LONG,
								   perProcWeightsRecv, recvCounts, recvDispls, MPI::LONG);

	// Free unnecessary memory
	#pragma omp parallel for
	for(unsigned long i = 0; i < numClusters; i++){
		free(neighborClustersMap[i]);
	}

	free(neighborClustersMap);
	free(sendStartVerticesWithAllBinRanges);
	free(sendDestinationVerticesWithAllBinRanges);
	free(sendWeightsWithAllBinRanges);

	free(sizesOfMaps);


	// TODO: need to take care of last proc having less than binRange Clusters (I think I did this part)
	map<unsigned long,long>** perProcNeighborClustersMap = (map<unsigned long,long>**) malloc(binRanges[rank] * sizeof(map<unsigned long,long>*));
	assert(perProcNeighborClustersMap != 0);

	#pragma omp parallel for
	for (unsigned long i = 0; i < binRanges[rank]; i++) {
		perProcNeighborClustersMap[i] = new map<unsigned long,long>();
	}

	//Create an array of locks for each map corresponding to each cluster
	omp_lock_t *perProcNeighClustLocks = (omp_lock_t *) malloc (binRanges[rank] * sizeof(omp_lock_t));
	assert(perProcNeighClustLocks != 0);

	  //Initialize locks
	#pragma omp parallel for
	for(unsigned long i = 0; i < binRanges[rank]; i++) {
		omp_init_lock(&perProcNeighClustLocks[i]);
	}

	// OKAY upto this point
	// size of below arrays is fishy

//	unsigned long* perProcStartVerticesIntermediate = (unsigned long*)malloc(numEdgesPerProcAfterBinWiseShuffling * sizeof(unsigned long));
//	unsigned long* perProcDestinationVerticesIntermediate = (unsigned long*)malloc(numEdgesPerProcAfterBinWiseShuffling * sizeof(unsigned long));
//	long* perProcWeightsIntermediate = (long *)malloc(numEdgesPerProcAfterBinWiseShuffling * sizeof(long));


	#pragma omp parallel for
	for(unsigned long i = 0; i < numEdgesPerProcAfterBinWiseShuffling; i++){
		unsigned long offset = binRangePrefixSum[rank];
		unsigned long currVertex = perProcStartVerticesRecv[i];
		unsigned long currVertexNeighbor = perProcDestinationVerticesRecv[i];
		assert(currVertex >= binRangePrefixSum[rank] && currVertex < binRangePrefixSum[rank + 1]);

		map<unsigned long, long>::iterator it;

		omp_set_lock(&perProcNeighClustLocks[currVertex-offset]);  // Locking the currVertex Map
		it = perProcNeighborClustersMap[currVertex-offset]->find(currVertexNeighbor); //TODO Proof read
		if(it != perProcNeighborClustersMap[currVertex-offset]->end()){
			it->second += perProcWeightsRecv[i];
		}
		else{
			(*(perProcNeighborClustersMap[currVertex-offset]))[currVertexNeighbor] = perProcWeightsRecv[i];
		}
		omp_unset_lock(&perProcNeighClustLocks[currVertex-offset]); // Unlocking the currVertex Map
	}

	#pragma omp parallel for
	for (unsigned long i = 0; i < binRanges[rank]; i++) {
		omp_destroy_lock(&perProcNeighClustLocks[i]);
	}
	free(perProcNeighClustLocks);

	sizesOfMaps = (unsigned long*)malloc((binRanges[rank] + 1) * sizeof(unsigned long));
	totalSize = 0;
	sizesOfMaps[0] = 0;

	#pragma omp parallel for reduction(+:totalSize)
	for(unsigned long i = 0; i < binRanges[rank]; i++){
		sizesOfMaps[i + 1] = perProcNeighborClustersMap[i]->size();
		totalSize += perProcNeighborClustersMap[i]->size();
	}

	assert(totalSize <= numEdgesPerProcAfterBinWiseShuffling);

	//TODO: Parallel In-Place Prefix Sum --- currently doing it sequentially
	// Prefix Sum of sizesOfMaps (not changing the variable name to save space)
	for(unsigned long i = 0; i < binRanges[rank]; i++){
		sizesOfMaps[i + 1] += sizesOfMaps[i];
	}

//	free(perProcStartVerticesIntermediate);
//	free(perProcDestinationVerticesIntermediate);
//	free(perProcWeightsIntermediate);

	unsigned long* perProcStartVertices = (unsigned long*)malloc(totalSize * sizeof(unsigned long));
	unsigned long* perProcDestinationVertices = (unsigned long*)malloc(totalSize * sizeof(unsigned long));
	long* perProcWeights = (long *)malloc(totalSize * sizeof(long));

	#pragma omp parallel for
	for(unsigned long i = 0; i < binRanges[rank]; i++){
		unsigned long startWriteIndex = sizesOfMaps[i];
		map<unsigned long, long>::iterator it;
		unsigned long j = 0;
		for(it = perProcNeighborClustersMap[i]->begin(); it != perProcNeighborClustersMap[i]->end(); ++it, j++){
			perProcStartVertices[startWriteIndex + j] = binRangePrefixSum[rank] + i;
			perProcDestinationVertices[startWriteIndex + j] = it->first;
			perProcWeights[startWriteIndex + j] = it->second;
		}
	}

	#pragma omp parallel for
	for(unsigned long i = 0; i < binRanges[rank]; i++){
		free(perProcNeighborClustersMap[i]);
	}

	free(perProcNeighborClustersMap);

	// Now every processor has equal number of vertices, which we need to repartition based on the number of edges
	unsigned long* unpartitionedNumEdgesOnEachProc = (unsigned long*)malloc((numProcs + 1) * sizeof(unsigned long));

	unpartitionedNumEdgesOnEachProc[0] = 0;
	// 1. Communicating #edgesPerProc currently to all procs
	MPI::COMM_WORLD.Allgather(&totalSize, 1, MPI::UNSIGNED_LONG, &unpartitionedNumEdgesOnEachProc[1], 1, MPI::UNSIGNED_LONG);

	// 2. Calculate Prefix sums of edgesOnEachProc
	//TODO: Parallel In-Place Prefix Sum --- currently doing it sequentially
	// Prefix Sum of edgesOnEachProc (not changing the variable name to save space)
	for(int i = 0; i < numProcs; i++){
		unpartitionedNumEdgesOnEachProc[i + 1] += unpartitionedNumEdgesOnEachProc[i];
	}

	free(sendCounts);
	free(sendDispls);
	free(recvCounts);
	free(recvDispls);

	sendCounts = (int*)calloc(numProcs, sizeof(int));
	sendDispls = (int*)calloc(numProcs, sizeof(int));
	recvCounts = (int*)calloc(numProcs, sizeof(int));
	recvDispls = (int*)calloc(numProcs, sizeof(int));

	// 3. Calculate sendCounts for each proc
	unsigned long totalEdges = unpartitionedNumEdgesOnEachProc[numProcs];

	unsigned long* edgesPerProc = (unsigned long*)malloc(numProcs * sizeof(unsigned long));
	unsigned long* edgesPerProcPrefixSum = (unsigned long*)malloc((numProcs + 1) * sizeof(unsigned long));

	calculateRanges(edgesPerProc, edgesPerProcPrefixSum, rank, numProcs, totalEdges);

//	unsigned long edgesPerProc = totalEdges / numProcs;

//	if(totalEdges % numProcs != 0) edgesPerProc++;

	int startProcToSend = 0, lastProcToSend = 0, countToSendToStartProc  = 0, countToSendToLastProc = 0, numProcFromStartToSendPlusOne = 0;

	// ********************** NEW CODE - START ************************ //

	startProcToSend = lookup(edgesPerProcPrefixSum, unpartitionedNumEdgesOnEachProc[rank], numProcs + 1) - 1;

	lastProcToSend = lookup(edgesPerProcPrefixSum, unpartitionedNumEdgesOnEachProc[rank + 1], numProcs + 1) - 1;

	if(startProcToSend == lastProcToSend){
		sendCounts[startProcToSend] = totalSize;
	}
	else{
		countToSendToStartProc = edgesPerProcPrefixSum[startProcToSend + 1] - unpartitionedNumEdgesOnEachProc[rank];
		countToSendToLastProc = totalSize - countToSendToStartProc - (edgesPerProcPrefixSum[lastProcToSend] - edgesPerProcPrefixSum[startProcToSend + 1]);
		sendCounts[startProcToSend] = countToSendToStartProc;
		sendCounts[lastProcToSend] = countToSendToLastProc;
		for(int i = startProcToSend + 1; i < lastProcToSend; i++){
				sendCounts[i] = edgesPerProc[i];
			}
	}
	// ********************** NEW CODE - END   ************************ //

	// ********************** OLD CODE - START   ************************ //

//	startProcToSend = (unpartitionedNumEdgesOnEachProc[rank] / edgesPerProc[rank]);
//	countToSendToStartProc = edgesPerProc - (unpartitionedNumEdgesOnEachProc[rank] % edgesPerProc);
//	if(totalSize - countToSendToStartProc <= 0){
//		lastProcToSend = startProcToSend;
//		countToSendToStartProc = totalSize;
//		countToSendToLastProc = totalSize;
//	}
//	else{
//		numProcFromStartToSendPlusOne = ceil((totalSize - countToSendToStartProc) / (long) edgesPerProc);
//		lastProcToSend = startProcToSend + numProcFromStartToSendPlusOne;
//		countToSendToLastProc = (totalSize - countToSendToStartProc) % edgesPerProc;
//	}
//
//
//	sendCounts[startProcToSend] = countToSendToStartProc;
//
//	// Need not use omp as it might cause overhead
////	#pragma omp parallel for
//	for(int i = startProcToSend + 1; i < lastProcToSend; i++){
//		sendCounts[i] = edgesPerProc;
//	}
//
//	sendCounts[lastProcToSend] = countToSendToLastProc;

	// ********************** OLD CODE - END   ************************ //


	// 4. Calculate sendDispls for each proc => prefix sum of sendCounts
	for(int i = 0; i < numProcs - 1; i++){
		sendDispls[i+1] = sendDispls[i] + sendCounts[i];
	}

	// 5. Calculate recvCounts for each proc

	// ********************** NEW CODE - START ************************ //
	MPI::COMM_WORLD.Alltoall(sendCounts, 1, MPI::INT,
	    						   recvCounts, 1, MPI::INT);
	// ********************** NEW CODE - END  ************************* //

	// ********************** OLD CODE - START  ************************* //
//	unsigned long startEdgeNumberToRecv = edgesPerProcPrefixSum[rank];
//	unsigned long endEdgeNumberToRecv = edgesPerProcPrefixSum[rank + 1] - 1;
//	int startProcToRecv = -1;
//	int endProcToRecv = -1;
//
//	// Need to write Binary Search Routine here
//	// Currently doing sequential search
//	for(int i = 1; i < numProcs; i++){
//		if(unpartitionedNumEdgesOnEachProc[i] > startEdgeNumberToRecv){
//			startProcToRecv = i - 1;
//		}
//		if(startProcToRecv != -1 && unpartitionedNumEdgesOnEachProc[i] > endEdgeNumberToRecv){
//			endProcToRecv = i - 1;
//			break;
//		}
//	}
//
//	int countToRecvFromStartProc = 0;
//	int countToRecvFromLastProc = 0;
//	int numProcsToRecvFromStartProcPlusOne = 0;
//	if(startProcToRecv == endProcToRecv){
//		recvCounts[startProcToRecv] = edgesPerProc[rank]; // TODO: edgesPerProc[startProcToRecv]
//	}
//	else{
//		countToRecvFromStartProc = unpartitionedNumEdgesOnEachProc[startProcToRecv + 1] - startEdgeNumberToRecv;
//		numProcsToRecvFromStartProcPlusOne = endProcToRecv - startProcToRecv;
//		countToRecvFromLastProc = endEdgeNumberToRecv - unpartitionedNumEdgesOnEachProc[endProcToRecv] + 1;
//		recvCounts[startProcToRecv] = countToRecvFromStartProc;
//		recvCounts[endProcToRecv] = countToRecvFromLastProc;
//		for(int i = startProcToRecv + 1; i < endProcToRecv; i++){
//			recvCounts[i] = unpartitionedNumEdgesOnEachProc[i + 1] - unpartitionedNumEdgesOnEachProc[i];
//		}
//	}

	// ********************** OLD CODE - START  ************************* //

	// 6. Calculate recvDispls for each proc
	// Prefix sum to calculate recvDispls
	for(int i = 0; i < numProcs - 1; i++){
		recvDispls[i+1] = recvDispls[i] + recvCounts[i];
	}



	// 7. Allocate memory for newGraph

	unsigned long numEdgesOnThisProc = edgesPerProc[rank];

	newG->startVertices = (unsigned long*)malloc(numEdgesOnThisProc * sizeof(unsigned long));
	newG->destinationVertices = (unsigned long*)malloc(numEdgesOnThisProc * sizeof(unsigned long));
	newG->weights = (long*)malloc(numEdgesOnThisProc * sizeof(long));

	MPI::COMM_WORLD.Alltoallv(perProcStartVertices, sendCounts, sendDispls, MPI::UNSIGNED_LONG,
							   newG->startVertices, recvCounts, recvDispls, MPI::UNSIGNED_LONG);

	MPI::COMM_WORLD.Alltoallv(perProcDestinationVertices, sendCounts, sendDispls, MPI::UNSIGNED_LONG,
							   newG->destinationVertices, recvCounts, recvDispls, MPI::UNSIGNED_LONG);

	MPI::COMM_WORLD.Alltoallv(perProcWeights, sendCounts, sendDispls, MPI::LONG,
							   newG->weights, recvCounts, recvDispls, MPI::LONG);

	newG->numOfVertices = numClusters;
	newG->numOfEdges = totalEdges;

	// Calculate vertexStartPointers
	// TODO: can this be done using openmp?
	unsigned long indexToVertexStartPointers = 1;
	unsigned long numVerticesOnThisProc = newG->startVertices[numEdgesOnThisProc - 1] - newG->startVertices[0] + 1;
	newG->vertexStartPointers = (unsigned long*)malloc(numVerticesOnThisProc * sizeof(unsigned long));
	newG->vertexStartPointers[0] = 0;
	for(unsigned long i = 1; i < numEdgesOnThisProc; i++){
		if(newG->startVertices[i] != newG->startVertices[i - 1]){
			newG->vertexStartPointers[indexToVertexStartPointers] = i;
			indexToVertexStartPointers++;
		}
		if(indexToVertexStartPointers == numVerticesOnThisProc) break;
	}

	free(sendCounts);
	free(sendDispls);
	free(recvCounts);
	free(recvDispls);
	free(perProcStartVertices);
	free(perProcDestinationVertices);
	free(perProcWeights);

	time2 = omp_get_wtime();

	time = time2 - time1;

	MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, &time, 1, MPI::DOUBLE, MPI::MAX);

	return time;
}

int lookup(unsigned long* arr, unsigned long target, int sizeOfArr){
	int lo = 0, hi = sizeOfArr - 1;
	int mid;
	   while(lo < hi){
		   mid = lo + (hi-lo)/2;
		   if(arr[mid] == target){
			   return mid + 1;
		   }
		   else if(arr[mid] < target){
			   lo = mid + 1;
		   }
		   else{
			   hi = mid;
		   }
	   }
	return lo;
}

void calculateStartAndEndIndices(unsigned long* startIndices, unsigned long* endIndices, int numProcs,
						unsigned long valueToDivideAmongProcs){
	unsigned long quotient = valueToDivideAmongProcs / numProcs;
	int remainder = valueToDivideAmongProcs % numProcs;

	unsigned long valuesOnEachProc[numProcs];
	for(int i = 0; i < numProcs; i++){
		if(i < remainder){
			valuesOnEachProc[i] = quotient + 1;
		}else{
			valuesOnEachProc[i] = quotient;
		}
	}
	startIndices[0] = 0;
	endIndices[numProcs - 1] = valueToDivideAmongProcs - 1;
	for(int i = 1; i < numProcs; i++){
		startIndices[i] = startIndices[i - 1] + valuesOnEachProc[i];
		endIndices[i - 1] = startIndices[i] - 1;
	}
}

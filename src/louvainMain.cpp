/*
 * louvainMain.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: osu8229
 */

#include <iostream>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include <Graph.h>

using namespace std;

int main(int argc, char** argv){

	MPI::Init(&argc, &argv);

	Graph *G;

// TODO: Write a function to create Graph from input File
//	buildGraphFromFile(&G);

	unsigned long numOfVerticesToMerge = 0;
	unsigned long *finalCommunities = (unsigned long *) malloc(G->numOfVertices * sizeof(unsigned long));

	numOfVerticesToMerge = vertexFollowing(G, finalCommunities);



	return 0;
}

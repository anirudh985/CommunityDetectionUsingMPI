/*
 * Graph.h
 *
 *  Created on: Apr 9, 2017
 *      Author: osu8229
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <cstdlib>

class Graph {
public:
	Graph();
	virtual ~Graph();

	unsigned long numOfVertices;
	unsigned long numOfEdges;
	unsigned long* vertexStartPointers;
	unsigned long* startVertices;
	unsigned long* destinationVertices;
	long* weights;

};

#endif /* GRAPH_H_ */

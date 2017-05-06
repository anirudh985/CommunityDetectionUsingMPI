/*
 * ReadFile.h
 *
 *  Created on: May 2, 2017
 *      Author: osu8229
 */

#ifndef READFILE_H_
#define READFILE_H_


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
#include "Graph.h"
#include "utilityFunctions.h"

void readFileInParallel(char* buffer, size_t size, int lineSize, int numberOfProcess, int processID, Graph* G);

void buildGraphFromFile(Graph* G);
#endif /* READFILE_H_ */

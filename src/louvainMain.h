/*
 * louvainMain.h
 *
 *  Created on: Apr 26, 2017
 *      Author: osu8229
 */

#ifndef LOUVAINMAIN_H_
#define LOUVAINMAIN_H_

#include <iostream>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include "Graph.h"
#include "louvainMultiPhase.h"
#include "ReadFile.h"

int main(int argc, char** argv);

void buildGraph(Graph *G, int rank);

#endif /* LOUVAINMAIN_H_ */

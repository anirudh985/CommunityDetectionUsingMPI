/*
 * allgatherv_Test.cpp
 *
 *  Created on: Apr 15, 2017
 *      Author: osu8229
 */

#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[])
{
//	sleep(30);
    int rank, size;

    MPI::Init(argc, argv);
    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);

    int sizeOfArray = ((rank + 2) * 10) / 7;
//    int sizeOfArray = 4;
    cout<<"SizeOfArray on Proc "<< rank << " is :" << sizeOfArray << endl;

    double val = 100 * rank;

    double* weights = (double *) malloc(sizeOfArray * sizeof(double));
    double* allWeights = (double *) malloc(18 * sizeof(double));
//    double weights[sizeOfArray];
//    double allWeights[18];
    for(int i = 0; i < sizeOfArray; i++){
            weights[i] = 0;
        }
    for(int i = 0; i < 18; i++){
        allWeights[i] = 0;
    }

    int recvCounts[] = {2, 4, 5, 7};
    int displs[] = {0, 2, 6, 11};

    for(int i = 0; i < sizeOfArray; i++){
    	weights[i] = (rank + 1)*100;
    }

    MPI::COMM_WORLD.Barrier();

    if(rank == 0){
		for(int i = 0; i < sizeOfArray; i++){
			cout << weights[i] <<endl;
		}
    }


    MPI::COMM_WORLD.Barrier();

    if(rank == 1){
    	for(int i = 0; i < sizeOfArray; i++){
    		cout << weights[i]<<endl;
    	}
    }


    MPI::COMM_WORLD.Barrier();

        if(rank == 2){
    		for(int i = 0; i < sizeOfArray; i++){
    			cout << weights[i]<<endl;
    		}
        }


        MPI::COMM_WORLD.Barrier();

        if(rank == 3){
        	for(int i = 0; i < sizeOfArray; i++){
        		cout << weights[i]<<endl;
        	}
        }


        MPI::COMM_WORLD.Barrier();

    try{
    	MPI::COMM_WORLD.Allgatherv(weights,
    							   sizeOfArray,
    							   MPI::DOUBLE,
    							   allWeights,
    							   recvCounts,
    							   displs,
    							   MPI::DOUBLE);
    }catch(MPI::Exception failure){
    	cout <<  failure.Get_error_string() << endl;
    	MPI::COMM_WORLD.Abort(1);
    }

    MPI::COMM_WORLD.Barrier();
//    MPI::COMM_WORLD.Allgather(&val, 1, MPI::DOUBLE, &allWeights, 4, MPI::DOUBLE);

    if(rank == 2){
    	for(int i = 0; i < 18; i++){
    		cout<< "allWeights["<<i<<"] : "<<allWeights[i] << endl;
    	}
    }

    MPI::COMM_WORLD.Barrier();

//    if(rank == 1){
//		for(int i = 0; i < 8; i++){
//			cout<< "allWeights["<<i<<"] : "<<allWeights[i] << endl;
//		}
//    }
//
//    MPI::COMM_WORLD.Barrier();
    cout<<"ABCD"<<rank<<endl;
    free(weights);
    free(allWeights);
    MPI::Finalize();
    return 0;
}

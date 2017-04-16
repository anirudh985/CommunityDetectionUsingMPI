#include "mpi.h"
#include <stdio.h>

#define MAX_PROCESSES 10

int main( int argc, char **argv )
{
    int rank, size, i,j;
    int table[MAX_PROCESSES][MAX_PROCESSES];
    int errors=0;
    int participants;
    int displs[MAX_PROCESSES];
    int recv_counts[MAX_PROCESSES];

    MPI_Init( &argc, &argv );
    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    /* A maximum of MAX_PROCESSES processes can participate */
    if ( MAX_PROCESSES % size || size > MAX_PROCESSES )
    {
        fprintf( stderr, "Number of processors must divide into %d\n", MAX_PROCESSES );fflush(stderr);
        MPI::COMM_WORLD.Abort(1);
    }
    participants = size;
    if ( (rank < participants) )
    {
        /* Determine what rows are my responsibility */
        int block_size = MAX_PROCESSES / participants;
        int begin_row = rank * block_size;
        int end_row = (rank+1) * block_size;
        int send_count = block_size * MAX_PROCESSES;

        /* Fill in the displacements and recv_counts */
        for (i=0; i<participants; i++) {
            displs[i] = i * block_size * MAX_PROCESSES;
            recv_counts[i] = send_count;
        }
        /* Paint my rows my color */
        for (i=begin_row; i<end_row ;i++)
            for (j=0; j<MAX_PROCESSES; j++)
                table[i][j] = rank + 10;

        /* Everybody gets the gathered data */
        MPI::COMM_WORLD.Allgatherv(&table[begin_row][0], send_count, MPI::INT,
                            &table[0][0], recv_counts, displs,
                            MPI::INT);
        /* Everybody should have the same table now.
            The entries are:
            Table[i][j] = (i/block_size) + 10;
        */
        for (i=0; i<MAX_PROCESSES;i++)
            if ( (table[i][0] - table[i][MAX_PROCESSES-1] !=0) )
                errors++;
        for (i=0; i<MAX_PROCESSES;i++) {
            for (j=0; j<MAX_PROCESSES;j++) {
                if (table[i][j] != (i/block_size) + 10) errors++;
            }
        }
        if (errors) {
            /* Print out table if there are any errors */
            for (i=0; i<MAX_PROCESSES;i++) {
                printf("\n");
                for (j=0; j<MAX_PROCESSES; j++)
                    printf(" %d",table[i][j]);
            }
            printf("\n");
            fflush(stdout);
        }
    }
    MPI::Finalize();
    if (errors)
    {
        printf( "[%d] done with ERRORS(%d)!\n", rank, errors );fflush(stdout);
    }
    return errors;
}

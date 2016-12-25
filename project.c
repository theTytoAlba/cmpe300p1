/* Student Name: Irmak Kavasoglu
 * Student Number: 2013400090
 * Compile Status: 
 * Program Status: 
 * Notes: 
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/**
 * Reads 200*200 input from file input.txt
 * Returns a int **, which is the 200*200 integer array that is read. 
 */
int **readInput() {
    // Prepare the reader.
    FILE *reader = fopen("input.txt", "r");

    // Allocate space.
    int **matrix = malloc(200 * sizeof(int *));
    for(int i = 0; i < 200; i++) {
        matrix[i] = malloc(200 * sizeof(int));
    }

    // Read data.
    for (int i = 0; i < 200; i++) {
        for (int j = 0; j < 200; j++) {
            fscanf(reader, "%d", &matrix[i][j]);
        }
    }
    // Close the reader.
    fclose(reader);

    printf("Data read from input.txt.\n");
    return matrix;
}

/**
 * Writes out 198*198 output from file input.txt
 */
void writeOutput(int **matrix) {
    // Prepare the writer.
    FILE *writer = fopen("output.txt", "w");

    // Write data.
    for (int i = 0; i < 198; i++) {
        for (int j = 0; j < 198; j++) {
            fprintf(writer, "%d ", matrix[i][j]);
        }
        fprintf(writer, "\n");
    }
    
    // Close the writer.
    fclose(writer);
    printf("Data written to output.txt.\n");
}

int main(int argc, char* argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

    int **matrix;
    // The boss will read the input.
    if (rank == 0) {
		printf("Reading the input.\n");
        matrix = readInput();
    } else {
        //DO NOTHING FOR NOW	
    }

    //TODO: Calculate

    // The boss will write the output.
    if (rank == 0) {
        printf("Writing the output.\n");
        writeOutput(matrix);
    }
        
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}

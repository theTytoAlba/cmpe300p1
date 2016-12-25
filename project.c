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
 * Allocates a 2D integer array of given size.
 */
int **allocate2DArray(int rows, int cols) {
    int *data = (int *) malloc(rows*cols*sizeof(int));
    int **array= (int **) malloc(rows*sizeof(int*));
    for (int i = 0; i < rows; i++) {
        array[i] = &(data[cols*i]);
    }
    return array;
}

/**
 * Calculates the average of a 3*3 matrix.
 */
int average(int window[3][3]){
    double sum = 0.0;
    int j=0, i=0;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            sum += window[i][j];
        }
    }
    return sum/9.0;
}

/**
 * Reads 200*200 input from file input.txt
 * Returns a int **, which is the 200*200 integer array that is read. 
 */
int **readInput() {
    // Prepare the reader.
    FILE *reader = fopen("input.txt", "r");

    // Allocate space.
    int **matrix = allocate2DArray(200, 200);

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
 * Writes out output.
 */
void writeOutput(int **matrix, int rows, int columns, char *name) {
    // Prepare the writer.
    FILE *writer = fopen(name, "w");

    // Write data.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            fprintf(writer, "%d ", matrix[i][j]);
        }
        fprintf(writer, "\n");
    }
    
    // Close the writer.
    fclose(writer);
    printf("Data written to %s.\n", name);
}

/**
 * Processor 0 will use boss mode.
 */
void bossMode(int size) {
    // The boss will read the input.
    printf("Reading the input.\n");
    int **matrix = readInput();
    
    // The boss will send the data to the slaves.
    int slaveRows = 200/(size-1);
    printf("Sending data of %d rows to %d slaves.\n", slaveRows, size-1);
    for (int i = 1; i < size; i++) {
        MPI_Send(&matrix[(i-1)*slaveRows][0], 200*slaveRows, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
}

/**
 * Processors 1, 2, ...,n will use slave mode. 
 */
void slaveMode(int rank, int size) {
    // Allocate space.
    int rows = 200/(size-1);
    int **matrixPart = allocate2DArray(rows, 200);
    
    // Receive matrix part.
    MPI_Recv(&matrixPart[0][0], rows*200, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Process %d has received the matrix part.\n", rank);
    
    // Print received matrix part.
    char name[11];
    sprintf(name, "output_%d.txt", rank);
    writeOutput(matrixPart, rows, 200, name);
}

int main(int argc, char* argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Processor 0 is the boss and the others are slaves.
    if (rank == 0) {
        bossMode(size);
    } else {
        slaveMode(rank, size);   
    }
        
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

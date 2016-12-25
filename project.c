/* Student Name: Irmak Kavasoglu
 * Student Number: 2013400090
 * Compile Status: 
 * Program Status: 
 * Notes: 
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define DEFAULT_TAG 0
#define REQUEST_ABOVE_TAG 1
#define REQUEST_BELOW_TAG 2
#define INCOMING_LINE_TAG 3
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
int smooth(int window[3][3]){
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
        MPI_Send(&matrix[(i-1)*slaveRows][0], 200*slaveRows, MPI_INT, i, DEFAULT_TAG, MPI_COMM_WORLD);
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
    MPI_Recv(&matrixPart[0][0], rows*200, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Process %d has received the matrix part.\n", rank);
    
    // Allocate space for smoothed matrix. First and last parts have one less line.
    int smoothedRows = rows;
    if (rank == 1 || rank == size-1) {
        smoothedRows -= 2;
    }
    
    int **smoothedPart = allocate2DArray(smoothedRows, 198);
    
    // Smooth the part.
    for (int smoothedRow = 0, currentRow = (rank == 1 ? 1 : 0); smoothedRow < smoothedRows; currentRow++, smoothedRow++) {
        for (int currCol = 1; currCol < 198; currCol++) {
            // Construct 3*3 matrix to be smoothed.
            int conv[3][3];
            // Prepare request variables.
            MPI_Request requestAbove, requestBelow, otherRequests;
            int waitingAbove = 0, waitingBelow = 0;
            // Put the top three. If we need info from above, send request.
            if (currentRow == 0) {
                // Request an array of 3, starting from current column-1.
                int startColumn = currCol - 1;
                MPI_Isend(&startColumn, 1, MPI_INT, rank-1, REQUEST_ABOVE_TAG, MPI_COMM_WORLD, &requestAbove);
                printf("I am %d and requesting above from %d for column %d\n", rank, rank-1, startColumn);
                waitingAbove = 1;
            } else {
                conv[0][0] = matrixPart[currentRow-1][currCol-1];
                conv[0][1] = matrixPart[currentRow-1][currCol];
                conv[0][2] = matrixPart[currentRow-1][currCol+1];
            }
            // Put the mid three
            conv[1][0] = matrixPart[currentRow][currCol-1];
            conv[1][1] = matrixPart[currentRow][currCol];
            conv[1][2] = matrixPart[currentRow][currCol+1];
            // Put the bottom three. If we need info from below, send request.
            
            if (currentRow == rows-1) {
                // Request an array of 3, starting from current column-1.
                int startColumn = currCol - 1;
                MPI_Isend(&startColumn, 1, MPI_INT, rank+1, REQUEST_BELOW_TAG, MPI_COMM_WORLD, &requestBelow);
                printf("I am %d and requesting below from %d\n", rank, rank+1);
                waitingBelow = 1;
            } else {
                conv[2][0] = matrixPart[currentRow+1][currCol-1];
                conv[2][1] = matrixPart[currentRow+1][currCol];
                conv[2][2] = matrixPart[currentRow+1][currCol+1];
            }

            // Wait until conv matrix has all needed info.
            while(waitingAbove || waitingBelow) {
                MPI_Status status;
                int received[3];
                MPI_Recv(&received, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == INCOMING_LINE_TAG) {
                    printf("Received incoming line from %d by %d\n",status.MPI_SOURCE, rank);
                    // Its the above line if it is from rank-1.
                    if (status.MPI_SOURCE == rank-1) {
                        conv[0][0] = received[0];
                        conv[0][1] = received[1];
                        conv[0][2] = received[2];
                        waitingAbove = 0;
                    }
                    // Its the below line if it is from rank+1.
                    if (status.MPI_SOURCE == rank+1) {
                        conv[2][0] = received[0];
                        conv[2][1] = received[1];
                        conv[2][2] = received[2];
                        waitingBelow = 0;
                    }
                } else if (status.MPI_TAG == REQUEST_ABOVE_TAG) {
                    printf("Received request above by %d, sending to %d\n", rank, rank+1);
                    // Processor below is making a request.
                    int startColumn = received[0];
                    int requested[3];
                    requested[0] = matrixPart[rows-1][startColumn];
                    requested[1] = matrixPart[rows-1][startColumn+1];
                    requested[2] = matrixPart[rows-1][startColumn+2];
                    MPI_Isend(&requested, 3, MPI_INT, rank+1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
                } else if (status.MPI_TAG == REQUEST_BELOW_TAG) {
                    printf("Received request below by %d, sending to %d\n", rank, rank-1);
                    // Processor above is making a request.
                    int startColumn = received[0];
                    int requested[3];
                    requested[0] = matrixPart[0][startColumn];
                    requested[1] = matrixPart[0][startColumn+1];
                    requested[2] = matrixPart[0][startColumn+2];
                    MPI_Isend(&requested, 3, MPI_INT, rank-1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
                }
            }

            // Now the conv matrix is complete, do the smoothing.
            smoothedPart[smoothedRow][currCol] = smooth(conv);
        }
    }

    // Print given part of the matrix.
    char name[11];
    sprintf(name, "output_%d.txt", rank);
    writeOutput(matrixPart, rows, 200, name);

    // Print smoothed version of it.
    char names[11];
    sprintf(names, "outputsm_%d.txt", rank);
    writeOutput(smoothedPart, smoothedRows, 198, names);

    //TODO: Boss should be acknowleged and send OK signal.
    printf("Now waiting for boss to tell me its OK to continue.\n");
    int bossSaidOk = 0;
    while (!bossSaidOk) {
        MPI_Status status;
        MPI_Request otherRequests;
        int received;
        MPI_Recv(&received, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == REQUEST_ABOVE_TAG) {
            printf("Received request above by %d, sending to %d\n", rank, rank+1);
            // Processor below is making a request.
            int startColumn = received;
            int requested[3];
            requested[0] = matrixPart[rows-1][startColumn];
            requested[1] = matrixPart[rows-1][startColumn+1];
            requested[2] = matrixPart[rows-1][startColumn+2];
            MPI_Isend(&requested, 3, MPI_INT, rank+1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
        } else if (status.MPI_TAG == REQUEST_BELOW_TAG) {
            printf("Received request below by %d, sending to %d\n", rank, rank-1);
            // Processor above is making a request.
            int startColumn = received;
            int requested[3];
            requested[0] = matrixPart[0][startColumn];
            requested[1] = matrixPart[0][startColumn+1];
            requested[2] = matrixPart[0][startColumn+2];
            MPI_Isend(&requested, 3, MPI_INT, rank-1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
        }
    }
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

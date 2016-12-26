/* Student Name: Irmak Kavasoglu
 * Student Number: 2013400090
 * Compile Status: Compiling
 * Program Status: Working
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define DEFAULT_TAG 0
#define REQUEST_ABOVE_TAG 1
#define REQUEST_BELOW_TAG 2
#define INCOMING_LINE_TAG 3
#define SMOOTHING_RESULT_TAG 4
#define THRESHOLDING_RESULT_TAG 5
#define BOSS_RECEIVED_SMOOTH 6
#define BOSS_RECEIVED_THRESHOLDED 7
#define SMOOTHING_PROCESS 8
#define THRESHOLDING_PROCESS 9
#define COLLECT_SMOOTHED 10
#define COLLECT_THRESHOLDED 11
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
 * Deallocates the allocated array.
 */
void deallocate2DArray(int **array) {
    free(*array);
    free(array); 
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
 * Returns 0 if the window does not exceed given threshold, 255 otherwise.
 */
int threshold(int window[3][3], int thresholdNumber) {
    int calc[4];
    // Calculate horizontal
    calc[0] = 0;
    for (int i = 0; i < 3; i++) {
        calc[0] -= window[0][i];
        calc[0] += 2*window[1][i];
        calc[0] -= window[2][i];
    }
    // Calculate vertical
    calc[1] = 0;
    for (int i = 0; i < 3; i++) {
        calc[1] -= window[i][0];
        calc[1] += 2*window[i][1];
        calc[1] -= window[i][2];
    }
    // Calculate +45
    calc[2] = 0;
    for (int i = 0; i < 3; i++) {
        calc[2] -= window[(i+1)%3][2-i];
        calc[2] += 2*window[i][2-i];
        calc[2] -= window[i][(4-i)%3];
    }
    // Calculate -45
    calc[3] = 0;
    for (int i = 0; i < 3; i++) {
        calc[3] -= window[(i+1)%3][i];
        calc[3] += 2*window[i][i];
        calc[3] -= window[i][(i+1)%3];
    }

    //threshold
    int result = 0;
    for (int i = 0; i < 4; i++) {
        if (calc[i] > thresholdNumber) {
            result = 255;
            break;
        }
    }
    return result;
}

/**
 * Reads 200*200 input from file input.txt to the given int **. 
 */
void readInput(char *input, int **matrix) {
    // Prepare the reader.
    FILE *reader = fopen(input, "r");

    // Read data.
    for (int i = 0; i < 200; i++) {
        for (int j = 0; j < 200; j++) {
            fscanf(reader, "%d", &matrix[i][j]);
        }
    }
    // Close the reader.
    fclose(reader);
    printf("Data read from %s.\n", input);
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
 * Can be function in two ways depending on the tag:
 * - Collector of the smoothed parts if tag is COLLECT_SMOOTHED
 * - Collector of the thresholded parts if tag is COLLECT_THRESHOLDED
 *
 * Does not return until it collects all the parts from every slave processor.
 * Puts the result in the given result matrice.
 */
void collectParts(int size, int tag, int **result) {
    MPI_Status status;

    int slaveRows = 200/(size-1);
    int columns;
    int **resultPart;
    // First and last processor's row number is less. Lack keeps number of missing lines.
    int lack;
    // Fill the variables depending on the tag.
    if (tag == COLLECT_SMOOTHED) {
        resultPart = allocate2DArray(slaveRows, 198);
        columns = 198;
        lack = 1;
    } else if (tag == COLLECT_THRESHOLDED) {
        resultPart = allocate2DArray(slaveRows, 196);
        columns = 196;
        lack = 2;
    }
    
    int partSize = slaveRows*columns;
    
    for (int i = 0; i < size-1; i++) {
        int resultTag = (tag == COLLECT_SMOOTHED) ? SMOOTHING_RESULT_TAG : THRESHOLDING_RESULT_TAG;
        MPI_Recv(&resultPart[0][0], partSize, MPI_INT, MPI_ANY_SOURCE, resultTag, MPI_COMM_WORLD, &status);
        //printf("here %d receiving %d.\n", i, status.MPI_SOURCE);
        // if it is the first or last one, take slaveRows-lack line.
        if (status.MPI_SOURCE == 1) {
            for (int j = 0; j < slaveRows-lack; j++) {
                for (int k = 0; k < columns; k++) {
                    result[j][k] = resultPart[j][k];
                }
            }
        } else if (status.MPI_SOURCE == size-1) {
            for (int j = 0; j < slaveRows-lack; j++) {
                for (int k = 0; k < columns; k++) {
                    int resultRowStart = (status.MPI_SOURCE-1)*slaveRows - lack;
                    result[resultRowStart + j][k] = resultPart[j][k];
                }
            }
        } else {
            for (int j = 0; j < slaveRows; j++) {
                for (int k = 0; k < columns; k++) {
                    int resultRowStart = (status.MPI_SOURCE-1)*slaveRows - lack;
                    // if there are 200 processors this is possible for 2nd and 199th processors.
                    if (resultRowStart < 0 || resultRowStart >= columns) {
                        continue;
                    }
                    result[resultRowStart + j][k] = resultPart[j][k];
                }
            }
        }
    }

    deallocate2DArray(resultPart);
}

/**
 * Processor 0 will use boss mode.
 * Distributes the matrix parts to slaves.
 * Collects and merges smoothed and thresholded images.
 */
void bossMode(int size, char *input, char *output) {
    // The boss will read the input.
    printf("Reading the input.\n");
    int **matrix = allocate2DArray(200, 200);
    readInput(input, matrix);
    
    // The boss will send the data to the slaves.
    int slaveRows = 200/(size-1);
    printf("Sending data of %d rows to %d slaves.\n", slaveRows, size-1);
    for (int i = 1; i < size; i++) {
        MPI_Send(&matrix[(i-1)*slaveRows][0], 200*slaveRows, MPI_INT, i, DEFAULT_TAG, MPI_COMM_WORLD);
    }

    // Deallocate matrix.
    deallocate2DArray(matrix);

    // Allocate space for smoothed image.
    int **smoothed = allocate2DArray(198, 198);

    // Receive the smoothed parts from everyone and merge. 
    collectParts(size, COLLECT_SMOOTHED, smoothed);
    printf("collected smoothed\n");

    // Write smoothed output if you want.
    // writeOutput(smoothed, 198, 198, "smoothed.txt");

    // Tell slaves to stop waiting.
    for (int i = 1; i < size; i++) {
        MPI_Send(&i, 1, MPI_INT, i, BOSS_RECEIVED_SMOOTH, MPI_COMM_WORLD);
    }

    // Deallocate.
    deallocate2DArray(smoothed);
    
    // Allocate space for thresholded image.
    int **thresholded = allocate2DArray(196, 196);    

    // Receive the thresholded parts from everyone and merge. 
    collectParts(size, COLLECT_THRESHOLDED, thresholded);
    printf("collected thresholded\n");

    // Tell slaves to stop waiting.
    for (int i = 1; i < size; i++) {
       MPI_Send(&i, 1, MPI_INT, i, BOSS_RECEIVED_THRESHOLDED, MPI_COMM_WORLD);
    }
    // Write output and deallocate.
    writeOutput(thresholded, 196, 196, output);
    deallocate2DArray(thresholded);
}

/**
 * Can be called for both smoothing and thresholding processes
 * by using tags SMOOTHING_PROCESS and THRESHOLDING_PROCESS.
 *
 * Performs the given process on the given matrix and returns the result.
 * Communicates with other processes for the missing window parts if needed.
 */
int **processPart(int rank, int size, int **matrixPart, int tag, int thresholdNumber) {
    MPI_Status status;
    MPI_Request otherRequests;

    // Find number of rows for the given matrix
    int rows = 200/(size-1), columns;
    int resultRows = rows;
    
    if (tag == SMOOTHING_PROCESS) {
        columns = 200;    
        if ((rank == 1) || (rank == (size-1))) {
            resultRows -= 1;
        }
    } else if (tag == THRESHOLDING_PROCESS) {
        columns = 198;
        if ((rank == 1) || (rank == (size-1))) {
            rows -= 1;
            resultRows -= 2;
            // If there are 200 processors this is possible.
            if (resultRows < 0) {
                resultRows = 0;
            }
        }
    }
    
    // Allocate space for result matrix. First and last parts have one less line.
    int **processResult = allocate2DArray(resultRows, columns-2);

    // Smooth the part.
    for (int resultRow = 0, currentRow = (rank == 1 ? 1 : 0); resultRow < resultRows; currentRow++, resultRow++) {
        for (int currCol = 1; currCol < columns-1; currCol++) {
            // Prepare request variables.
            MPI_Request requestAbove, requestBelow;
            int waitingAbove = 0, waitingBelow = 0;
            
            // Construct 3*3 matrix to be processed.
            int conv[3][3];
            
            // Put the top three. If we need info from above, send request.
            if (currentRow == 0) {
                int startColumn = currCol - 1;
                MPI_Isend(&startColumn, 1, MPI_INT, rank-1, REQUEST_ABOVE_TAG, MPI_COMM_WORLD, &requestAbove);
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
                waitingBelow = 1;
            } else {
                conv[2][0] = matrixPart[currentRow+1][currCol-1];
                conv[2][1] = matrixPart[currentRow+1][currCol];
                conv[2][2] = matrixPart[currentRow+1][currCol+1];
            }

            // Wait until conv matrix has all needed info.
            while(waitingAbove || waitingBelow) {
                int received[3];
                MPI_Recv(&received, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == INCOMING_LINE_TAG) {
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
                    // Processor below is making a request.
                    int startColumn = received[0];
                    int requested[3];
                    requested[0] = matrixPart[rows-1][startColumn];
                    requested[1] = matrixPart[rows-1][startColumn+1];
                    requested[2] = matrixPart[rows-1][startColumn+2];
                    MPI_Isend(&requested, 3, MPI_INT, rank+1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
                } else if (status.MPI_TAG == REQUEST_BELOW_TAG) {
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
            if (tag == SMOOTHING_PROCESS) {
                processResult[resultRow][currCol-1] = smooth(conv);
            } else if (tag == THRESHOLDING_PROCESS) {
                processResult[resultRow][currCol-1] = threshold(conv, thresholdNumber);
            }
        }
    }
    return processResult;
}

/**
 * Waits until the boss gives the signal that is the same as the given tag.
 * If it receives requests form other slaves, sends them the parts they want.
 */
void waitForBoss(int rank, int size, int **currentPart, int tag) {
    MPI_Status status;
    MPI_Request otherRequests;

    int rows = 200/(size-1);

    // First and last processor's row number is less. Lack keeps number of missing lines.
    int lack = 0;
    if ((rank == 1 || rank == (size-1)) && tag == COLLECT_THRESHOLDED) {
        lack = 1;
    }

    if ((rank == 1) || (rank == (size-1))) {
      //  processedRows -= lack;
    }

    int bossSaidOk = 0;
    while (!bossSaidOk) {
        int startColumn;
        MPI_Recv(&startColumn, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == REQUEST_ABOVE_TAG) {
            // Processor below is making a request.
            int requested[3];
            requested[0] = currentPart[rows-lack-1][startColumn];
            requested[1] = currentPart[rows-lack-1][startColumn+1];
            requested[2] = currentPart[rows-lack-1][startColumn+2];
            MPI_Isend(&requested, 3, MPI_INT, rank+1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
        } else if (status.MPI_TAG == REQUEST_BELOW_TAG) {
            // Processor above is making a request.
            int requested[3];
            requested[0] = currentPart[0][startColumn];
            requested[1] = currentPart[0][startColumn+1];
            requested[2] = currentPart[0][startColumn+2];
            MPI_Isend(&requested, 3, MPI_INT, rank-1, INCOMING_LINE_TAG, MPI_COMM_WORLD, &otherRequests);
        } else if (status.MPI_TAG == tag) {
            bossSaidOk = 1;
        }
    }
}

/**
 * Processors 1, 2, ...,n will use slave mode. 
 */
void slaveMode(int rank, int size, int thresholdNumber) {
    MPI_Status status;
    MPI_Request otherRequests;

    // Allocate space.
    int rows = 200/(size-1);
    int **matrixPart = allocate2DArray(rows, 200);
    
    // Receive matrix part.
    MPI_Recv(&matrixPart[0][0], rows*200, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // Find the number of rows for this processor's part.
    int smoothedRows = rows;
    if ((rank == 1) || (rank == (size-1))) {
        smoothedRows -= 1;
    }

    // Calculate smoothed part, send the result to the boss and wait for its approval before continuing.
    int **smoothedPart = processPart(rank, size, matrixPart, SMOOTHING_PROCESS, thresholdNumber);
    MPI_Isend(&smoothedPart[0][0], smoothedRows*198, MPI_INT, 0, SMOOTHING_RESULT_TAG, MPI_COMM_WORLD, &otherRequests);
    waitForBoss(rank, size, matrixPart, BOSS_RECEIVED_SMOOTH);

    // After boss said continue, slave no longer needs the original matrix part.
    deallocate2DArray(matrixPart);

    // Find the number of rows for this process.
    int thresholdRows = rows;
    if ((rank == 1) || (rank == (size-1))) {
        thresholdRows -= 2;
        // This is possible if there are 200 processors.
        if (thresholdRows < 0) {
            thresholdRows = 0;
        }
    }
    int **thresholdPart = processPart(rank, size, smoothedPart, THRESHOLDING_PROCESS, thresholdNumber); 
    MPI_Isend(&thresholdPart[0][0], thresholdRows*196, MPI_INT, 0, THRESHOLDING_RESULT_TAG, MPI_COMM_WORLD, &otherRequests);
    waitForBoss(rank, size, smoothedPart, BOSS_RECEIVED_THRESHOLDED);

    // Deallocate remaining matrices.
    deallocate2DArray(smoothedPart);
    deallocate2DArray(thresholdPart);
}


/**
 * Program starts from main function.
 * Main function redirects processors to the relevant function,
 * depending on if they are slaves or the boss.
 */
int main(int argc, char* argv[])
{
    // Get rank and size.
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Processor 0 is the boss and the others are slaves.
    if (rank == 0) {
        // Send input and output file names to the boss.
        char *input = argv[1];
        char *output = argv[2];
        bossMode(size, input, output);
    } else {
        // Send threshold number to each slave.
        int thresholdNumber = atoi(argv[3]);
        slaveMode(rank, size, thresholdNumber);   
    }
        
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

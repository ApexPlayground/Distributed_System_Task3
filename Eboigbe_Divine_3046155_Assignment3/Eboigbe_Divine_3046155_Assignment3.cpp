#include <iostream>
#include <mpi.h>
#include <string>

// Helper methods

void printArray(int array[2][2]) {
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << array[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int calculateInvDeterminant(int matrix[2][2]) {
    int determinant = (matrix[0][0] * matrix[1][1] - (matrix[0][1] * matrix[1][0])) % 26;
    if (determinant < 0) determinant += 26;

    int inverseDeter = -1; 
    for (int i = 0; i < 26; ++i) {
        if ((determinant * i) % 26 == 1) {
            inverseDeter = i;
            break;
        }
    }
    return inverseDeter;
}

void calculateInverse(int matrix[2][2], int inverseDeterminant, int inverseMatrix[2][2]) {
    inverseMatrix[0][0] = (matrix[1][1] * inverseDeterminant) % 26;
    inverseMatrix[0][1] = (-matrix[0][1] * inverseDeterminant) % 26;
    if (inverseMatrix[0][1] < 0) inverseMatrix[0][1] += 26;

    inverseMatrix[1][0] = (-matrix[1][0] * inverseDeterminant) % 26;
    if (inverseMatrix[1][0] < 0) inverseMatrix[1][0] += 26;

    inverseMatrix[1][1] = (matrix[0][0] * inverseDeterminant) % 26;

    if ((matrix[0][0] * inverseMatrix[0][0] + matrix[0][1] * inverseMatrix[1][0]) % 26 == 1 &&
        (matrix[0][0] * inverseMatrix[0][1] + matrix[0][1] * inverseMatrix[1][1]) % 26 == 0 &&
        (matrix[1][0] * inverseMatrix[0][0] + matrix[1][1] * inverseMatrix[1][0]) % 26 == 0 &&
        (matrix[1][0] * inverseMatrix[0][1] + matrix[1][1] * inverseMatrix[1][1]) % 26 == 1) {
        std::cout << "The inverse matrix is correct" << std::endl;
    }
    else {
        std::cout << "The inverse matrix is wrong" << std::endl;
    }
}

void decodeBlock(int block[2], int inverseMatrix[2][2]) {
    int* result = new int[2];

    result[0] = (inverseMatrix[0][0] * block[0] + inverseMatrix[0][1] * block[1]) % 26;
    result[1] = (inverseMatrix[1][0] * block[0] + inverseMatrix[1][1] * block[1]) % 26;

    for (int i = 0; i < 2; ++i) {
        if (result[i] < 0) {
            result[i] += 26;
        }
    }

    block[0] = result[0];
    block[1] = result[1];
}

// Coordinator method
void coordinator() {
    int matrix[2][2] = { {0, 0}, {0, 0} };
    int inverseMatrix[2][2];

    // Read matrix from the user
    std::cout << "Enter the 2x2 matrix elements, one by one:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << "Enter element at position (" << i << "," << j << "): ";

            while (!(std::cin >> matrix[i][j])) {
                // Clear input buffer in case of invalid input
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cout << "Invalid input. Please enter a valid integer: ";
            }
        }
    }

    std::cout << "Matrix input completed." << std::endl;
    std::cout << "Matrix:" << std::endl;
    printArray(matrix);

    // Read the message from console or command line
    std::string ciphertext;
    std::cout << "Enter the ciphertext: ";
    std::cin >> ciphertext;



    // Convert message to uppercase
    for (int i = 0; i < ciphertext.length(); ++i) {
        ciphertext[i] = std::toupper(ciphertext[i]);
    }

    // Convert characters to numbers
    int* numbers = new int[ciphertext.length() + 1];

    for (int i = 0; i < ciphertext.length(); ++i) {
        numbers[i] = ciphertext[i] - 'A';
    }
    numbers[ciphertext.length()] = 0;

    int blockSize = 2;

   
    int localBlock[2];
    MPI_Scatter(numbers, blockSize, MPI_INT, localBlock, blockSize, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the local block in the coordinator
    std::cout << "Coordinator - Local block: " << localBlock[0] << " " << localBlock[1] << std::endl;

    //// Decode the local block using the decodeBlock helper method
    //decodeBlock(localBlock, inverseMatrix);

    //// Print the decoded block (you can modify this part as needed)
    //std::cout <<" Decoded block: " << localBlock[0] << " " << localBlock[1] << std::endl;
    //

    // Print the elements of the numbers array
    std::cout << "Interger array elements: ";
    for (int i = 0; i < ciphertext.length(); ++i) {
        std::cout << numbers[i] << " ";
    }
    std::cout << std::endl;

    

    int inverseDeterminant = calculateInvDeterminant(matrix);
    std::cout << "Inverse determinant: " << inverseDeterminant << std::endl;

    /*std::cout << " " << std::endl;*/

    calculateInverse(matrix, inverseDeterminant, inverseMatrix);

    std::cout << "Inverse matrix:" << std::endl;
    printArray(inverseMatrix);

    // Broadcast the inverse matrix to all other nodes
    MPI_Bcast(&inverseMatrix[0][0], 4, MPI_INT, 0, MPI_COMM_WORLD);

    // Add MPI Barrier
    MPI_Barrier(MPI_COMM_WORLD);

    std::cout << std::endl;

    delete[] numbers;
}

// Participant method
void participant() {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int inverseMatrix[2][2] = { {0, 0}, {0, 0} };



    // Add MPI Barrier
    MPI_Barrier(MPI_COMM_WORLD);

    // Receive the inverse matrix broadcasted by the coordinator
    MPI_Bcast(&inverseMatrix[0][0], 4, MPI_INT, 0, MPI_COMM_WORLD);

    std::cout << "Rank " << world_rank << " Received Inverse matrix:" << std::endl;
    printArray(inverseMatrix);

    // Receive the local block from the coordinator
    int localBlock[2];
    MPI_Scatter(nullptr, 0, MPI_INT, localBlock, 2, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the received local block in the participant
    std::cout << "Rank " << world_rank << " - Received local block: " << localBlock[0] << " " << localBlock[1] << std::endl;

    

    // Decode the local block using the decodeBlock helper method
    decodeBlock(localBlock, inverseMatrix);

    // Print the decoded block (you can modify this part as needed)
    std::cout << "Rank " << world_rank << " Decoded block: " << localBlock[0] << " " << localBlock[1] << std::endl;
}

int main(int argc, char** argv) {
    int world_size, world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
        coordinator();
    }
    else {
        participant();
    }

    MPI_Finalize();
    return 0;
}

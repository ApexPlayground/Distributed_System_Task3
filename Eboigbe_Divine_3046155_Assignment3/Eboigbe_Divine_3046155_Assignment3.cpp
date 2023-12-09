#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <string>

// Helper method to print an array to the console
void printArray(const int* array, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}

// Helper method to decode a block using an inverse matrix
void decodeBlock(const int* block, const int* inverseMatrix, int* decodedBlock) {
    // Perform matrix multiplication and modulo 26
    decodedBlock[0] = (block[0] * inverseMatrix[0] + block[1] * inverseMatrix[1]) % 26;
    decodedBlock[1] = (block[0] * inverseMatrix[2] + block[1] * inverseMatrix[3]) % 26;

    // Ensure the result is non-negative
    if (decodedBlock[0] < 0) {
        decodedBlock[0] += 26;
    }
    if (decodedBlock[1] < 0) {
        decodedBlock[1] += 26;
    }
}

// Helper method to calculate the modulo 26 inverse determinant
int calculateInvDeterminant(int determinant) {
    for (int i = 0; i < 26; ++i) {
        if ((determinant * i) % 26 == 1) {
            return i;
        }
    }
    // If no inverse exists, return an error value (e.g., -1)
    return -1;
}

// Helper method to calculate the modulo 26 inverse of a 2x2 matrix
void calculateInverse(const int* matrix, int* inverseMatrix) {
    int determinant = (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % 26;

    // Ensure the determinant is non-negative
    if (determinant < 0) {
        determinant += 26;
    }

    // Calculate the inverse determinant
    int invDeterminant = calculateInvDeterminant(determinant);

    // Calculate the elements of the inverse matrix
    inverseMatrix[0] = (matrix[3] * invDeterminant) % 26;
    inverseMatrix[1] = (-matrix[1] * invDeterminant) % 26;
    inverseMatrix[2] = (-matrix[2] * invDeterminant) % 26;
    inverseMatrix[3] = (matrix[0] * invDeterminant) % 26;

    // Ensure the elements are non-negative
    for (int i = 0; i < 4; ++i) {
        if (inverseMatrix[i] < 0) {
            inverseMatrix[i] += 26;
        }
    }
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank and size of each process
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Rank 0 acts as the coordinator, while other ranks are participants
    if (rank == 0) {
        // Coordinator code here
        std::cout << "Coordinator (Rank " << rank << ") initialized.\n";

        // a) Coordinator reads the matrix, calculates the inverse determinant, and broadcasts it.
        int matrix[4];
        std::cout << "Enter a 2x2 matrix: ";
        for (int i = 0; i < 4; ++i) {
            std::cin >> matrix[i];
        }

        int determinant = (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % 26;
        if (determinant < 0) {
            determinant += 26;
        }

        int invDeterminant = calculateInvDeterminant(determinant);
        MPI_Bcast(&invDeterminant, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // b) Coordinator calculates the inverse matrix (A-1).
        int invMatrix[4];
        calculateInverse(matrix, invMatrix);

        // c) Coordinator checks the correctness of the inverse matrix.
        // (Optional: Print the inverse matrix for verification)
        std::cout << "Inverse Matrix: ";
        printArray(invMatrix, 4);

        // d) Coordinator broadcasts the inverse matrix to all other nodes.
        MPI_Bcast(invMatrix, 4, MPI_INT, 0, MPI_COMM_WORLD);

        // e) Coordinator reads the ciphertext and distributes it among all nodes.
        std::string ciphertext;
        std::cout << "Enter the ciphertext: ";
        std::cin >> ciphertext;
        int cipherLength = ciphertext.length();
        MPI_Bcast(&cipherLength, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int* cipherArray = new int[cipherLength];
        for (int i = 0; i < cipherLength; ++i) {
            cipherArray[i] = ciphertext[i] - 'A'; // Convert character to corresponding number
        }
        MPI_Scatter(cipherArray, 1, MPI_INT, cipherArray, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // g) Coordinator collects decrypted blocks from participants.
        int* decryptedBlocks = new int[cipherLength];
        MPI_Gather(decryptedBlocks, 1, MPI_INT, decryptedBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // h) Coordinator converts the decrypted message to characters.
        std::string decryptedMessage;
        for (int i = 0; i < cipherLength; ++i) {
            decryptedMessage += static_cast<char>(decryptedBlocks[i] + 'A'); // Convert number to character
        }

        // i) Coordinator outputs the deciphered text.
        std::cout << "Deciphered Text: " << decryptedMessage << std::endl;

        delete[] cipherArray;
        delete[] decryptedBlocks;
    }
    else {
        // Participant code here
        std::cout << "Participant (Rank " << rank << ") initialized.\n";

        // a) Participants receive the inverse determinant broadcasted by the coordinator.
        int invDeterminant;
        MPI_Bcast(&invDeterminant, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // b) Participants receive the inverse matrix broadcasted by the coordinator.
        int invMatrix[4];
        MPI_Bcast(invMatrix, 4, MPI_INT, 0, MPI_COMM_WORLD);

        // e) Participants receive the length of the ciphertext broadcasted by the coordinator.
        int cipherLength;
        MPI_Bcast(&cipherLength, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // e) Participants receive their portion of the ciphertext from the coordinator.
        int* receiveBuffer = new int[cipherLength / size];
        MPI_Scatter(nullptr, 0, MPI_INT, receiveBuffer, cipherLength / size, MPI_INT, 0, MPI_COMM_WORLD);

        // Use the receiveBuffer in the decoding process
        int* decodedBlocks = new int[cipherLength / size * 2];  // Each element of receiveBuffer corresponds to 2 decoded values
        decodeBlock(receiveBuffer, invMatrix, decodedBlocks);

        // g) Participants send their decrypted blocks to the coordinator.
        MPI_Gather(decodedBlocks, cipherLength / size * 2, MPI_INT, nullptr, 0, MPI_INT, 0, MPI_COMM_WORLD);

        delete[] receiveBuffer;
        delete[] decodedBlocks;
    }


    // Finalize MPI
    MPI_Finalize();

    return 0;
}

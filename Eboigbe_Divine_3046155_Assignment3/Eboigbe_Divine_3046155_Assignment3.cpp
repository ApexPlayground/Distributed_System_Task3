#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <string>

// Helper method to print an array
void printArray(int* array, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}

// Helper method to calculate the modulo 26 inverse determinant
int calculateInvDeterminant(int a, int b, int c, int d) {
    // Calculate the determinant
    int determinant = (a * d - b * c) % 26;

    // Ensure the determinant is positive
    determinant = (determinant + 26) % 26;

    // Find the modulo 26 inverse from the given table
    for (int i = 0; i < 26; ++i) {
        if ((determinant * i) % 26 == 1) {
            return i;
        }
    }

    // Handle the case where no inverse exists (shouldn't happen for Hill Cipher)
    return -1;
}

// Helper method to calculate the modulo 26 inverse of a 2x2 matrix
void calculateInverse(int a, int b, int c, int d, int& invA, int& invB, int& invC, int& invD) {
    // Calculate the determinant
    int determinant = (a * d - b * c) % 26;

    // Ensure the determinant is positive
    determinant = (determinant + 26) % 26;

    // Find the modulo 26 inverse of the determinant
    int invDeterminant = 0;
    for (int i = 0; i < 26; ++i) {
        if ((determinant * i) % 26 == 1) {
            invDeterminant = i;
            break;
        }
    }

    // Calculate the entries of the inverse matrix
    invA = (d * invDeterminant) % 26;
    invB = ((-b) * invDeterminant) % 26;
    invC = ((-c) * invDeterminant) % 26;
    invD = (a * invDeterminant) % 26;

    // Ensure all entries are positive
    invA = (invA + 26) % 26;
    invB = (invB + 26) % 26;
    invC = (invC + 26) % 26;
    invD = (invD + 26) % 26;
}


// Helper method to decode a block using the inverse matrix
void decodeBlock(int block[2], int invA, int invB, int invC, int invD) {
    int x = block[0];
    int y = block[1];

    // Perform matrix multiplication with the inverse matrix
    int decodedX = (invA * x + invB * y) % 26;
    int decodedY = (invC * x + invD * y) % 26;

    // Ensure the result is positive
    decodedX = (decodedX + 26) % 26;
    decodedY = (decodedY + 26) % 26;

    // Update the block with the decoded values
    block[0] = decodedX;
    block[1] = decodedY;
}

// Helper method to convert a character to its corresponding number
int charToNumber(char c) {
    if ('A' <= c && c <= 'Z') {
        return c - 'A';
    }
    else if ('a' <= c && c <= 'z') {
        return c - 'a';
    }
    else {
        // Handle non-alphabetic characters
        return -1;
    }
}




int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank == 0) {
        // Coordinator's code
        std::cout << "I am the coordinator with rank " << world_rank << std::endl;

        // Read the matrix from the user
        int matrix[4];  // Assuming a 2x2 matrix

        if (world_rank == 0) {
            std::cout << "Enter the elements of the 2x2 matrix:" << std::endl;
            for (int i = 0; i < 4; ++i) {
                std::cin >> matrix[i];
            }
        }

        // Broadcast the matrix to all participants
        MPI_Bcast(matrix, 4, MPI_INT, 0, MPI_COMM_WORLD);

        // Calculate the modulo 26 inverse determinant
        int invDeterminant = calculateInvDeterminant(matrix[0], matrix[1], matrix[2], matrix[3]);

        // Calculate the modulo 26 matrix inverse
        int invA, invB, invC, invD;
        calculateInverse(matrix[0], matrix[1], matrix[2], matrix[3], invA, invB, invC, invD);

        // Output the result
        if (invDeterminant != -1) {
            std::cout << "Mod 26 Inverse Determinant: " << invDeterminant << std::endl;
            std::cout << "Mod 26 Matrix Inverse (A^(-1)):" << std::endl;
            std::cout << invA << " " << invB << std::endl;
            std::cout << invC << " " << invD << std::endl;

            // Check A * A^(-1) = I
            int identity[4] = { 1, 0, 0, 1 };  // Identity matrix

            int result[4];
            result[0] = (matrix[0] * invA + matrix[1] * invC) % 26;
            result[1] = (matrix[0] * invB + matrix[1] * invD) % 26;
            result[2] = (matrix[2] * invA + matrix[3] * invC) % 26;
            result[3] = (matrix[2] * invB + matrix[3] * invD) % 26;

            // Output the result of the check
            std::cout << "A * A^(-1) = I Check:" << std::endl;
            std::cout << result[0] << " " << result[1] << std::endl;
            std::cout << result[2] << " " << result[3] << std::endl;

            // Broadcast the inverse matrix to all participants
            MPI_Bcast(&invA, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&invB, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&invC, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&invD, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Read the ciphertext from the console
            std::string ciphertext;
            std::cout << "Enter the ciphertext:" << std::endl;
            std::cin >> ciphertext;

            // Convert each character into its corresponding number and save as an array of integers
            int* ciphertextNumbers = new int[ciphertext.length()];
            for (size_t i = 0; i < ciphertext.length(); ++i) {
                char c = ciphertext[i];
                int number = charToNumber(c);
                if (number != -1) {
                    ciphertextNumbers[i] = number;
                }
            }

            // Output the array of integers
            std::cout << "Ciphertext as array of integers:" << std::endl;
            printArray(ciphertextNumbers, ciphertext.length());

            // Split the array into blocks of size 2
            int blockSize = 2;
            int numBlocks = ciphertext.length() / blockSize;
            int remainder = ciphertext.length() % blockSize;

            // Determine the size of the array for each participant
            int participantBlockSize = numBlocks / (world_size - 1);
            int participantRemainder = numBlocks % (world_size - 1);

            // Send the size of the array for each participant
            MPI_Bcast(&participantBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&participantRemainder, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // Send the array to each participant
            for (int i = 1; i < world_size; ++i) {
                int startIdx = (i - 1) * participantBlockSize * blockSize;
                int endIdx = startIdx + participantBlockSize * blockSize - 1;

                // Adjust for remainder
                if (i <= participantRemainder) {
                    startIdx += (i - 1) * blockSize;
                    endIdx += blockSize;
                }
                else {
                    startIdx += participantRemainder * blockSize;
                    endIdx += participantRemainder * blockSize;
                }

                // Send the size of the block
                int blockSizeToSend = endIdx - startIdx + 1;
                MPI_Send(&blockSizeToSend, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

                // Send the block of integers to the participant
                MPI_Send(&ciphertextNumbers[startIdx], blockSizeToSend, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

            // Receive the size of the received block
            int receivedBlockSize;
            MPI_Recv(&receivedBlockSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Receive the block of integers
            int* receivedBlock = new int[receivedBlockSize];
            MPI_Recv(receivedBlock, receivedBlockSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            // Decode the received block using the inverse matrix
            for (int i = 0; i < blockSize; i += 2) {
                int block[2] = { receivedBlock[i], receivedBlock[i + 1] };
                decodeBlock(block, invA, invB, invC, invD);
                receivedBlock[i] = block[0];
                receivedBlock[i + 1] = block[1];
            }

            // Output the decoded block for each participant
            std::cout << "Decoded block for participant " << world_rank << ": ";
            printArray(receivedBlock, blockSize);

            // Clean up the dynamically allocated array
            delete[] receivedBlock;

            // Receive the size of the received block for each participant
            int receivedParticipantBlockSize;
            MPI_Bcast(&receivedParticipantBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


            // Allocate space to store all decrypted blocks
            int* decryptedBlocks = new int[ciphertext.length()];

            // Collect decrypted blocks from each participant
            for (int i = 1; i < world_size; ++i) {
                // Receive the size of the block
                int blockSize;
                MPI_Recv(&blockSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Receive the block of decrypted integers
                int* receivedBlock = new int[blockSize];
                MPI_Recv(receivedBlock, blockSize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Copy the received block to the appropriate position in the array
                int startIdx = (i - 1) * participantBlockSize;
                for (int j = 0; j < blockSize; ++j) {
                    decryptedBlocks[startIdx + j] = receivedBlock[j];
                }

                // Clean up the dynamically allocated array
                delete[] receivedBlock;
            }

            // Convert decrypted integers to characters
            std::string decipheredText = "";
            for (size_t i = 0; i < ciphertext.length(); ++i) {
                int decryptedNumber = decryptedBlocks[i];
                char decryptedChar = (decryptedNumber >= 0 && decryptedNumber < 26) ? decryptedNumber + 'A' : '?';

                decipheredText += decryptedChar;
            }

            // Output the deciphered text
            std::cout << "Deciphered Text: " << decipheredText << std::endl;

            // Clean up the dynamically allocated array
            delete[] decryptedBlocks;

            MPI_Finalize();
            return 0;



        }
        else {
            std::cout << "No inverse exists for the given matrix." << std::endl;
        }

       

    }
    else {
        // Participants' code
        std::cout << "I am a participant with rank " << world_rank << std::endl;

        // Receive the broadcasted inverse matrix from the coordinator
        int invA, invB, invC, invD;
        MPI_Bcast(&invA, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&invB, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&invC, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&invD, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Your participant code here
    }

    MPI_Finalize();
    return 0;
}

#ifndef __OUTPUT_H
#define __OUTPUT_H

#include "structure.h"

/**
 * @file output.h
 * @brief Output Header File
 *
 * This file contains function for handling various output-related operations.
 */


/**
 * @brief Create a directory with the given input name in the "../results" folder.
 *
 * @param input The input filename.
 * @return A dynamically allocated string representing the created directory name.
 */
char* createDir(char *);

/**
 * @brief Extract the basename from the given input string.
 *
 * @param in The input string.
 * @return A dynamically allocated string representing the extracted basename.
 */
char* getBasename (char *);

/**
 * @brief Write the contents of a list to the console.
 *
 * @param l The list to be written.
 */
void LST_write(List_t*);

/**
 * @brief Write the information of a molecule to the console.
 *
 * @param m The molecule to be written.
 */
void MOL_write(Molecule_t*);

/**
 * @brief Write the information of a shell to the console.
 *
 * @param s The shell to be written.
 */
void SHL_write(Shell_t*);

/**
 * @brief Write the information of a graph to the console.
 *
 * @param g The graph to be written.
 */
void GPH_write(Graph_t*) ;

/**
 * @brief Write the contents of a molecule to a .mol2 file.
 *
 * @param output The name of the output file.
 * @param m The molecule to be written.
 */
void MOL_writeMol2(char*, Molecule_t*);

/**
 * @brief Write the contents of a shell to a .mol2 file.
 *
 * @param output The name of the output file.
 * @param s The shell to be written.
 */
void SHL_writeMol2(char*, Shell_t*);

/**
 * @brief Write the main output files for a given main structure.
 *
 * @param inputFile The input file name.
 * @param m The main structure to be written.
 */
void writeMainOutput(char*, Main_t*);

/**
 * @brief Write the output files for a given shell.
 *
 * @param inputFile The input file name.
 * @param s The shell to be written.
 * @param tailleMocInit The initial size of the cage.
 */
void writeShellOutput(char* InputFile, Shell_t* s, int tailleMocInit);

#endif

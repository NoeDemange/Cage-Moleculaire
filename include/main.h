#ifndef _MAIN_H
#define _MAIN_H

/**
 * @file main.h
 * @brief Main Header File
 *
 * This file contains the definition of the Options_t structure and function
 * related to program usage, executing R script files, and setting the working directory.
 */

/**
 * @brief Structure to store command-line options.
 *
 * This structure holds command-line options that can be passed to the program. 
 * It includes fields for input filename, alpha value, maximum size, and maximum results.
 */
typedef struct {
  char*         input; /**< Input filename */
  double				alpha; /**< Aphashape parameter value */
	int					sizeMax; /**< Maximum size */
  int      maxResults; /**< Maximum results */
} Options_t;

/**
 * @brief Print the usage of the program.
 *
 * This function prints information about how to use the program, including the available
 * command-line options and their descriptions.
 */
void usage();

/**
 * @brief Load and execute an R script file.
 *
 * This function takes the name of an R script file and executes it in R.
 *
 * @param name Name of the R script file to execute.
 */
void source(const char*);

/**
 * @brief Set the working directory for R scripts.
 *
 * This function sets the working directory for R scripts to the specified directory.
 *
 * @param dir Directory to set as the working directory for R scripts.
 */
void setWorkingDirectory(char*);

#endif
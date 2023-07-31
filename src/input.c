#include "input.h"
#include <string.h>

/**
 * @file input.c
 * @brief This file contains functions that handle the reading of molecule data from input files.
 */

/**************************************/
/* INITIALISATION MOLÉCULE ************/
/**************************************/

/**
* Retrieves the file containing the molecule data.
* Must have an .xyz extension.
*
* @param inputname Name of the file containing the molecule data.
* @param m Address of the molecule.
*/
Molecule_t* readInput_xyz(char* inputname) {
	FILE* filestream = NULL;
	int size, ret;
	Molecule_t* m;
  float Ex = 0; //esperance
  float Ey = 0;
  float Ez = 0;
	
	filestream = fopen(inputname, "r");

  if(!filestream) {
    fprintf(stderr, "The file %s could not be open for reading.\n", inputname);
    exit(EXIT_FAILURE);
  }
	
	ret = fscanf(filestream, "%d", &size);
	m = MOL_create(size);
	
	for (int i = 0; i < size(m); i++) {
		ret = fscanf(filestream, "%s %f %f %f", symbol(atom(m,i)), 
			&atomX(atom(m,i)), &atomY(atom(m,i)),	&atomZ(atom(m,i)));
      Ex += atomX(atom(m,i));
      Ey += atomY(atom(m,i));
      Ez += atomZ(atom(m,i));
	}

  Ex/=size(m);
  Ey/=size(m);
  Ez/=size(m);

  for (int i = 0; i < size(m); i++) {
      atomX(atom(m,i)) -= Ex;
      atomY(atom(m,i)) -= Ey;
      atomZ(atom(m,i)) -= Ez;
	}

	
	fclose(filestream);
	
	if (ret < 0) {
		fprintf(stderr, "An error occured while reading %s.\n", inputname);
		exit(EXIT_FAILURE);
	}

	return m;
}

/**
* Retrieves the covalent radius of atoms.
* They are stored in the rdc.dat file.
*
* @param m Address of the molecule.
*/
void readCovalence(Molecule_t* m) {
  FILE* filestream = NULL;
  int i, j, number, ret;
  Atom_t* atoms;
  
  filestream = fopen("resources/rdc.dat", "r");

  if(!filestream) {
    fprintf(stderr, "The file resources/rdc.dat could not be open for reading.\n");
    exit(EXIT_FAILURE);
  }
  
  ret = fscanf(filestream, "%d", &number);
  atoms = malloc(number*sizeof(Atom_t));

  // Retrieves the covalent radius.
  for (i = 0; i < number; i++) {
    ret = fscanf(filestream, "%s %u", symbol(atoms+i), &radius(atoms+i));
  }
  
  for (i = 0; i < size(m); i++) {
    // Find the corresponding symbol for each atom.
    for (j = 0; strcmp(symbol(atom(m,i)), symbol(atoms + j)) && j < number; j++);

    if (number <= j) {
      fprintf(stderr, "The %s atom is not referenced.\n", symbol(atom(m,i)));
      exit(EXIT_FAILURE);
    }
    else {
      radius(atom(m,i)) = radius(atoms + j);
    }
  }

  fclose(filestream);

  free(atoms);

  if (ret < 0) {
    fprintf(stderr, "An error occured while reading resources/rdc.dat.\n");
    exit(EXIT_FAILURE);
  }
}

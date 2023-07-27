#include <stdlib.h>
#include <time.h>

#include "structure.h"
#include "initialization.h"
#include "voxelization.h"
#include "pathFinding.h"

#define nb_instances 1000000

int main(int argc, char** argv) {
	time_t start = time(NULL);

    //Main_t* m = MN_create();
	//substrat(m) = initMolecule("./demos/substrates/ADENOS10.xyz"); //charge le substrat depuis un fichier .xyz et la centre en (0,0,0)
	time_t startVoxel = time(NULL);
	VOXELGRID voxelGrid = initVoxelGrid();//voxelization(substrat(m));
	time_t endVoxel = time(NULL);
	long seconds = (long) difftime(endVoxel, startVoxel);
	int hours = seconds / 3600;
	seconds -= hours * 3600;
	int minutes = seconds / 60;
	seconds -= minutes * 60;
	printf("\nVoxel time : %d hour(s) %d minute(s) %ld second(s)\n\n", hours, minutes, seconds);
	
	Point3D st = {0,0,0}; Point3D Pend = {GRID_SIZE-1,GRID_SIZE-1, GRID_SIZE-1};

	/*time_t startDij = time(NULL);
	float Ddist = dijkstra(st, Pend, voxelGrid);
	time_t endDij = time(NULL);
	seconds = (long) difftime(endDij, startDij);	
	hours = seconds / 3600;
	seconds -= hours * 3600;
	minutes = seconds / 60;
	seconds -= minutes * 60;
	printf("Dijkstra time : %d hour(s) %d minute(s) %ld second(s)\n", hours, minutes, seconds);
	printf("dijkstra dist : %f\n\n",Ddist);*/

	time_t startA = time(NULL);
	VMap*** vMap = VMap_alloc();
	NodeHeap nodeHeap = NH_initAlloc();
	float Adist = 0.0;
	for(int i = 0; i<nb_instances; i++){
		Adist += aStarPathfinding(st, Pend, voxelGrid, vMap, nodeHeap);
	}
	VMap_free(vMap);
	NH_free(nodeHeap);
	time_t endA = time(NULL);
	seconds = (long) difftime(endA, startA);	
	hours = seconds / 3600;
	seconds -= hours * 3600;
	minutes = seconds / 60;
	seconds -= minutes * 60;
	printf("A* time : %d hour(s) %d minute(s) %ld second(s)\n", hours, minutes, seconds);
	printf("A* dist mean : %f\n",Adist/nb_instances);

	freeVoxelGrid(voxelGrid);


	time_t end = time(NULL);
	seconds = (long) difftime(end, start);	
	hours = seconds / 3600;
	seconds -= hours * 3600;
	minutes = seconds / 60;
	seconds -= minutes * 60;
	printf("\nExecution time : %d hour(s) %d minute(s) %ld second(s)\n", hours, minutes, seconds);

    return 0;
}
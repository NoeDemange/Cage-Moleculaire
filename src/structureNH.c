#include "structure.h"

/**
 * @file structureNH.c
 * @brief NodeHeap and VMap Data Structure Implementation
 *
 * This file contains the implementation of the NodeHeap and VMap data structures
 * and related operations. The NodeHeap represents a priority queue, and VMap is a
 * 3D array of mapping structures used for storing distance g and node indices in the NodeHeap.
 */

/**************************************/
/* NODE HEAP **************************/
/**************************************/

/**
 * @brief Initializes and allocates memory for a new NodeHeap.
 *
 * This function initializes a new NodeHeap structure, sets its size to 0, and allocates memory
 * for the node array based on the GRID_SIZE constant. The function returns the initialized NodeHeap.
 *
 * @return The initialized NodeHeap.
 */
NodeHeap NH_initAlloc(){
    NodeHeap nodeHeap;
    nodeHeap.size = 0;
    nodeHeap.node = malloc(((GRID_SIZE*GRID_SIZE*GRID_SIZE)+1) * sizeof(Node));
    return nodeHeap;
}

/**
 * @brief Frees the memory allocated for a NodeHeap.
 *
 * This function frees the memory allocated for the node array in the NodeHeap.
 *
 * @param nodeHeap The NodeHeap to be freed.
 */
void NH_free(NodeHeap nodeHeap){
    free(nodeHeap.node);
}

/**
 * @brief Swaps two nodes in the NodeHeap and updates VMap indices.
 *
 * This function swaps two nodes at the given indices (i and j) in the NodeHeap.
 * It also updates the corresponding VMap indices after the swap.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param i Index of the first node to be swapped.
 * @param j Index of the second node to be swapped.
 * @param vMap Pointer to the 3D array of VMap structures.
 */
void swap(NodeHeap* nodeHeap, int i, int j, VMap*** vMap){
    Node tmp = nodeHeap->node[i];
    nodeHeap->node[i] = nodeHeap->node[j];
    nodeHeap->node[j] = tmp;
    vMap[nodeHeap->node[i].point.x][nodeHeap->node[i].point.y][nodeHeap->node[i].point.z].indexHeap = i;
    vMap[nodeHeap->node[j].point.x][nodeHeap->node[j].point.y][nodeHeap->node[j].point.z].indexHeap = j;
}

/**
 * @brief Maintains the binary heap property of the NodeHeap.
 *
 * This function is used to maintain the binary heap property of the NodeHeap after
 * a node's priority has been decreased. It recursively adjusts the heap if necessary.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param i Index of the node to be adjusted.
 * @param vMap Pointer to the 3D array of VMap structures.
 */
void stack(NodeHeap* nodeHeap, int i, VMap*** vMap){
    int l, r, min;
    l = 2*i;
    r = (2*i)+1;
    if(l<=nodeHeap->size && (nodeHeap->node[l].f) < (nodeHeap->node[i].f))
        min = l;
    else
        min = i;

    if(r<=nodeHeap->size && (nodeHeap->node[r].f) < (nodeHeap->node[min].f))
        min = r;
    //printf("l : %f; r : %f; i : %f; min : %d\n", nodeHeap->node[l].f, nodeHeap->node[r].f, nodeHeap->node[i].f, min);
    if(min != i){
        swap(nodeHeap,i,min,vMap);
        stack(nodeHeap, min, vMap);
    }
}

/**
 * @brief Inserts a new node into the NodeHeap.
 *
 * This function inserts a new node into the NodeHeap and adjusts the heap to maintain its binary property.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param p The new node to be inserted.
 * @param vMap Pointer to the 3D array of VMap structures.
 */
void NH_insert(NodeHeap* nodeHeap, Node p, VMap*** vMap){
    nodeHeap->size++;
    int i = nodeHeap->size;
    while(i>1 && nodeHeap->node[i/2].f > p.f){
        nodeHeap->node[i] = nodeHeap->node[i/2];
        vMap[nodeHeap->node[i].point.x][nodeHeap->node[i].point.y][nodeHeap->node[i].point.z].indexHeap = i;
        i = i/2;
    }
    nodeHeap->node[i] = p;
    vMap[nodeHeap->node[i].point.x][nodeHeap->node[i].point.y][nodeHeap->node[i].point.z].indexHeap = i;
}

/**
 * @brief Extracts the node with the minimum f value from the NodeHeap.
 *
 * This function extracts the node with the minimum f value (the root) from the NodeHeap,
 * adjusts the heap, and returns the extracted node.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param vMap Pointer to the 3D array of VMap structures.
 * @return The node with the minimum f value.
 */
Node NH_extractMin(NodeHeap* nodeHeap, VMap*** vMap){
    Node min;
    if(nodeHeap->size<1){
        printf("Empty Heap !!");
        exit(EXIT_FAILURE);
    }
    min = nodeHeap->node[1];
    vMap[nodeHeap->node[1].point.x][nodeHeap->node[1].point.y][nodeHeap->node[1].point.z].indexHeap = -1;
    nodeHeap->node[1] = nodeHeap->node[nodeHeap->size];
    vMap[nodeHeap->node[1].point.x][nodeHeap->node[1].point.y][nodeHeap->node[1].point.z].indexHeap = 1;
    nodeHeap->size--;
    stack(nodeHeap, 1, vMap);
    return min;
}

/**
 * @brief increases the priority of a node in the NodeHeap.
 *
 * This function increases the priority (decreases f value) of a node in the NodeHeap
 * and adjusts the heap to maintain its binary property.
 *
 * @param nodeHeap Pointer to the NodeHeap structure.
 * @param vMap Pointer to the 3D array of VMap structures.
 * @param point The 3D coordinates of the node to be updated.
 * @param distG The new g (distance from the start node) value for the node.
 */
void NH_increase_priority(NodeHeap* nodeHeap, VMap*** vMap, Point3D point, float distG){
    int i = vMap[point.x][point.y][point.z].indexHeap;
    nodeHeap->node[i].g = distG;
    nodeHeap->node[i].f = nodeHeap->node[i].g+nodeHeap->node[i].h;
    while(i>1 && nodeHeap->node[i/2].f > nodeHeap->node[i].f){
        swap(nodeHeap, i, i/2, vMap);
    }
}

/**
 * @brief Prints the priorities of nodes in the NodeHeap.
 *
 * This function prints the priorities (f values) of all nodes in the NodeHeap.
 * It is used for debugging and testing purposes.
 *
 * @param nodeHeap The NodeHeap to be printed.
 */
void NH_print(NodeHeap nodeHeap){
    for(int i = 1; i<=nodeHeap.size; i++){
        printf("%f, ",nodeHeap.node[i].f);
    }
    printf("\n");
}

/**
 * @brief Allocates memory for the VMap data structure.
 *
 * This function allocates memory for the VMap 3D array, which is used for storing
 * node indices in the NodeHeap. The function returns the pointer to the allocated VMap.
 *
 * @return Pointer to the allocated VMap 3D array.
 */
VMap*** VMap_alloc(){
    VMap*** vMap = (VMap***)malloc(GRID_SIZE * sizeof(VMap**));
    for (int x = 0; x < GRID_SIZE; x++) {
        vMap[x] = (VMap**)malloc(GRID_SIZE * sizeof(VMap*));
        for (int y = 0; y < GRID_SIZE; y++) {
            vMap[x][y] = (VMap*)malloc(GRID_SIZE * sizeof(VMap));
        }
    }
    return vMap;
}

/**
 * @brief Frees the memory allocated for the VMap data structure.
 *
 * This function frees the memory allocated for the VMap 3D array.
 *
 * @param vMap Pointer to the VMap 3D array to be freed.
 */
void VMap_free(VMap*** vMap){
    // Free dynamically allocated memory
    for (int x = 0; x < GRID_SIZE; x++) {
        for (int y = 0; y < GRID_SIZE; y++) {
            free(vMap[x][y]);
        }
        free(vMap[x]);
    }
    free(vMap);
}
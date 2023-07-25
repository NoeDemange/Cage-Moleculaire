#include "structure.h"

/**************************************/
/* NODE HEAP **************************/
/**************************************/

NodeHeap NH_initAlloc(int size){
    NodeHeap nodeHeap;
    nodeHeap.size = 0;
    nodeHeap.node = malloc((size+1) * sizeof(Node));
    return nodeHeap;
}

void NH_free(NodeHeap nodeHeap){
    free(nodeHeap.node);
}

void swap(NodeHeap* nodeHeap, int i, int j, VMap*** vMap){
    Node tmp = nodeHeap->node[i];
    nodeHeap->node[i] = nodeHeap->node[j];
    nodeHeap->node[j] = tmp;
    vMap[nodeHeap->node[i].point.x][nodeHeap->node[i].point.y][nodeHeap->node[i].point.z].indexHeap = i;
    vMap[nodeHeap->node[j].point.x][nodeHeap->node[j].point.y][nodeHeap->node[j].point.z].indexHeap = j;
}

//garantir le maintien de la propriété de tas binaire.
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

//inserer un nouvel élément dans l’ensemble.
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

//get the min node
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

void NH_decrease_priority(NodeHeap* nodeHeap, VMap*** vMap, Point3D point, float distG){
    int i = vMap[point.x][point.y][point.z].indexHeap;
    nodeHeap->node[i].g = distG;
    nodeHeap->node[i].f = nodeHeap->node[i].g+nodeHeap->node[i].h;
    while(i>1 && nodeHeap->node[i/2].f > nodeHeap->node[i].f){
        swap(nodeHeap, i, i/2, vMap);
    }
}

void NH_print(NodeHeap nodeHeap){
    for(int i = 1; i<=nodeHeap.size; i++){
        printf("%f, ",nodeHeap.node[i].f);
    }
    printf("\n");
}
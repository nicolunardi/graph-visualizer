#include <assert.h>  // assert
#include <math.h>    // fabs
#include <stdlib.h>  // malloc, realloc
#include <string.h>  // strcmp

#include "graph.h"

// size to increment the vertex array by when full
#define V_SIZE_BUFFER 100

// graph ADT
typedef struct Edge {
    struct Vertex *dest;
    size_t weight;
    struct Edge *prev;
    struct Edge *next;
} Edge;

typedef struct Edge *edge;

typedef struct Vertex {
    string data;
    size_t inboundE;
    size_t outboundE;
    double oldRank;
    double pageRank;
    struct List *edges;
    // for dijkstra
    struct Vertex *prev;  // holds the previous vertex towards source
    size_t dist;          // holds distance to source
    bool visited;
} Vertex;

typedef struct Vertex *vertex;

typedef struct Graph_Repr {
    // max size of the adjlist array
    size_t maxSize;
    // keeps track of the next index of adjlist where a vertex can be
    // inserted. Used in order to not go through the entire adjList
    // as maxSize can be a lot larger than the number of Vertices in
    // the list
    size_t hiIndex;
    size_t nV;
    size_t nE;
    vertex *adjList;
    // for dijkstra
    size_t nVisited;  // when nVisited = nV end dijkstra
} Graph_Repr;

// ADT for adjacency list
typedef struct List {
    size_t size;
    struct Edge *head;
    struct Edge *tail;
} List;

typedef struct List *list;

// ***list ADT functions***

static list list_create(void);
static void list_destroy(list, bool);
static bool list_is_empty(list l);
// adds an edge to the list if not already in the list. returns
// true if edge added else false
static bool list_add(graph G, list l, vertex dest, size_t weight);
static size_t list_remove(graph G, list l, vertex target);
static bool list_contains(list l, vertex target);
// updateInbound is a bool that indicates if the inbound value
// of the edges' destination vertex is to be updated. Should be
// false when the function is called during graph destruction as the
// dest vertex may have already been freed and updating inbound value
// is unnecessary
static void freeNodes(edge head, bool updateInbound);
static edge createEdgeNode(vertex dest, size_t weight);
static void list_append(list l, vertex dest, size_t weight);
static size_t list_shrink(list l);
static size_t list_shift(list l);

// ***helper functions***

static vertex createVertex(string data);
// returns a vertex if found else NULL
static vertex getVertex(graph G, string target);
// returns an edge from a list if it exists, else NULL
static edge getEdge(list l, string target);
// increases the size of the adjList by the amount determined by
// V_SIZE_BUFFER
static void increaseAdjListSize(graph G);
// print the vertices of a graph
static void printVertices(graph G, FILE *file);
// print the edges of a graph
static void printEdges(list l, vertex src, FILE *file);
// removes all edges inbound on a vertex. Used when a vertex is removed
static void removeInboundEdges(graph G, vertex target);

// ***pagerank helper functions***

// initialize the pagerank of all the vertices
static void initPageRank(graph G);
// update the old rank of all vertices
static void updateOldRank(graph G);
// calculate the sinkRank for all V with no outbound edges
static double getSinkRank(graph G, double damping);
// updates the pagerank for each V and updates the highest difference
// between pagerank and oldrank
static void updatePageRank(
    graph G, double damping, double sink, double *hiDiff);
// used for qsort in graph_viewrank
static int qSortCompFunc(const void *V1, const void *V2);

// ***Dijkstra helper functions***

// initialize all the vertex dist and prev. sets source dist to 0
static void dijkstraInitValues(graph G, vertex source);
// finds the V with the min distance from source not already visited
static vertex getMinDistV(graph G);
// updates the distances from source of all adjacent vertices
static void updateNeighbourDist(vertex src);
// prints the path from source to dest recursively
static void printPath(vertex V);

graph graph_create(void) {
    graph newG = malloc(sizeof *newG);
    if (!newG) {
        // return null if malloc failed
        fprintf(stderr, "%s\n", MEMORY_ALLOC_ERROR);
        return newG;
    }
    newG->hiIndex = 0;
    newG->nV = 0;
    newG->nE = 0;
    newG->adjList = malloc(
        V_SIZE_BUFFER * sizeof(*newG->adjList));
    // set the max size of the graph to the buffer
    newG->maxSize = V_SIZE_BUFFER;
    // initialize the adj list to null
    for (size_t i = 0; i < newG->maxSize; i++) {
        newG->adjList[i] = NULL;
    }
    newG->nVisited = 0;
    return newG;
}

void graph_destroy(graph G) {
    assert(G);
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        // if the vertex exists destroy its adjlist
        if (V) {
            list_destroy(V->edges, false);
            free(V->data);
            free(V);
        }
    }
    free(G->adjList);
    free(G);
    return;
}

void graph_show(graph G, FILE *file) {
    if (!file) file = stdout;
    printVertices(G, file);
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) {
            printEdges(V->edges, V, file);
        }
    }
}

void graph_add_vertex(graph G, string data) {
    // if the vertex already exists do nothing and return
    if (graph_has_vertex(G, data)) return;
    // if the vertex array is full, increase its size
    if (G->maxSize == G->hiIndex) {
        increaseAdjListSize(G);
    }
    vertex newV = createVertex(data);
    G->adjList[G->hiIndex] = newV;
    G->hiIndex++;
    G->nV++;
    return;
}

bool graph_has_vertex(graph G, string data) {
    if (!G) return false;
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        // ensure vertex at index i not empty and that data matches
        if (V && !strcmp(V->data, data)) {
            return true;
        }
    }
    return false;
}

void graph_remove_vertex(graph G, string data) {
    assert(G);
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V && !strcmp(V->data, data)) {
            removeInboundEdges(G, V);
            list_destroy(V->edges, true);
            free(V->data);
            free(V);
            G->adjList[i] = NULL;
            G->nV--;
        }
    }
    return;
}

size_t graph_vertices_count(graph G) {
    if (!G) return 0;
    return G->nV;
}

void graph_add_edge(graph G, string src, string dest, size_t weight) {
    assert(G);
    // ensure both vertices exist if not create them
    if (!graph_has_vertex(G, src)) graph_add_vertex(G, src);
    if (!graph_has_vertex(G, dest)) graph_add_vertex(G, dest);
    vertex srcV = getVertex(G, src);
    vertex destV = getVertex(G, dest);
    // list_add will only add an edge if it doesnt already exist
    bool edgeAdded = list_add(G, srcV->edges, destV, weight);
    if (edgeAdded) srcV->outboundE = srcV->edges->size;
    return;
}

bool graph_has_edge(graph G, string src, string dest) {
    if (!G || !graph_has_vertex(G, src) || !graph_has_vertex(G, dest)) {
        return false;
    }
    vertex srcV = getVertex(G, src);
    vertex destV = getVertex(G, dest);
    // if either of src or dest vertices dont exit return false
    if (!srcV || !destV) return false;
    return list_contains(srcV->edges, destV);
}

size_t graph_remove_edge(graph G, string src, string dest) {
    if (!G || !graph_has_vertex(G, src) || !graph_has_vertex(G, dest)) {
        return 0;
    }
    size_t weight = 0;
    vertex srcV = getVertex(G, src);
    vertex destV = getVertex(G, dest);
    if (srcV && destV) {
        weight = list_remove(G, srcV->edges, destV);
        srcV->outboundE = srcV->edges->size;
    }
    return weight;
}

void graph_set_edge(graph G, string src, string dest, size_t weight) {
    if (!G) return;
    vertex srcV = getVertex(G, src);
    if (srcV) {
        edge targetE = getEdge(srcV->edges, dest);
        if (targetE) targetE->weight = weight;
    }
    return;
}

size_t graph_get_edge(graph G, string src, string dest) {
    if (!G) return 0;
    size_t weight = 0;
    vertex srcV = getVertex(G, src);
    if (srcV) {
        edge targetE = getEdge(srcV->edges, dest);
        if (targetE) weight = targetE->weight;
    }
    return weight;
}

size_t graph_edges_count(graph G, string target) {
    if (!G) return 0;
    size_t outboundE = 0;
    vertex targetV = getVertex(G, target);
    if (targetV) outboundE = targetV->outboundE;
    return outboundE;
}

void graph_pagerank(graph G, double damping, double delta) {
    // the number of V in the graph
    size_t N = G->nV;
    // initiates all the pagerank of every V to 1/N, oldrank is already
    // set to 0 during V creation
    initPageRank(G);
    // hiDiff holds the highest difference between pagerank and old
    // rank of any V. It is updated within updatePageRank
    double hiDiff = 1.0 / N;
    while (hiDiff > delta) {
        // set to 0 in order to find highest difference, will be
        // apparent within updatePageRank
        hiDiff = 0.0;
        // set all the oldranks of each V to their pagerank
        updateOldRank(G);
        // calculates the sink rank from all the V with no outbound
        // edges
        double sinkRank = getSinkRank(G, damping);
        // the rest of the equation as per the pseudocode, passing in
        // the adress of hiDiff to modifiy it within
        updatePageRank(G, damping, sinkRank, &hiDiff);
    }
}

void graph_viewrank(graph G, FILE *file) {
    // sort the adjList by pagerank and alpha
    qsort(G->adjList, G->hiIndex, sizeof(vertex), qSortCompFunc);
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) fprintf(file, "%s (%.3f)\n", V->data, V->pageRank);
    }
    return;
}

void graph_shortest_path(graph G, string source) {
    vertex srcV = getVertex(G, source);
    // init all the values to their starting values
    dijkstraInitValues(G, srcV);
    while (G->nVisited < G->nV) {
        // get the V with the shortest dist to source that hasnt yet
        // been visited
        vertex minV = getMinDistV(G);
        // update the shortest distance for all neighbours of minV
        updateNeighbourDist(minV);
        minV->visited = true;
        G->nVisited++;
    }
    return;
}

void graph_view_path(graph G, string destination) {
    vertex destV = getVertex(G, destination);
    printPath(destV);
    printf("\n");
}

// ***list ADT functions***

static list list_create(void) {
    list newList = malloc(sizeof *newList);
    if (!newList) {
        fprintf(stderr, "%s\n", MEMORY_ALLOC_ERROR);
        return newList;
    }
    newList->size = 0;
    newList->head = newList->tail = NULL;
    return newList;
}

static void list_destroy(list l, bool updateInbound) {
    assert(l);
    if (l->head) {
        freeNodes(l->head, updateInbound);
    }
    free(l);
    return;
}

static bool list_is_empty(list l) {
    assert(l);
    return l->size == 0;
}

static bool list_add(graph G, list l, vertex dest, size_t weight) {
    bool edgeAdded = false;
    if (!list_contains(l, dest)) {
        list_append(l, dest, weight);
        dest->inboundE++;
        edgeAdded = true;
        G->nE++;
    }
    return edgeAdded;
}

static size_t list_remove(graph G, list l, vertex target) {
    if (list_is_empty(l)) return 0;
    size_t weight = 0;
    if (!strcmp(l->head->dest->data, target->data)) {
        // target node is the head of list
        weight = list_shift(l);
        G->nE--;
    } else if (!strcmp(l->tail->dest->data, target->data)) {
        // target is the tail of the list
        weight = list_shrink(l);
        G->nE--;
    } else {
        edge temp = l->head;
        while (temp && (strcmp(temp->dest->data, target->data))) {
            temp = temp->next;
        }
        if (!temp) return 0;  // gone through whole list and not found

        temp->prev->next = temp->next;
        temp->next->prev = temp->prev;
        weight = temp->weight;
        free(temp);
        target->inboundE--;
        l->size--;
        G->nE--;
    }
    return weight;
}

static bool list_contains(list l, vertex target) {
    assert(l);
    edge temp = l->head;
    while (temp) {
        // found match
        if (!strcmp(temp->dest->data, target->data)) return true;
        temp = temp->next;
    }
    return false;
}

static void freeNodes(edge head, bool updateInbound) {
    if (head->next) freeNodes(head->next, updateInbound);
    // when an edge is removed, lower the inbound counter for the
    // edge's destination vertex.
    if (updateInbound) head->dest->inboundE--;
    free(head);
}
static edge createEdgeNode(vertex dest, size_t weight) {
    edge newNode = malloc(sizeof *newNode);
    if (!newNode) {
        fprintf(stderr, "%s\n", MEMORY_ALLOC_ERROR);
        exit(0);
    }
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = newNode->prev = NULL;
    return newNode;
}

static void list_append(list l, vertex dest, size_t weight) {
    edge newNode = createEdgeNode(dest, weight);
    if (list_is_empty(l)) {
        l->head = l->tail = newNode;
    } else {
        l->tail->next = newNode;
        newNode->prev = l->tail;
        l->tail = newNode;
    }
    l->size++;
    return;
}

static size_t list_shrink(list l) {
    size_t weight = l->tail->weight;
    l->tail->dest->inboundE--;
    if (l->head == l->tail) {
        // single node case
        free(l->tail);
        l->head = l->tail = NULL;
    } else {
        l->tail = l->tail->prev;
        free(l->tail->next);
        l->tail->next = NULL;
    }
    l->size--;
    return weight;
}

static size_t list_shift(list l) {
    size_t weight = l->head->weight;
    l->head->dest->inboundE--;
    if (l->head == l->tail) {
        // if theres only one node
        free(l->head);
        l->head = l->tail = NULL;
    } else {
        l->head = l->head->next;  // shift head across
        free(l->head->prev);
        l->head->prev = NULL;
    }
    l->size--;
    return weight;
}

// ***helper functions***

static vertex createVertex(string data) {
    vertex newV = malloc(sizeof *newV);
    if (!newV) {
        fprintf(stderr, "%s\n", MEMORY_ALLOC_ERROR);
        exit(0);
    }
    newV->data = strdup(data);
    newV->inboundE = 0;
    newV->outboundE = 0;
    newV->pageRank = 0.0;
    newV->oldRank = 0.0;
    newV->edges = list_create();
    newV->prev = NULL;
    newV->dist = INF;
    newV->visited = false;
    return newV;
}

static vertex getVertex(graph G, string target) {
    vertex targetV = NULL;
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V && !strcmp(V->data, target)) {
            targetV = V;
            break;
        }
    }
    return targetV;
}

static edge getEdge(list l, string target) {
    edge targetE = l->head;
    while (targetE && (strcmp(targetE->dest->data, target))) {
        targetE = targetE->next;
    }
    return targetE;
}

static void increaseAdjListSize(graph G) {
    size_t newSize = V_SIZE_BUFFER + G->maxSize;
    vertex *temp = realloc(G->adjList, sizeof(*temp) * newSize);
    if (!temp) {
        fprintf(stderr, "%s\n", MEMORY_ALLOC_ERROR);
        exit(0);
    }
    G->adjList = temp;
    G->maxSize += V_SIZE_BUFFER;
    // set new space to NULL
    for (size_t i = G->hiIndex; i < G->maxSize; i++) {
        G->adjList[i] = NULL;
    }
    return;
}

static void printVertices(graph G, FILE *file) {
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) {
            fprintf(file, "%s\n", V->data);
        }
    }
    return;
}

static void printEdges(list l, vertex src, FILE *file) {
    edge temp = l->head;
    while (temp) {
        fprintf(
            file,
            "%s %s %zu\n",
            src->data, temp->dest->data, temp->weight);
        temp = temp->next;
    }
    return;
}

static void removeInboundEdges(graph G, vertex target) {
    // traverse through all edges and rmeove all edges that have
    // target as a dest
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        // if the vertex exists and has at least 1 outbound edge
        if (V && V->outboundE) {
            list_remove(G, V->edges, target);
            V->outboundE = V->edges->size;
        }
    }
    return;
}

static void initPageRank(graph G) {
    // set the pagerank of all V to 1/N, oldrank already set to 0
    // during vertex creation
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) {
            V->pageRank = 1.0 / G->nV;
        }
    }
    return;
}

static void updateOldRank(graph G) {
    // set the oldrank of each V to the pagerank
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) {
            V->oldRank = V->pageRank;
        }
    }
    return;
}

static double getSinkRank(graph G, double damping) {
    double sinkRank = 0.0;
    // for all V in G that have no outbound edges
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V && !V->outboundE) {
            sinkRank += (damping * (V->oldRank / G->nV));
        }
    }
    return sinkRank;
}

static void updatePageRank(
    graph G, double damping, double sink, double *hiDiff) {
    // for all the vertices V in G
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) {
            V->pageRank = (sink + ((1.0 - damping) / G->nV));
            // for all the vertices I in G that have and outbound edge
            // to V
            for (size_t j = 0; j < G->hiIndex; j++) {
                vertex I = G->adjList[j];
                if (I && (getEdge(I->edges, V->data))) {
                    V->pageRank +=
                        ((damping * I->oldRank) / I->outboundE);
                }
            }
            // update hiDiff
            double diff = fabs(V->pageRank - V->oldRank);
            if (diff > *hiDiff) *hiDiff = diff;
        }
    }
    return;
}

static int qSortCompFunc(const void *V1, const void *V2) {
    // cast the void pointer to a pointer to vertex for manipulation
    vertex *myV1 = (vertex *)V1;
    vertex *myV2 = (vertex *)V2;
    if (!(*myV1) || !(*myV2)) return 0;
    int compValue;
    // compare the pagerank first and then alphabetically if the
    // ranks are the same. Used compValue to store the strcmp value
    // in order to not call strcmp twice
    if ((*myV1)->pageRank > (*myV2)->pageRank) {
        return -1;
    } else if ((*myV1)->pageRank < (*myV2)->pageRank) {
        return 1;
    } else if ((compValue = strcmp((*myV1)->data, (*myV2)->data)) < 0) {
        return -1;
    } else if (compValue > 0) {
        return 1;
    } else return 0;
}

static void dijkstraInitValues(graph G, vertex source) {
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V) {
            V->prev = NULL;
            V->dist = INF;
            V->visited = false;
        }
    }
    source->dist = 0;
}

static vertex getMinDistV(graph G) {
    size_t min = INF;
    vertex minV = NULL;
    for (size_t i = 0; i < G->hiIndex; i++) {
        vertex V = G->adjList[i];
        if (V && !V->visited) {
            if (V->dist <= min) {
                min = V->dist;
                minV = V;
            }
        }
    }
    return minV;
}

static void updateNeighbourDist(vertex src) {
    // src has no neighbours
    if (!src->outboundE) return;
    size_t newDist = src->dist + 1;
    edge temp = src->edges->head;
    while (temp) {
        if (newDist < temp->dest->dist) {
            temp->dest->dist = newDist;
            temp->dest->prev = src;
        }
        temp = temp->next;
    }
    return;
}

static void printPath(vertex V) {
    if (V && !V->prev) {
        printf("%s", V->data);
    } else if (V && V->prev) {
        printPath(V->prev);
        printf(" -> %s", V->data);
    }
    return;
}

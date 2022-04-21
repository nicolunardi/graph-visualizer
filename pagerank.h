#ifndef PAGERANK_H
#define PAGERANK_H

#include "graph.h"

void graph_pagerank(graph G, double damping, double delta);
void graph_viewrank(graph G, FILE *file);

#endif

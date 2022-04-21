
<h1 align="center">Graph Visualizer</h1>

<!-- Status -->

<!-- <h4 align="center"> 
	🚧  Graph-visualizer 🚀 Under construction...  🚧
</h4> 

<hr> -->

<p align="center">
  <a href="#dart-About">About</a> &#xa0; | &#xa0;
  <a href="#sparkles-features">Features</a> &#xa0; | &#xa0;
  <a href="#rocket-technologies">Technologies</a> &#xa0; | &#xa0;
  <a href="#white_check_mark-requirements">Requirements</a> &#xa0; | &#xa0;
  <a href="#checkered_flag-starting-out">Starting Out</a> &#xa0; | &#xa0;
  <a href="#memo-license">License</a> &#xa0; | &#xa0;
  <a href="https://github.com/nicolunardi" target="_blank">Author</a>
</p>

<br>

## :dart: About ##

Graph Visualizer is a library written entirely in C that allows for the creating of graphs using adjacency lists. Once the graphs are created, you can add/remove/edit vertices (weighted or not) as well as edges. The adjacency list representation of the graph can then be printed out to a specified file.  

The library also has functions to perform a simple Pagerank algorithm on the graph as well as calculate and display the shortest path using Dijkstra's algorithm

## :sparkles: Features ##

[:white_check_mark:] Create graphs using adjacency lists;\
[:white_check_mark:] Add/remove/edit Vertices;\
[:white_check_mark:] Add/edit/remove edges;\
[:white_check_mark:] Calculate **pagerank**;\
[:white_check_mark:] Calculate shortest path with **Dijkstra's algorithm**;\
[:x:] GUI;

## :rocket: Technologies ##

The library is entirely written in C

## :white_check_mark: Requirements ##

none

## :checkered_flag: Starting Out ##

So you want to start creating graphs?

```bash
# Clone the library
$ git clone https://github.com/nicolunardi/graph-visualizer

# Access it
$ cd graph-visualizer
```

First of all, let's create an empty graph.

```C
// import the necessary files
#include "graph.h"

int main(int argc, char const *argv[]) {
    // create the graph using graph_create from graph.h
    // This will allocate the required memory for a new graph and return
    // return a pointer to the new graph. 
    graph G = graph_create();

    ...  
}
```

Now lets add some vertices to the graph.

```C
// import the necessary files
#include "graph.h"

int main(int argc, char const *argv[]) {
    ...

    // We add some vertices using the graph_add_vertex function from graph.h.

    graph_add_vertex(G, "0");
    graph_add_vertex(G, "1");
    graph_add_vertex(G, "2");
    graph_add_vertex(G, "3");
    graph_add_vertex(G, "4");
    graph_add_vertex(G, "5");
    graph_add_vertex(G, "6");
    graph_add_vertex(G, "7");

    ...  
}
```

Now lets add some edges.

```C
// import the necessary files
#include "graph.h"

int main(int argc, char const *argv[]) {
    ...

    // We add some edges using the graph_add_edge function from graph.h.
    // The function takes 4 parameters; the graph, the first and second vertices,
    // and the edge weight.

    graph_add_edge(G, "0", "1", 1);
    graph_add_edge(G, "1", "2", 1);
    graph_add_edge(G, "2", "3", 1);
    graph_add_edge(G, "2", "6", 1);
    graph_add_edge(G, "2", "7", 1);
    graph_add_edge(G, "3", "4", 1);
    graph_add_edge(G, "7", "5", 1);
    graph_add_edge(G, "6", "5", 1);
    graph_add_edge(G, "4", "5", 1);
    graph_add_edge(G, "5", "4", 1);

    ...  
}
```

Then we can generate a shortest path table for a source vertex and display shortest
path to a destination vertex. We must import "dijkstra.h"

```C
// import the necessary files
#include "graph.h"
#include "dijkstra.h" // import dijkstra.h to be able to use shortest path functions

int main(int argc, char const *argv[]) {
    ...

    // Use graph_shortest_path from dijkstra.h to generate the shortest path table
    // from source vertex ("0" in this case).

    graph_shortest_path(G, "0");

    // Display the shortest path from vertex "0" to vertex "4"
    graph_view_path(G, "4"); // 0 -> 1 -> 2 -> 3 -> 4

    ...
}
```

Once done with a graph, you must de-allocate any memory used.

```C
int main(int argc, char const *argv[]) {
    ...

    graph_destroy(G);

    return 0
}

```

Other functions from "graph.h"

```C

// Will print to a provided file the vertices and their associated edges as 
// well as weights
void graph_show (graph G, FILE *file);

// Check to see if a graph has a particular vertex
bool graph_has_vertex (graph G, string vertex);

// remove a vertex from a graph
void graph_remove_vertex (graph G, string vertex);

// Show the number of vertices in a graph
size_t graph_vertices_count (graph G);

// Check if a graph has a particular edge
bool graph_has_edge (graph G, string vertex1, string vertex2);

// Remove a particular edge
size_t graph_remove_edge (graph G, string vertex1, string vertex2);

// Change the weight of a particular edge
void graph_set_edge (graph G, string vertex1, string vertex2, size_t weight);

// Get the weight of a particular edge
size_t graph_get_edge (graph G, string vertex1, string vertex2);

// Display the degree of a vertex
size_t graph_edges_count (graph G, string vertex);


```

Functions from "pagerank.h"

```C
// Generate the pagerank for a graph
void graph_pagerank(graph G, double damping, double delta);

// Display all vertices in order of pagerank
void graph_viewrank(graph G, FILE *file);
```

## :memo: License ##

This project is under license from MIT

Made with :heart: by <a href="https://github.com/nicolunardi" target="_blank">Nico Lunardi</a>

&#xa0;

<a href="#top">Back to top</a>

#ifndef CONN_COMP_IGRAPH_H
#define CONN_COMP_IGRAPH_H

#include <igraph/igraph.h>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

// Function to compute connected components of an igraph
int conn_comp_igraph(vector<vector<int> > dataset);

bool createGraphFromDataset(const vector<vector<int> >& dataset, igraph_t* graph);

#endif // CONN_COMP_IGRAPH_H

#include "conn_comp_igraph.h"

int conn_comp_igraph(vector<vector<int> > dataset) {

    // Create a graph from vector dataset
    igraph_t graph;
    if (!createGraphFromDataset(dataset, &graph)) {
        cerr << "Failed to create graph from dataset." << endl;
        return 1;
    }

    // Print the number of vertices and edges in the graph
    cout << "Number of vertices: " << igraph_vcount(&graph) << endl;
    cout << "Number of edges: " << igraph_ecount(&graph) << endl;

    // conn comp
    igraph_vector_int_t result;
    igraph_vector_int_t membership; // Vector to store membership
    igraph_integer_t num_components;

    // Initialize igraph objects
    igraph_vector_int_init(&result, 0);
    igraph_vector_int_init(&membership, 0);

    // Find connected components
    igraph_connected_components(&graph, &result, NULL, &num_components, IGRAPH_WEAK);

    // Print the connected components
    cout << "Connected components:" << endl;
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&result); i++) {
        igraph_integer_t component_id = VECTOR(result)[i];
        cout << "Vertex " << (i + 1) << " belongs to component " << (component_id + 1) << endl;
    }

    // Number of components
    cout << "Number of connected components: " << num_components << endl;

    // Clean up
    igraph_vector_int_destroy(&result);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&graph);

    return 0;
}

// Function to create igraph from vector of vectors of edges
bool createGraphFromDataset(const vector<vector<int> >& dataset, igraph_t* graph) {
    // Initialize the graph
    igraph_empty(graph, 0, false);

    // Determine the number of vertices in the graph
    size_t numVertices = dataset.size();

    // Add vertices to the graph
    igraph_vector_t vertices;
    igraph_vector_init_range(&vertices, 0, numVertices - 1); // Initialize with sequential vertex IDs
    igraph_add_vertices(graph, numVertices, NULL);

    // Add edges from the dataset
    for (size_t i = 0; i < numVertices; ++i) {
        for (int neighbor : dataset[i]) {
            // Check if the neighbor ID is within the range of vertices
            if (neighbor >= 0 && neighbor < numVertices) {
                // Add the edge to the graph
                igraph_add_edge(graph, i, neighbor);
            } else {
                cerr << "Error: Invalid vertex ID in dataset." << endl;
                igraph_vector_destroy(&vertices);
                return false;
            }
        }
    }

    igraph_vector_destroy(&vertices);

    return true;
}
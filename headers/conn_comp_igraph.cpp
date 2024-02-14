#include "conn_comp_igraph.h"

/* function to compute connected components of an igraph

*/
vector<vector<int> > conn_comp_igraph(vector<vector<int> > dataset, int n, igraph_integer_t &num_components) {
    vector<vector<int> > conncomps;
    conncomps.resize(dataset.size());

    // adact graph for igraph
	for (int i = 0; i < dataset.size(); i++) {

		// delete row < minpts
		if (dataset[i][0] == 0){
			dataset.erase(dataset.begin() + i); // Erase row
			i--;
		} else {
			dataset[i].erase(dataset[i].begin()); // Erase the first item (is Minpts?)
		}
	}

	// cout << "dataset post -----------" << endl;
	// for (const auto& row : dataset) {
	// 	cout << "- " ;
	// 	for (int elem : row) {
	// 		cout << elem << " ";
	// 	}
	// 	cout << endl;
	// }
    // cout << endl;

    // Create a graph from vector dataset
    igraph_t graph;
    if (!createGraphFromDataset(dataset, &graph, n)) {
        cerr << "Failed to create graph from dataset." << endl;
        return conncomps;
    }

    // Print the number of vertices and edges in the graph
    cout << "Number of vertices: " << igraph_vcount(&graph) << endl;
    cout << "Number of edges: " << igraph_ecount(&graph) << endl;

    // conn comp
    igraph_vector_int_t result;
    igraph_vector_int_t membership; // Vector to store membership

    // Initialize igraph objects
    igraph_vector_int_init(&result, 0);
    igraph_vector_int_init(&membership, 0);

    // Find connected components
    igraph_connected_components(&graph, &result, NULL, &num_components, IGRAPH_STRONG);

    // extracting the connected components
    for (int i = 0; i < igraph_vcount(&graph); ++i) {
        // add the component it belongs to
        conncomps[i].push_back(VECTOR(result)[i]);

        // add the items of the vertex
        for (int j = 0; j < dataset[i].size(); ++j) {
            conncomps[i].push_back(dataset[i][j]);
        }
    }

    // cout << "conncomp:" << endl;
	// for (const auto& row : conncomps) {
    //     cout << "- ";
	// 	for (int elem : row) {
	// 		cout << elem << " ";
	// 	}
	// 	cout << endl << endl;
	// }

    // Number of components
    cout << "Number of connected components: " << num_components << endl;

    // Clean up
    igraph_vector_int_destroy(&result);
    igraph_vector_int_destroy(&membership);
    igraph_destroy(&graph);

    return conncomps;
}

/* function to adapt the graph in a proper data structure

*/
bool createGraphFromDataset(const vector<vector<int> >& dataset, igraph_t* graph, int n) {
    // Initialize the graph
    igraph_empty(graph, 0, false);

    // Determine the number of vertices in the graph
    int numVertices = dataset.size();
    

    // Initialization failed
    if (numVertices == 0){
        cerr << "Error: MinPts did not take elements, please decrease the value or increate the upper bound of epsilon." << endl;
        return false;
    }

    // Add vertices to the graph
    igraph_vector_t vertices;
    igraph_vector_init_range(&vertices, 0, numVertices - 1); // Initialize with sequential vertex IDs
    igraph_add_vertices(graph, n, NULL);

    // Add edges from the dataset
    for (size_t i = 0; i < numVertices; i++) {
        for (int neighbor : dataset[i]) {
            // cout << "i " << i << " neighbor " << neighbor << " n " << n << endl;

            // cout << "dataset[i]:" << endl;
            // for (int i : dataset[i]) {
            //     cout << i << " ";
            // }
            // cout << endl;

            // cout << "Number of vertices: " << igraph_vcount(graph) << endl;
            // cout << "Number of edges: " << igraph_ecount(graph) << endl;

            // Check if the neighbor ID is within the range of vertices
            if (neighbor >= 0 && neighbor <= n) {
                // cout << "0 <= {" << neighbor << "} <= " << n << endl;
                // Add the edge to the graph
                // cout << "pre" << endl;
                // cout << "inserting i: " << i << "; neighbor: " << neighbor << "; numVertices: " << numVertices << endl;
                igraph_add_edge(graph, i, neighbor);
                // cout << "post" << endl;
            } else {
                cerr << "Error: Invalid vertex ID in dataset." << endl;
                igraph_vector_destroy(&vertices);
                return false;
            }
        }
    }

    cout << "numVertices: " << numVertices << endl;
    return true;
}
#include "GWU.h"

int main()
{
    std::vector<InputGraph> input_graphs;
    
    std::ifstream file("test.txt");
    if(!file.is_open()) 
    {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    int mode = -1;

    std::string line;
    while (std::getline(file, line)) {
        if (line[0] == 't') {
            mode += 1;
            input_graphs.push_back(InputGraph());
            continue;
        } else if (mode == -1) {
            int vertex1, vertex2, edge;
            double utility;
            std::istringstream iss(line);
            iss >> vertex1 >> vertex2 >> edge >> utility;
            external_utility[std::make_tuple(vertex1, vertex2, edge)] = utility;
        } else if (line[0] == 'v') {
            char _;
            int vertex, label;
            std::istringstream iss(line);
            iss >> _ >> vertex >> label;
            input_graphs[mode].vertices.push_back(std::make_pair(vertex, label));
            input_graphs[mode].dfs_order[vertex] = -1;
        } else if (line[0] == 'e') {
            char _;
            int vertex1, vertex2, label;
            double utility;
            std::istringstream iss(line);
            iss >> _ >> vertex1 >> vertex2 >> label >> utility;

            if (input_graphs[mode].adjacency.find(vertex1) != input_graphs[mode].adjacency.end()) {
                input_graphs[mode].adjacency[vertex1].push_back(vertex2);
            } else {
                input_graphs[mode].adjacency[vertex1] = {vertex2};
            }

            if (input_graphs[mode].adjacency.find(vertex2) != input_graphs[mode].adjacency.end()) {
                input_graphs[mode].adjacency[vertex2].push_back(vertex1);
            } else {
                input_graphs[mode].adjacency[vertex2] = {vertex1};
            }

            input_graphs[mode].edges.push_back(std::make_tuple(vertex1, vertex2, label, utility));
        }
    }

    file.close();

    // Print the contents of input_graphs and external_utility
    // for (const auto& graph : input_graphs) {
    //     // Assuming you have a print method in the InputGraph class
    //     std::cout<<"-----------------"<<std::endl;
    //     graph.print();
    // }

    // Print the contents of external_utility
    // for (const auto& entry : external_utility) {
    //     const auto& key = entry.first;
    //     double value = entry.second;
    //     std::cout << "(" << std::get<0>(key) << ", " << std::get<1>(key) << ", " << std::get<2>(key) << "): " << value << std::endl;
    // }

    std::vector<Graph> graphs;

    for (auto& inputGraph : input_graphs) {
        for (const auto& v : inputGraph.vertices) {
            if (inputGraph.dfs_order[v.first] == -1) 
            {
                inputGraph.dfs_visit(v.first);
            }
        }

        Graph graph;
        for (const auto& v : inputGraph.vertices) {
            graph.vertices[inputGraph.dfs_order[v.first]] = v.second;
        }

        for (const auto& e : inputGraph.edges) {
            int v0_dfs = inputGraph.dfs_order.at(std::get<0>(e));
            int v1_dfs = inputGraph.dfs_order.at(std::get<1>(e));

            graph.adjacency[v0_dfs].push_back(std::make_tuple(v1_dfs, std::get<2>(e), std::get<3>(e)));
            graph.adjacency[v1_dfs].push_back(std::make_tuple(v0_dfs, std::get<2>(e), std::get<3>(e)));
        }

        graphs.push_back(graph);
    }

    // for (const auto& graph : graphs) {
    //     // Assuming you have a print method in the InputGraph class
    //     std::cout<<"-----------------"<<std::endl;
    //     graph.print();
    // }
    double min_util = 0.0;

    for (auto& graph : graphs) {
        min_util += graph.graph_utility(external_utility);
        std::cout << graph.graph_utility(external_utility) << std::endl;
    }

    min_util *= 1.0 / 3.0;

    // std::cout << "----------------------------------" << std::endl;
    // std::cout << "Min Utility: " << min_util << std::endl;

    GSpan({}, graphs, min_util, 0);

    // system("Pause");
    return 0;
}
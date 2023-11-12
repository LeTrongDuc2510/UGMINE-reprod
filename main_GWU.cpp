#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <fstream>
#include <sstream>

// Example usage:
// g++ -std=c++11 mainG.cpp -o mainG
// ./mainG

// Hash function for std::tuple
namespace std {
    template <>
    struct hash<std::tuple<std::string, std::string, int>> {
        size_t operator()(const std::tuple<std::string, std::string, int>& key) const {
            auto hash1 = std::hash<std::string>{}(std::get<0>(key));
            auto hash2 = std::hash<std::string>{}(std::get<1>(key));
            auto hash3 = std::hash<int>{}(std::get<2>(key));

            // Combine the hash values
            size_t seed = hash1;
            seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            return seed;
        }
    };

    template <>
    struct hash<std::tuple<int, int, int>> {
        size_t operator()(const std::tuple<int, int, int>& key) const {
            size_t hash_val = 0;
            hash_val ^= std::hash<int>{}(std::get<0>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<int>{}(std::get<1>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<int>{}(std::get<2>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            return hash_val;
        }
    };

    template <>
    struct equal_to<std::tuple<int, int, int>> {
        bool operator()(const std::tuple<int, int, int>& lhs, const std::tuple<int, int, int>& rhs) const {
            return lhs == rhs;
        }
    };
}

class InputGraph {
public:
    std::vector<std::pair<int, int>> vertices;
    std::vector<std::tuple<int, int, int, int>> edges;
    std::unordered_map<int, std::vector<int>> adjacency;
    std::unordered_map<int, int> dfs_order;
    int time;

    InputGraph() : time(0) {}

    void dfs_visit(int u) 
    {
        dfs_order[u] = time;
        time++;

        if (adjacency.find(u) == adjacency.end()) {
            return;
        }

        for (int v : adjacency[u]) {
            if (dfs_order[v] == -1) {
                dfs_visit(v);
            }
        }
    }

    void print() const {
        std::cout << "Vertices: ";
        for (std::pair<int, int> v : vertices) {
            std::cout << "(" << v.first << ", " << v.second << ") ";
        }

        std::cout << "\nEdges: ";
        for (auto edge : edges) {
            std::cout << "(" << std::get<0>(edge) << ", " << std::get<1>(edge) << ", " << std::get<2>(edge) << ", " << std::get<3>(edge) << ") ";
        }

        std::cout << "\nDFS order: ";
        for (const auto& pair : dfs_order) {
            std::cout << "(" << pair.first << ": " << pair.second << ") ";
        }

        std::cout << "\nAdjacency List: ";
        for (const auto& pair : adjacency) {
            std::cout << pair.first << ": ";
            for (int v : pair.second) {
                std::cout << v << " ";
            }
            std::cout << "| ";
        }

        std::cout << "\n";
    }
};

class Graph {
public:
    std::unordered_map<int, int> vertices;
    std::unordered_map<int, std::vector<std::tuple<int, int, double>>> adjacency;
    double utility;

    Graph() : utility(-999) {}

    double graph_utility(const std::unordered_map<std::tuple<int, int, int>, double>& external_utility) 
    {
        if (utility > -999) 
        {
            return utility;
        }

        double total_utility = 0.0;

        for (const auto& u_entry : adjacency) {
            int u = u_entry.first;

            for (const auto& e : u_entry.second) {
                int v, edge_label;
                double internal_utility;
                std::tie(v, edge_label, internal_utility) = e;

                if (u > v) {
                    continue;
                }

                int u_label = vertices[u];
                int v_label = vertices[v];

                if (u_label <= v_label) {
                    total_utility += internal_utility * external_utility.at(std::make_tuple(u_label, v_label, edge_label));
                } else {
                    total_utility += internal_utility * external_utility.at(std::make_tuple(v_label, u_label, edge_label));
                }
            }
        }

        utility = total_utility;
        return total_utility;
    }

    void print() const {
        std::cout << "Vertices: ";
        for (const auto& entry : vertices) {
            std::cout << "(" << entry.first << ": " << entry.second << ") ";
        }

        std::cout << "\nAdjacency List: ";
        for (const auto& u_entry : adjacency) {
            int u = u_entry.first;
            std::cout << u << ": ";
            for (const auto& e : u_entry.second) {
                int v, edge_label;
                double internal_utility;
                std::tie(v, edge_label, internal_utility) = e;
                std::cout << "(" << v << ", " << edge_label << ", " << internal_utility << ") ";
            }
            std::cout << "| ";
        }

        std::cout << "\n";
    }
};

int main()
{
    std::vector<InputGraph> input_graphs;
    std::unordered_map<std::tuple<int, int, int>, double> external_utility;

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


    // system("Pause");
    return 0;
}









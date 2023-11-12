#ifndef __GWU_H__
#define __GWU_H__

#include "main.h"


// Hash function for std::tuple
namespace std {

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
    struct hash<std::tuple<int, int, int, int, int>> {
        size_t operator()(const std::tuple<int, int, int, int, int>& key) const {
            size_t hash_val = 0;
            hash_val ^= std::hash<int>{}(std::get<0>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<int>{}(std::get<1>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<int>{}(std::get<2>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<int>{}(std::get<3>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= std::hash<int>{}(std::get<4>(key)) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            return hash_val;
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

// DECLARE GLOBAL VARIABLES
static std::unordered_map<std::tuple<int, int, int>, double> external_utility;
static int hup = 0;
static int candidates = 0;

// DECLARE FUNCTIONS
std::pair<int, std::vector<int>> RightMostPath(const std::vector<std::tuple<int, int, int, int, int>>& code);
double get_utility(const std::vector<std::tuple<int, int, int, int, int>>& code,
                   const Graph& graph,
                   const std::unordered_map<int, int>& isomorphism);

std::unordered_map<std::tuple<int, int, int, int, int>, std::pair<double, double>> RightMostExtensions(const std::vector<std::tuple<int, int, int, int, int>>& code,
                    std::vector<Graph>& graphs);

Graph buildGraph(const std::vector<std::tuple<int, int, int, int, int>>& code);

std::tuple<int, int, int, int, int> minTuple(const std::tuple<int, int, int, int, int>& tuple1,
                                            const std::tuple<int, int, int, int, int>& tuple2);

std::tuple<int, int, int, int, int> minExtension(const std::vector<std::tuple<int, int, int, int, int>>& tuples);

bool isCannonical(const std::vector<std::tuple<int, int, int, int, int>>& code);

std::vector<std::unordered_map<int, int>> subgraphIsomorphisms(const std::vector<std::tuple<int, int, int, int, int>>& code, const Graph& graph);

void GSpan(
    const std::vector<std::tuple<int, int, int, int, int>>& code,
    std::vector<Graph>& graphs,
    double min_util,
    time_t t
);

#endif // __GWU_H__
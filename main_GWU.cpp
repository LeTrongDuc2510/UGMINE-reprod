#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>

// Example usage:
// g++ -std=c++11 mainG.cpp -o mainG
// ./mainG

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

std::unordered_map<std::tuple<int, int, int>, double> external_utility;

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


std::pair<int, std::vector<int>> RightMostPath(const std::vector<std::tuple<int, int, int, int, int>>& code) {
    std::unordered_map<int, int> adj;
    int ur = 0;

    for (const auto& c : code) {
        ur = std::max(ur, std::get<0>(c));
        ur = std::max(ur, std::get<1>(c));

        if (std::get<1>(c) > std::get<0>(c)) {
            adj[std::get<1>(c)] = std::get<0>(c);
        }
    }

    std::vector<int> result{ur};
    int u = ur;

    while (u != 0) {
        u = adj[u];
        result.push_back(u);
    }

    std::reverse(result.begin(), result.end());

    return {ur, result};
}


double get_utility(const std::vector<std::tuple<int, int, int, int, int>>& code,
                   const Graph& graph,
                   const std::unordered_map<int, int>& isomorphism) {
    double result = 0.0;

    for (const auto& c : code) {
        double ext_utility = 0.0;
        int u, v, u_label, v_label, edge_label;
        std::tie(u, v, u_label, v_label, edge_label) = c;

        if (u_label < v_label) {
            ext_utility = external_utility[std::make_tuple(u_label, v_label, edge_label)];
        } else {
            ext_utility = external_utility[std::make_tuple(v_label, u_label, edge_label)];
        }

        int iso_u = isomorphism.at(u);
        int iso_v = isomorphism.at(v);

        double int_utility = 0.0;

        for (const auto& e : graph.adjacency.at(iso_u)) {
            if (std::get<0>(e) == iso_v && std::get<1>(e) == edge_label) {
                int_utility = std::get<2>(e);
            }
        }

        result += ext_utility * int_utility;
    }

    return result;
}


std::unordered_map<std::tuple<int, int, int, int, int>, std::pair<double, double>>
RightMostExtensions(const std::vector<std::tuple<int, int, int, int, int>>& code,
                    std::vector<Graph>& graphs) {

    std::unordered_map<std::tuple<int, int, int, int, int>, std::pair<double, double>> result;


    for (size_t i = 0; i < graphs.size(); ++i) 
    {
        Graph& graph = graphs[i];
        std::unordered_map<std::tuple<int, int, int, int, int>, double> temp_result;

        if(code.empty())
        {
            for(auto u = graph.adjacency.begin(); u != graph.adjacency.end(); u++)
            {
                for (const auto& e : u->second)
                {
                    int v, edge_label;
                    double internal_utility;
                    std::tie(v, edge_label, internal_utility) = e;

                    int a = u->first;
                    int u_label = graph.vertices.at(a);
                    int v_label = graph.vertices.at(v);

                    double utility = 0.0;

                    if (u_label < v_label) {
                        utility = internal_utility * external_utility[{u_label, v_label, edge_label}];
                    } else {
                        utility = internal_utility * external_utility[{v_label, u_label, edge_label}];
                    }

                    auto key = std::make_tuple(0, 1, u_label, v_label, edge_label);

                    if (temp_result.find(key) != temp_result.end()) {
                        temp_result[key] = std::max(utility, temp_result[key]);
                    } else {
                        temp_result[key] = utility;
                    }
                }
            }

        } else {
            auto isomorphisms = subgraphIsomorphisms(code, graph);
            std::pair<int, std::vector<int>> temp;
            temp = RightMostPath(code);
            int u = temp.first;
            std::vector<int> R = temp.second;
            for (const auto& isomorphism : isomorphisms) 
            {
                for (int v : R) 
                {
                    if (u == v) { continue; }

                    int iso_u = isomorphism.at(u);
                    int iso_v = isomorphism.at(v);

                    // std::vector<std::tuple<int, int, double>> graph.adjacency.at(iso_u);
                    for (auto& e : graph.adjacency.at(iso_u))
                    {
                        if (std::get<0>(e) != iso_v) { continue; }

                        int edge_label = std::get<1>(e);
                        bool exists = false;

                        for (const auto& c : code) 
                        {
                            if ((std::get<0>(c) == u && std::get<1>(c) == v && std::get<3>(c) == edge_label) ||
                                (std::get<0>(c) == v && std::get<1>(c) == u && std::get<3>(c) == edge_label)) {
                                exists = true;
                                break;
                            }
                        }

                        if (!exists) {
                            std::vector<std::tuple<int, int, int, int, int>> new_code = code;
                            new_code.emplace_back(u, v, graph.vertices.at(iso_u), graph.vertices.at(iso_v), edge_label);
                            double utility = get_utility(new_code, graph, isomorphism);

                            auto key = std::make_tuple(u, v, graph.vertices.at(iso_u), graph.vertices.at(iso_v), edge_label);

                            if (temp_result.find(key) != temp_result.end()) {
                                temp_result[key] = std::max(utility, temp_result[key]);
                            } else {
                                temp_result[key] = utility;
                            }
                        }

                    }
                    int ur = u;
                    for (int u : R)
                    {
                        iso_u = isomorphism.at(u);
                        for (auto& e : graph.adjacency.at(iso_u))
                        {
                            int iso_v, edge_label, int_utility;
                            std::tie(iso_v, edge_label, int_utility) = e;
                            bool flag = false;
                            for (auto& su_iso : isomorphism)
                            {
                                if (su_iso.second == iso_v)
                                {
                                    flag = true;
                                    break;
                                }
                            }
                            if (flag) { continue; }

                            int u_label = graph.vertices.at(iso_u);
                            int v_label = graph.vertices.at(iso_v);

                            std::vector<std::tuple<int, int, int, int, int>> new_code = code;
                            new_code.push_back(std::make_tuple(u, ur + 1, u_label, v_label, edge_label));

                            std::unordered_map<int, int> isomorphism_ = isomorphism;
                            isomorphism_[ur + 1] = iso_v;

                            double utility = get_utility(new_code, graph, isomorphism_);
                            std::tuple<int, int, int, int, int> key = std::make_tuple(u, ur + 1, u_label, v_label, edge_label);

                            if (temp_result.find(key) != temp_result.end()) {
                                temp_result[key] = std::max(utility, temp_result[key]);
                            } else 
                            {
                                temp_result[key] = utility;
                            }

                        }
                    }
                }
            }
        }


        for(auto a = temp_result.begin(); a!= temp_result.end(); a++)
        {
            int x = std::get<0>(a->first);
            int y = std::get<1>(a->first);
            int u_label = std::get<2>(a->first);
            int v_label = std::get<3>(a->first);
            int edge_label = std::get<4>(a->first);
            double utility = a->second;
            auto key = std::make_tuple(x , y, u_label, v_label, edge_label);
            if (result.find(key) != result.end()) 
            {
                std::pair<double, double> temp;
                temp = result[key];
                double existing_utility = temp.first;
                double gwu = temp.second;
                // auto [existing_utility, gwu] = result[key];
                result[key] = {utility + existing_utility, gwu + graph.graph_utility(external_utility)};
            } else {
                result[key] = {utility, graph.graph_utility(external_utility)};
            }
        }        

    }

    return result;
}


Graph buildGraph(const std::vector<std::tuple<int, int, int, int, int>>& code) {
    Graph graph;

    for (const auto& tuple : code) {
        int u, v, u_label, v_label, edge_label;
        std::tie(u, v, u_label, v_label, edge_label) = tuple;

        graph.vertices[u] = u_label;
        graph.vertices[v] = v_label;

        if (graph.adjacency.find(u) != graph.adjacency.end()) {
            graph.adjacency[u].push_back(std::make_tuple(v, edge_label, 0.0));
        } else {
            graph.adjacency[u] = {std::make_tuple(v, edge_label, 0.0)};
        }

        if (graph.adjacency.find(v) != graph.adjacency.end()) {
            graph.adjacency[v].push_back(std::make_tuple(u, edge_label, 0.0));
        } else {
            graph.adjacency[v] = {std::make_tuple(u, edge_label, 0.0)};
        }
    }

    return graph;
}


std::tuple<int, int, int, int, int> minTuple(const std::tuple<int, int, int, int, int>& tuple1,
                                            const std::tuple<int, int, int, int, int>& tuple2) 
{
    int u1, v1, u1_label, v1_label, edge1label;
    std::tie(u1, v1, u1_label, v1_label, edge1label) = tuple1;

    int u2, v2, u2_label, v2_label, edge2label;
    std::tie(u2, v2, u2_label, v2_label, edge2label) = tuple2;

    if (u1 == u2 && v1 == v2) {
        if (u1_label < u2_label) {
            return tuple1;
        } else if (u1_label > u2_label) {
            return tuple2;
        } else if (v1_label < v2_label) {
            return tuple1;
        } else if (v1_label > v2_label) {
            return tuple2;
        } else if (edge1label < edge2label) {
            return tuple1;
        }
        return tuple2;
    } else {
        if (u1 < v1 && u2 < v2) {  // both forward edge
            if (v1 < v2) {
                return tuple1;
            } else if (v1 == v2 && u1 > u2) {
                return tuple1;
            }
            return tuple2;
        }
        if (u1 > v1 && u2 > v2) {  // both backward edge
            if (u1 < u2) {
                return tuple1;
            } else if (u1 == u2 && v1 < v2) {
                return tuple1;
            }
            return tuple2;
        }
        if (u1 < v1 && u2 > v2) {  // tuple1 forward tuple2 backward
            if (v1 <= u2) {
                return tuple1;
            }
            return tuple2;
        }
        if (u1 > v1 && u2 < v2) {  // tuple1 backward tuple2 forward
            if (u1 < v2) {
                return tuple1;
            }
            return tuple2;
        }
    }
    // Return a default value or handle the case when no condition is met
    return std::make_tuple(0, 0, 0, 0, 0);
}

std::tuple<int, int, int, int, int> minExtension(const std::vector<std::tuple<int, int, int, int, int>>& tuples) 
{
    if (tuples.empty()) {
        // Handle the case when the input vector is empty
        return std::make_tuple(0, 0, 0, 0, 0);  // Adjust the default values accordingly
    }

    std::tuple<int, int, int, int, int> result = tuples.front();

    for (const auto& t : tuples) 
    {
        result = minTuple(result, t);
    }

    return result;
}

bool isCannonical(const std::vector<std::tuple<int, int, int, int, int>>& code) {
    Graph graph = buildGraph(code);
    std::vector<std::tuple<int, int, int, int, int>> c;
    std::vector<Graph> temp;
    temp.push_back(graph);
    for (std::size_t i = 0; i < code.size(); ++i) 
    {
        std::vector<std::tuple<int, int, int, int, int>> extension_vec;
        for(auto& tam : RightMostExtensions(c, temp))
        {
            extension_vec.push_back(tam.first);
        }
        std::tuple<int, int, int, int, int> extension;
        extension = minExtension(extension_vec);

        
        if (minTuple(extension, code[i]) != code[i]) {
            return false;  // prune not min DFS code
        }
        c.push_back(extension);
    }

    return true;
}


std::vector<std::unordered_map<int, int>> subgraphIsomorphisms(const std::vector<std::tuple<int, int, int, int, int>>& code, const Graph& graph) {
    std::vector<std::unordered_map<int, int>> isomorphisms;
    int l0 = std::get<2>(code[0]);
    
    for (const auto& vertex : graph.vertices) {
        // Find vertices with the same label as the first vertex in code
        if (vertex.second == l0) {
            isomorphisms.push_back({{0, vertex.first}});
        }
    }

    for (const auto& tuple : code) {
        int u, v, u_label, v_label, edge_label;
        std::tie(u, v, u_label, v_label, edge_label) = tuple;

        std::vector<std::unordered_map<int, int>> temp_isomorphisms;

        for (const auto& isomorphism : isomorphisms) {
            if (v > u) {
                int iso_u = isomorphism.at(u);
                try {
                    auto _ = graph.adjacency.at(iso_u);
                } catch (const std::out_of_range&) {
                    continue;
                }

                for (const auto& e : graph.adjacency.at(iso_u)) {
                    int iso_v, iso_edge_label, _;
                    std::tie(iso_v, iso_edge_label, _) = e;

                    // Check if iso_v is in the isomorphism map
                    auto iso_v_iter = std::find_if(isomorphism.begin(), isomorphism.end(), [iso_v](const auto& kv) { return kv.second == iso_v; });

                    if (iso_v_iter == isomorphism.end() && graph.vertices.at(iso_v) == v_label && edge_label == iso_edge_label) {
                        std::unordered_map<int, int> new_iso = isomorphism;
                        new_iso[v] = iso_v;
                        temp_isomorphisms.push_back(new_iso);
                    }
                }
            } else {
                int iso_u = isomorphism.at(u);
                int iso_v = isomorphism.at(v);

                for (const auto& e : graph.adjacency.at(iso_u)) {
                    int c_iso_v, c_iso_edge_label, _;
                    std::tie(c_iso_v, c_iso_edge_label, _) = e;

                    if (c_iso_v == iso_v && edge_label == c_iso_edge_label) {
                        temp_isomorphisms.push_back(isomorphism);
                    }
                }
            }
        }
        isomorphisms = temp_isomorphisms;
    }

    return isomorphisms;
}

int hup = 0;
int candidates = 0;

void GSpan(
    const std::vector<std::tuple<int, int, int, int, int>>& code,
    std::vector<Graph>& graphs,
    double min_util,
    time_t t
) {

    std::unordered_map<std::tuple<int, int, int, int, int>, std::pair<double, double>> extensions = RightMostExtensions(code, graphs);
    for(auto& item : extensions)
    {
        std::tuple<int, int, int, int, int> key = item.first;
        std::pair<double, double> value = item.second;
        double utility = value.first;
        double gwu = value.second;

        std::vector<std::tuple<int, int, int, int, int>> new_code = code;
        new_code.push_back(key);

        if (isCannonical(new_code) && (gwu>min_util))
        {
            if (utility > min_util) {
                for(auto& it : new_code)
                {
                    std::cout<<std::get<0>(it)<<" "<<std::get<1>(it)<<" "<<std::get<2>(it)<<" "<<std::get<3>(it)<<" "<<std::get<4>(it)<<std::endl;
                }
                std::cout<<utility<<" "<<gwu<<std::endl;
                std::cout<<isCannonical(new_code)<<std::endl;
                std::cout<<"-----------------"<<std::endl;
                hup++;
            }
            GSpan(new_code, graphs, min_util, t);
        }
    }

}

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









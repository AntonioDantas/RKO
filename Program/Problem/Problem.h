// *******************************************************************
//      file with specific functions to solve the TSP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"
#include <ranges>

// Global Variables
extern int n; // size of the vector solution
int nodeCount;
int veichesCount;

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------

// struct with node informations
struct TNode
{
    int id;
    double c; // Cost of security
    double p; // Probability
    double e; // Energy of UAV
};

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ----------

static std::vector<std::vector<double>> dist;  // matrix with Euclidean distance
static std::vector<std::vector<double>> distNorm;  // matrix with distance normalized
static std::vector<TNode> node;                // vector of nodes

//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------
// find max
double find_max_matrix(const std::vector<std::vector<double>>& m, size_t nRowsCols)
{
    double maxV = -std::numeric_limits<double>::infinity();
    size_t n = nRowsCols;
    for (size_t i = 0; i < n && i < m.size(); ++i) {
        for (size_t j = 0; j < n && j < m[i].size(); ++j) {
            double val = m[i][j];
            if (!std::isfinite(val)) continue;
            if (std::fabs(val) <= 0.0) continue;
            if (val > maxV) maxV = val;
        }
    }
    return maxV;
}

// normaliza in-place
void normalize_matrix(std::vector<std::vector<double>>& m, size_t nRowsCols)
{
    auto maxV = find_max_matrix(m, nRowsCols);

    if (!std::isfinite(maxV)) {
        return;
    }

    if (maxV <= 0.0) {
        return;
    }

    size_t n = nRowsCols;
    for (size_t i = 0; i < n && i < m.size(); ++i) {
        for (size_t j = 0; j < n && j < m[i].size(); ++j) {
            double val = m[i][j];
            if (!std::isfinite(val)) continue;
            if (std::fabs(val) <= 0.0) continue;
            m[i][j] = val / maxV;
        }
    }
}

// normalize NodeType from nodes
template<typename NodeType>
void normalize_member(std::vector<NodeType>& nodes, double NodeType::* member, size_t nCount)
{
    double maxV = -std::numeric_limits<double>::infinity();

    size_t n = std::min(nCount, nodes.size());
    for (size_t i = 0; i < n; ++i) {
        double val = nodes[i].*member;
        if (!std::isfinite(val)) continue;
        if (val > maxV) maxV = val;
    }

    if (!std::isfinite(maxV)) return;
    if (maxV <= 0.0) return;

    for (size_t i = 0; i < n; ++i) {
        double val = nodes[i].*member;
        if (!std::isfinite(val)) continue;
        nodes[i].*member = (val) / maxV;
    }
}

// print matrix for debug
void print_matrix(std::vector<std::vector<double>> m, size_t n) {
    for (size_t i = 0; i < n && i < m.size(); ++i) {
        for (size_t j = 0; j < n && j < m[i].size(); ++j) {
            std::cout << std::fixed << std::setprecision(3) << m[i][j] << " ";
        }
        std::cout << "\n";
    }
}

// print nodes detail for debug
void print_nodes(std::vector<TNode> nodes, size_t n) {
    for (size_t i = 0; i < n && i < nodes.size(); ++i) {
        std::cout << "Node " << i
                  << " | id=" << nodes[i].id
                  << " p=" << nodes[i].p
                  << " c=" << nodes[i].c
                  << " energy=" << nodes[i].e
                  << "\n";
    }
}

/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[])
{
    nodeCount = 0;
    veichesCount = 0;
    char name[200] = "../Instances/";
    strcat(name, nameTable);

    FILE *arq;
    arq = fopen(name, "r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    // read instance head
    char temp[100];
    fgets(temp, sizeof(temp), arq);
    //printf("\n%s", temp);

    fscanf(arq, "%d", &nodeCount);
    //printf("Number of nodes: %d\n", nodeCount);

    fscanf(arq, "%s", temp);
    //printf("\n%s", temp);

    fscanf(arq, "%d", &veichesCount);
    //printf("Number of nodes: %d\n", nodeCount);

    // read node informations
    fscanf(arq, "%s", temp);
    //printf("\n%s\n", temp);
    node.clear();
    TNode nodeTemp;

    while (node.size() < nodeCount && !feof(arq))
    {
        fscanf(arq, "%d %lf %lf %lf", &nodeTemp.id, &nodeTemp.c, &nodeTemp.p, &nodeTemp.e);
        //printf("Node read: id=%d c=%.3lf p=%.3lf energy=%d\n", nodeTemp.id, nodeTemp.c, nodeTemp.p, nodeTemp.e);
        node.push_back(nodeTemp);
    }

    // read matrix of distances
    dist.clear();
    distNorm.clear();
    dist.resize(nodeCount, std::vector<double>(nodeCount));
    distNorm.resize(nodeCount, std::vector<double>(nodeCount));
    
    fscanf(arq, "%s", temp);
    //printf("\n%s\n", temp);

    for (int i = 0; i < nodeCount; i++) {
        for (int j = 0; j < nodeCount; j++) {
            if (fscanf(arq, "%lf", &dist[i][j]) != 1) {
                printf("Erro de leitura no elemento [%d][%d]\n", i, j);
            }
            distNorm[i][j] = dist[i][j];
        }
    }

    fclose(arq);
    
    normalize_matrix(distNorm, nodeCount);
    normalize_member(node, &TNode::p, nodeCount);
    normalize_member(node, &TNode::c, nodeCount);

    //std::cout << "Matriz:\n";
    //print_matrix(dist, nodeCount);

    //std::cout << "\nNodes:\n";
    //print_nodes(node, nodeCount);

    n = nodeCount + (nodeCount - veichesCount); // size of the vector solution
}

double calculate_route_distance(const std::vector<int>& route, const std::vector<std::vector<double>>& dist) {
    double total_distance = 0.0;
    
    // Uma rota com 0 ou 1 ponto tem distância 0.
    if (route.size() < 2) {
        return 0.0;
    }

    // Itera pelos pares de nós consecutivos na rota
    for (size_t i = 0; i < route.size() - 1; ++i) {
        int from_node = route[i];
        int to_node = route[i + 1];
        total_distance += dist[from_node][to_node];
    }
    
    return total_distance;
}

/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solutions into problem solutions
*************************************************************************************/
double Decoder(TSol& s)
{
    // penality for infeasible solutions
    const double PENALTY_FACTOR = 99999.0; 

    // Clean the last one solution
    s.routes.clear();
    s.ofv = 0.0;
    s.f1 = 0.0; 
    s.f2 = 0.0; 

    // Variables
    int last_vehicle_id = -1;
    std::vector<int> temp_target_pool;
    std::vector<int> processing_order(nodeCount);
    std::vector<int> extra_visits_pool;

    // base on rk order
    std::iota(processing_order.begin(), processing_order.end(), 0);
    std::sort(processing_order.begin(), processing_order.end(),
              [&](int a, int b) { return s.rk[a] < s.rk[b]; });

    // Empty routes
    for (int i = 0; i < nodeCount; ++i) {
        if (node[i].e > 0) { 
            s.routes.push_back({i});
        }
    }

    //Group nodes in routes
    for (int node_id : processing_order) {
        if (node[node_id].e == 0) { 
            temp_target_pool.push_back(node_id);
        } else {
            auto it = std::find_if(s.routes.begin(), s.routes.end(), 
                                   [&](const std::vector<int>& route){ return route[0] == node_id; });
            if (it != s.routes.end()) {
                it->insert(it->end(), temp_target_pool.begin(), temp_target_pool.end());
            }
            temp_target_pool.clear();
            last_vehicle_id = node_id;
        }
    }

    // insert pool of targets at the last vehicle
    if (!temp_target_pool.empty() && last_vehicle_id != -1) {
        auto it = std::find_if(s.routes.begin(), s.routes.end(), 
                               [&](const std::vector<int>& route){ return route[0] == last_vehicle_id; });
        if (it != s.routes.end()) {
            it->insert(it->end(), temp_target_pool.begin(), temp_target_pool.end());
        }
    }
    
    // Add extra visits to the pool based on rk values
    int count_extra_visits = 0;
    for (int i = 0; i < nodeCount; ++i) {
        if (node[i].e == 0) {
            double revisit_key = s.rk[nodeCount + count_extra_visits];
            int extra_visits = floor(revisit_key * VMAX);
            for (int k = 0; k < extra_visits; ++k) {
                extra_visits_pool.push_back(i);
            }
            count_extra_visits++;
        }
    }

    //Order extra visits by rk value
    std::sort(extra_visits_pool.begin(), extra_visits_pool.end(),
          [&](int a, int b) { return s.rk[a] < s.rk[b]; });

    // Distribute extra visits to vehicles
    if (!extra_visits_pool.empty()) {
        for (size_t i = 0; i < extra_visits_pool.size(); ++i) {
            int target_to_add = extra_visits_pool[i];
            
            auto it = std::min_element(s.routes.begin(), s.routes.end(),
                [&dist](const std::vector<int>& routeA, const std::vector<int>& routeB) {
                    return calculate_route_distance(routeA, dist) < calculate_route_distance(routeB, dist);
                });

            if (it != s.routes.end()) {
                if (it->back() != target_to_add) {
                    it->push_back(target_to_add);
                }
            }
        }
    }

    //validate and calculate objective functions
    double total_energy_overflow = 0.0;
    for (const auto& route : s.routes) {
        if (route.size() <= 1) continue;

        double current_distance = 0.0;
        double current_distance_norm = 0.0;
        for (size_t i = 1; i < route.size(); ++i) {
            int from_node = route[i-1];
            int to_node = route[i];
            current_distance += dist[from_node][to_node];
            current_distance_norm += distNorm[from_node][to_node];
            s.f2 += node[to_node].p;
        }
        s.f1 += current_distance_norm;

        if (current_distance > node[route[0]].e) {
            total_energy_overflow += (current_distance - node[route[0]].e);
        }
    }

    // Objective function with penalty if necessary
    s.ofv = (ALPHA * s.f1) - ((1.0 - ALPHA) * s.f2);
    if (total_energy_overflow > 0) {
        s.ofv += PENALTY_FACTOR * total_energy_overflow;
    }

    if (debug && print) {
        printf("\nConstructed Routes:\n");
        for (const auto& route : s.routes) {
            printf("Vehicle %d -> ", node[route[0]].id);
            for (size_t i = 1; i < route.size(); ++i) {
                printf("%d ", node[route[i]].id);
            }
            printf("\n");
        }
    }
    
    if (!debug && print)
    {
        for (const auto& route : s.routes) {
            fprintf(arqSol, "%d ", node[route[0]].id);
            for (size_t i = 1; i < route.size(); ++i) {
                fprintf(arqSol, "%d ", node[route[i]].id);
            }
        }
    }

    return s.ofv;
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem()
{
    // specific problem
    dist.clear();
    node.clear();
}

#endif
// *******************************************************************
//      file with specific functions to solve the TSP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"
#include <ranges>

// Global Variables
extern int n; // size of the vector solution

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------

// struct with node informations
struct TNode
{
    int id;
    double x;
    double y;
    double p;
};
struct TVehicle
{
    int id;
    double x;
    double y;
    double e;
};

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ----------

static std::vector<std::vector<double>> dist;  // matrix with Euclidean distance
static std::vector<std::vector<double>> distV; // matrix with Euclidean distance for vehicles
static std::vector<TNode> node;                // vector of nodes
static std::vector<TVehicle> vehicle;          // vector of vehicles
static double maxDistance;                     // max of distances
static double minDistance;                     // min of distances

//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------

/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
*************************************************************************************/
void ReadData(char nameTable[])
{
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

    // => read data

    // read instance head
    char temp[100];
    fgets(temp, sizeof(temp), arq);

    // read node informations
    int nAux = 0;
    node.clear();
    TNode nodeTemp;

    while (!feof(arq))
    {
        fscanf(arq, "%d %lf %lf %lf", &nodeTemp.id, &nodeTemp.x, &nodeTemp.y, &nodeTemp.p);
        // printf("\n %d %lf %lf %lf", nodeTemp.id, nodeTemp.x, nodeTemp.y, nodeTemp.p);
        node.push_back(nodeTemp);

        nAux++;
    }
    fclose(arq);

    // calculate the euclidean distance
    dist.clear();
    dist.resize(nAux, std::vector<double>(nAux));

    // for points
    maxDistance = 0.0;
    minDistance = INFINITY;

    for (int i = 0; i < nAux; i++)
    {
        for (int j = i; j < nAux; j++)
        {
            double dx = node[j].x - node[i].x;
            double dy = node[j].y - node[i].y;

            dist[i][j] = dist[j][i] = sqrt(dx * dx + dy * dy);

            if (maxDistance < dist[i][j])
            {
                maxDistance = dist[i][j];
            }

            if (minDistance > dist[i][j])
            {
                minDistance = dist[i][j];
            }

            // printf("%.6f ", dist[i][j]);
        }
        // printf("\n");
    }

    n = nAux;
}

/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solutions into problem solutions
*************************************************************************************/
double Decoder(TSol s)
{
    // create an initial list of candidates
    std::vector<int> sC(n);
    for (int i = 0; i < n; i++)
    {
        sC[i] = i;
    }

    // printf("\n Decoder Antes");
    // for (int i=0; i<n; i++)
    //     printf("\n %d %f %d", sC[i], s.rk[i], node[sC[i]].id);

    // sort the problem vector based on the values in the rk vector
    std::sort(sC.begin(), sC.end(), [&s](int i1, int i2)
              { return s.rk[i1] < s.rk[i2]; });

    // printf("\n Decoder Depois");
    // for (int i=0; i<n; i++)
    //     printf("\n %d %f %d", sC[i], s.rk[i], node[sC[i]].id);

    // problem solution
    std::vector<int> sol;
    s.ofv = 0.0;

    int lastNode = -1;
    int current = 0;
    int currentVehicle = -1;
    double currentDistance = 0;

    while (sol.size() < n)
    {
        if (current > n - 1)
        {
            current = 0;
        }

        // printf("\n current %d lastNode %d currentVehicle %d nodeid %d", current, lastNode, currentVehicle, node[sC[current]].id);

        if (node[sC[current]].id > 1000) // isVehicle
        {
            currentDistance = 0;
            currentVehicle = lastNode = current;
            sol.push_back(sC[current]);
            current++;
            continue;
        }

        if (currentVehicle == -1)
        {
            current++;
            continue;
        }

        float distance = dist[sC[lastNode]][sC[current]];
        currentDistance += distance;

        if (currentDistance > node[sC[currentVehicle]].p)
        {
            return maxDistance * n;
        }
        else
        {
            s.ofv += distance;
            sol.push_back(sC[current]);
        }

        lastNode = current;
        current++;
    }

    // print the solution in the screen
    if (debug && print)
    {
        printf("\n");
        for (int i = 0; i < n; i++)
            printf("%d ", node[sol[i]].id);
    }

    // print the solution in a file
    if (!debug && print)
    {
        for (int i = 0; i < n; i++)
            fprintf(arqSol, "%d ", sol[i]);
    }

    // return the objective function value
    // printf("\n objective function value %f \n\n", s.ofv);
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
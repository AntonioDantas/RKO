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

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  ----------

static std::vector<std::vector<double>> dist;  // matrix with Euclidean distance
static std::vector<TNode> node;                // vector of nodes

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
    for (int i = 0; i < nAux; i++)
    {
        for (int j = i; j < nAux; j++)
        {
            double dx = node[j].x - node[i].x;
            double dy = node[j].y - node[i].y;

            dist[i][j] = dist[j][i] = sqrt(dx * dx + dy * dy);
        }
    }
    
    double minDist = HUGE_VAL;
    double maxDist = -HUGE_VAL;
    for (int i = 0; i < nAux; i++) {
        for (int j = 0; j < nAux; j++) {
            minDist = std::min(minDist, dist[i][j]);
            maxDist = std::max(maxDist, dist[i][j]);
        }
    }

    double range = maxDist - minDist;
    //printf("\nMin Dist: %lf, Max Dist: %lf, Range Dist: %lf", minDist, maxDist, range);
    if (range > 0) {
        for (int i = 0; i < nAux; i++) {
            for (int j = 0; j < nAux; j++) {
                dist[i][j] = (dist[i][j] - minDist) / range;
            }
        }
    }

    double minP = HUGE_VAL;
    double maxP = -HUGE_VAL;
    
    for (int i = 0; i < nAux; i++) {
        if (node[i].id <= 1000) { 
            minP = std::min(minP, node[i].p);
            maxP = std::max(maxP, node[i].p);
        }
    }

    double rangeP = maxP - minP;
    //printf("\nMin P: %lf, Max P: %lf, Range P: %lf", minP, maxP, rangeP);
    if (rangeP > 0) {
        for (int i = 0; i < nAux; i++) {
            if (node[i].id <= 1000) {
                node[i].p = (node[i].p - minP) / rangeP;
            }
        }
    }

    n = nAux;
}

/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solutions into problem solutions
*************************************************************************************/
double Decoder(TSol& s)
{
    // create an initial list of candidates
    std::vector<int> sC(n);
    for (int i = 0; i < n; i++)
    {
        sC[i] = i;
    }

    // sort the problem vector based on the values in the rk vector
    std::sort(sC.begin(), sC.end(), [&s](int i1, int i2)
              { return s.rk[i1] < s.rk[i2]; });

    // problem solution
    std::vector<int> sol;
    s.ofv = 0.0;
    s.f1 = 0.0;
    s.f2 = 0.0;

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

        s.f1 += (ALPHA * distance);
        s.f2 += ((1 - ALPHA) * node[sC[current]].p);

        s.ofv += ((ALPHA * distance) - ((1 - ALPHA) * node[sC[current]].p));
        sol.push_back(sC[current]);

        if (currentDistance > node[currentVehicle].p)
        {
            return s.ofv * 9999;
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
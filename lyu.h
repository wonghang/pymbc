#pragma once

#include "bitset.h"
#include "bipartite_graph.h"
 
int MBC(bipartite_graph *G, biclique *C, uint64_t tau_u, uint64_t tau_v);
int MBC_star(bipartite_graph *G, biclique *C, uint64_t tau_u, uint64_t tau_v, int max_iter, unsigned int optimize);

biclique *InitMBC_Star(bipartite_graph *G, int byV);
biclique *InitMBC_Greedy(bipartite_graph *G, int init_iter);
biclique *InitMBC_Prune(bipartite_graph *G);
biclique *InitMBC_Best(bipartite_graph *G, int init_iter);

int lyu_random_init(unsigned int seed);

#pragma once

#include <stdint.h>
#include "bitset.h"

typedef struct {
  uint64_t M;
  uint64_t N;
  bitset *U;
  bitset *V;
} bipartite_graph;

typedef struct {
  uint64_t M;
  uint64_t N;
  bitset U;
  bitset V;
} biclique;

bipartite_graph *bipartite_graph_alloc(uint64_t M,uint64_t N);
bipartite_graph *bipartite_graph_from_C_bool_array(uint64_t M,uint64_t N,bitset U, bitset V, int8_t **edge);
void bipartite_graph_free(bipartite_graph *G);
bipartite_graph *bipartite_graph_dup(bipartite_graph *G);
bipartite_graph *bipartite_graph_subgraph(bipartite_graph *G, bitset U, bitset V);
int bipartite_graph_copy(bipartite_graph *D, bipartite_graph *S);
void bipartite_graph_transpose(bipartite_graph *G);

uint64_t bipartite_graph_U_length(bipartite_graph *G);
uint64_t bipartite_graph_V_length(bipartite_graph *G);
uint64_t bipartite_graph_get_U_array(bipartite_graph *G, uint64_t *out);
uint64_t bipartite_graph_get_V_array(bipartite_graph *G, uint64_t *out);
uint64_t bipartite_graph_get_U_bitset(bipartite_graph *G, bitset out);
uint64_t bipartite_graph_get_V_bitset(bipartite_graph *G, bitset out);
int bipartite_graph_get_U_at(bipartite_graph *G, uint64_t pos, uint64_t *out);
int bipartite_graph_get_V_at(bipartite_graph *G, uint64_t pos, uint64_t *out);
int bipartite_graph_random_U(bipartite_graph *G, uint64_t *out, struct random_data *buf);
int bipartite_graph_random_V(bipartite_graph *G, uint64_t *out, struct random_data *buf);

uint64_t bipartite_graph_deg_u(bipartite_graph *G, uint64_t u);
uint64_t bipartite_graph_deg_v(bipartite_graph *G, uint64_t v);
uint64_t bipartite_graph_max_deg_u(bipartite_graph *G);
uint64_t bipartite_graph_max_deg_v(bipartite_graph *G);

int bipartite_graph_remove_u(bipartite_graph *G, uint64_t u);
int bipartite_graph_remove_v(bipartite_graph *G, uint64_t v);
int bipartite_graph_is_complete(bipartite_graph *G);

#ifndef NDEBUG
void bipartite_graph_dump(bipartite_graph *G);
#endif

biclique *biclique_alloc(uint64_t M,uint64_t N);
void biclique_free(biclique *B);
void biclique_clear(biclique *B);
uint64_t biclique_U_length(biclique *B);
uint64_t biclique_V_length(biclique *B);
uint64_t biclique_size(biclique *B);
int biclique_copy(biclique *D, biclique *S);
void biclique_transpose(biclique *B);
void biclique_add_U(biclique *B, uint64_t u);
void biclique_add_V(biclique *B, uint64_t v);
void biclique_set_U(biclique *B, bitset U);
void biclique_set_V(biclique *B, bitset V);

#ifndef NDEBUG
int biclique_is_biclique(bipartite_graph *G, biclique *B);
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifndef NDEBUG
#include <stdio.h>
#endif

#include "bipartite_graph.h"

bipartite_graph *bipartite_graph_alloc(uint64_t M,uint64_t N)
{
  bipartite_graph *G;
  assert(M > 0);
  assert(N > 0);
  
  G = malloc(sizeof(bipartite_graph));
  G->U = NULL;
  G->V = NULL;
  
  G->M = M;
  G->N = N;
  
  G->U = malloc(sizeof(bitset)*M);
  if(G->U == NULL) goto fail;
  memset(G->U, 0, sizeof(bitset)*M);
  
  G->V = malloc(sizeof(bitset)*N);
  if(G->V == NULL) goto fail;
  memset(G->V, 0, sizeof(bitset)*N);
  
  return G;
  
 fail:
  bipartite_graph_free(G);
  return NULL;
}
bipartite_graph *bipartite_graph_from_C_bool_array(uint64_t M,uint64_t N,bitset U,bitset V, int8_t **edge)
{
  bipartite_graph *G;
  bitset tmp;
  size_t M_ndata, N_ndata;
  uint64_t u,v,k;

  G = bipartite_graph_alloc(M,N);
  if(G == NULL) goto fail;
  
  M_ndata = bitset_required_capacity(M);
  N_ndata = bitset_required_capacity(N);
#define HAS_EDGE(i,j) edge[(i)][(j)]
  
  // scan and construct for U direction
  for(u=0;u<M;u++) {
    if(U != NULL && !bitset_is_member(U,u,M_ndata)) continue;

    tmp = NULL;
    for(k=0;k<N;k++) {
      if(bitset_is_member(V,k,N_ndata) && HAS_EDGE(u,k)) {
	if(tmp == NULL) {
	  tmp = bitset_alloc(N_ndata);
	  if(tmp == NULL) goto fail;
	  bitset_clear(tmp,N_ndata);
	}
	bitset_add(tmp,k,N_ndata);
      }
    }
    if(tmp != NULL) G->U[u] = tmp;
  }
  
  // scan and construct for V direction
  for(v=0;v<N;v++) {
    if(V != NULL && !bitset_is_member(V,v,N_ndata)) continue;
    
    tmp = NULL;
    for(k=0;k<M;k++) {
      if(bitset_is_member(U,k,M_ndata) && HAS_EDGE(k,v)) {
	if(tmp == NULL) {
	  tmp = bitset_alloc(M_ndata);
	  if(tmp == NULL) goto fail;
	  bitset_clear(tmp,M_ndata);
	}
	bitset_add(tmp,k,M_ndata);
      }
    }
    if(tmp != NULL) G->V[v] = tmp;
  }
#undef HAS_EDGE
  
  return G;
  
 fail:
  bipartite_graph_free(G);
  return NULL;
}
void bipartite_graph_free(bipartite_graph *G)
{
  uint64_t k;
  
  if(G) {
    if(G->U) {
      for(k=0;k<G->M;k++) bitset_free(G->U[k]);
      free(G->U);
    }
    if(G->V) {
      for(k=0;k<G->N;k++) bitset_free(G->V[k]);
      free(G->V);
    }
    free(G);
  }
}
bipartite_graph *bipartite_graph_dup(bipartite_graph *G)
{
  bipartite_graph *NG;
  uint64_t u,v;
  size_t M_ndata,N_ndata;

  NG = bipartite_graph_alloc(G->M,G->N);
  if(NG == NULL) goto fail;

  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);
  
  for(u=0;u<G->M;u++) {
    if(G->U[u]) {
      NG->U[u] = bitset_alloc(N_ndata);
      if(NG->U[u] == NULL) goto fail;
      bitset_copy(NG->U[u],G->U[u],N_ndata);
    }
  }

  for(v=0;v<G->N;v++) {
    if(G->V[v]) {
      NG->V[v] = bitset_alloc(M_ndata);
      if(NG->V[v] == NULL) goto fail;
      bitset_copy(NG->V[v],G->V[v],M_ndata);
    }
  }

  return NG;

 fail:
  bipartite_graph_free(NG);
  return NULL;
}
bipartite_graph *bipartite_graph_subgraph(bipartite_graph *G, bitset U, bitset V)
{
  bipartite_graph *NG;
  uint64_t u,v;
  size_t M_ndata,N_ndata;

  NG = bipartite_graph_alloc(G->M,G->N);
  if(NG == NULL) goto fail;

  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);
  
  for(u=0;u<G->M;u++) {
    if(G->U[u] && bitset_is_member(U,u,M_ndata)) {
      NG->U[u] = bitset_alloc(N_ndata);
      if(NG->U[u] == NULL) goto fail;
      bitset_copy(NG->U[u],G->U[u],N_ndata);
    }
  }

  for(v=0;v<G->N;v++) {
    if(G->V[v] && bitset_is_member(V,v,N_ndata)) {
      NG->V[v] = bitset_alloc(M_ndata);
      if(NG->V[v] == NULL) goto fail;
      bitset_copy(NG->V[v],G->V[v],M_ndata);
    }
  }

  return NG;

 fail:
  bipartite_graph_free(NG);
  return NULL;
}

int bipartite_graph_copy(bipartite_graph *D, bipartite_graph *S)
{
  uint64_t u,v;
  uint64_t M,N;
  size_t M_ndata,N_ndata;
  
  if(D->M != S->M || D->N != S->N) return -1;

  M = S->M;
  N = S->N;
  
  M_ndata = bitset_required_capacity(M);
  N_ndata = bitset_required_capacity(N);
  
  for(u=0;u<M;u++) {
    if(S->U[u]) {
      if(D->U[u] == NULL) {
	D->U[u] = bitset_alloc(N_ndata);
	if(D->U[u] == NULL) return -1;
      }
      bitset_copy(D->U[u],S->U[u],N_ndata);
    }
    else {
      if(D->U[u]) {
	bitset_free(D->U[u]);
	D->U[u] = NULL;
      }
    }
  }

  for(v=0;v<N;v++) {
    if(S->V[v]) {
      if(D->V[v] == NULL) {
	D->V[v] = bitset_alloc(M_ndata);
	if(D->V[v] == NULL) return -1;
      }
      bitset_copy(D->V[v],S->V[v],M_ndata);
    }
    else {
      if(D->V[v]) {
	bitset_free(D->V[v]);
	D->V[v] = NULL;
      }
    }
  }
  return 0;
}
void bipartite_graph_transpose(bipartite_graph *G)
{
  uint64_t tmpI;
  bitset *tmpS;

  if(G == NULL) return;

  tmpI = G->M;
  G->M = G->N;
  G->N = tmpI;

  tmpS = G->U;
  G->U = G->V;
  G->V = tmpS;
}
uint64_t bipartite_graph_U_length(bipartite_graph *G)
{
  uint64_t c,i;
  c = 0;
  for(i=0;i<G->M;i++) {
    if(G->U[i] != NULL) c++;
  }
  return c;
}
uint64_t bipartite_graph_V_length(bipartite_graph *G)
{
  uint64_t c,i;
  c = 0;
  for(i=0;i<G->N;i++) {
    if(G->V[i] != NULL) c++;
  }
  return c;
}
uint64_t bipartite_graph_get_U_array(bipartite_graph *G, uint64_t *out)
{
  uint64_t c,i;
  c = 0;
  for(i=0;i<G->M;i++) {
    if(G->U[i] != NULL) *out++ = i;
  }
  return c;
}
uint64_t bipartite_graph_get_V_array(bipartite_graph *G, uint64_t *out)
{
  uint64_t c,i;
  c = 0;
  for(i=0;i<G->N;i++) {
    if(G->V[i] != NULL) *out++ = i;
  }
  return c;
}
uint64_t bipartite_graph_get_U_bitset(bipartite_graph *G, bitset out)
{
  uint64_t c,i;
  size_t ndata;
  c = 0;

  ndata = bitset_required_capacity(G->M);
  bitset_clear(out,ndata);
  for(i=0;i<G->M;i++) {
    if(G->U[i] != NULL) bitset_add(out, i, ndata);
  }
  return c;
}
uint64_t bipartite_graph_get_V_bitset(bipartite_graph *G, bitset out)
{
  uint64_t c,i;
  size_t ndata;
  c = 0;

  ndata = bitset_required_capacity(G->N);
  bitset_clear(out,ndata);
  for(i=0;i<G->N;i++) {
    if(G->V[i] != NULL) bitset_add(out, i, ndata);
  }
  return c;
}
int bipartite_graph_get_U_at(bipartite_graph *G, uint64_t pos, uint64_t *out)
{
  uint64_t c,i;
  c = 0;
  for(i=0;i<G->M;i++) {
    if(G->U[i] != NULL) {
      if(c == pos) {
	*out = i;
	return 1;
      }
      c++;
    }
  }
  return 0;
}
int bipartite_graph_get_V_at(bipartite_graph *G, uint64_t pos, uint64_t *out)
{
  uint64_t c,i;
  c = 0;
  for(i=0;i<G->N;i++) {
    if(G->V[i] != NULL) {
      if(c == pos) {
	*out = i;
	return 1;
      }
      c++;
    }
  }
  return 0;
}
int bipartite_graph_random_U(bipartite_graph *G, uint64_t *out, struct random_data *buf)
{
  uint64_t length;
  int32_t randval;
  uint64_t pos;

  length = bipartite_graph_U_length(G);
  if(length == 0) return 0;
  if(length > RAND_MAX) return -1;
  if(random_r(buf,&randval)<0) return -1;

  pos = ((uint64_t)randval)%length;
  return bipartite_graph_get_U_at(G,pos,out);
}
int bipartite_graph_random_V(bipartite_graph *G, uint64_t *out, struct random_data *buf)
{
  uint64_t length;
  int32_t randval;
  uint64_t pos;

  length = bipartite_graph_V_length(G);
  if(length == 0) return 0;
  if(length > RAND_MAX) return -1;
  if(random_r(buf,&randval)<0) return -1;

  pos = ((uint64_t)randval)%length;
  return bipartite_graph_get_V_at(G,pos,out);
}

uint64_t bipartite_graph_deg_u(bipartite_graph *G, uint64_t u)
{
  assert(u >= 0);
  if(G->U[u] == NULL) return 0;
  else return bitset_length(G->U[u],bitset_required_capacity(G->N));
}
uint64_t bipartite_graph_deg_v(bipartite_graph *G, uint64_t v)
{
  assert(v >= 0);
  if(G->V[v] == NULL) return 0;
  else return bitset_length(G->V[v],bitset_required_capacity(G->M));
}
uint64_t bipartite_graph_max_deg_u(bipartite_graph *G)
{
  uint64_t r = 0, c;
  uint64_t u;
  for(u=0;u<G->M;u++) {
    c = bipartite_graph_deg_u(G, u);
    if(c > r) r = c;
  }
  return r;
}
uint64_t bipartite_graph_max_deg_v(bipartite_graph *G)
{
  uint64_t r = 0, c;
  uint64_t v;
  for(v=0;v<G->N;v++) {
    c = bipartite_graph_deg_v(G, v);
    if(c > r) r = c;
  }
  return r;
}
int bipartite_graph_remove_u(bipartite_graph *G, uint64_t u)
{
  bitset_iterator it;
  size_t M_ndata, N_ndata;
  int c;
  
  if(G->U[u] == NULL) return 0;
  else {
    M_ndata = bitset_required_capacity(G->M);
    N_ndata = bitset_required_capacity(G->N);
    bitset_iterator_begin(&it, G->U[u], N_ndata);
    c = 0;
    while(bitset_iterator_next(&it)) {
      if(G->V[it.current]) {
	bitset_remove(G->V[it.current], u, M_ndata);
	c += 1;
      }
    }
    bitset_free(G->U[u]);
    G->U[u] = NULL;
    return c;
  }
}
int bipartite_graph_remove_v(bipartite_graph *G, uint64_t v)
{
  bitset_iterator it;
  size_t M_ndata, N_ndata;
  int c;
  
  if(G->V[v] == NULL) return 0;
  else {
    M_ndata = bitset_required_capacity(G->M);
    N_ndata = bitset_required_capacity(G->N);
    bitset_iterator_begin(&it, G->V[v], M_ndata);
    c = 0;
    while(bitset_iterator_next(&it)) {
      if(G->U[it.current]) {
	bitset_remove(G->U[it.current], v, N_ndata);
	c += 1;
      }
    }
    bitset_free(G->V[v]);
    G->V[v] = NULL;
    return c;
  }
}
int bipartite_graph_is_complete(bipartite_graph *G)
{
  uint64_t u,v;
  uint64_t U_length,V_length;
  size_t M_ndata,N_ndata;

  if(G == NULL) return 0;

  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);
  U_length = bipartite_graph_U_length(G);
  V_length = bipartite_graph_V_length(G);
  for(u=0;u<G->M;u++) {
    if(G->U[u] == NULL) continue;
    if(bitset_length(G->U[u],N_ndata) != V_length) return 0;
  }
  for(v=0;v<G->N;v++) {
    if(G->V[v] == NULL) continue;
    if(bitset_length(G->V[v],M_ndata) != U_length) return 0;
  }
  return 1;
}
#ifndef NDEBUG
void bipartite_graph_dump(bipartite_graph *G)
{
  uint64_t k;
  char *buf;
  size_t buflen = (G->M > G->N) ? G->M : G->N;

  buf = malloc(buflen+2);
  printf("--- U (%lu) ---\n", G->M);
  for(k=0;k<G->M;k++) {
    if(G->U[k]) {
      bitset_to_string(G->U[k], buf, G->N);
    }
    else {
      memset(buf, '0', G->N);
      buf[G->N] = '*';
      buf[G->N+1] = '\0';
    }
    printf("%lu: %s\n", k, buf);
  }
  printf("--- V (%lu) ---\n", G->N);
  for(k=0;k<G->N;k++) {
    if(G->V[k]) {
      bitset_to_string(G->V[k], buf, G->M);
    }
    else {
      memset(buf, '0', G->M);
      buf[G->M] = '*';
      buf[G->M+1] = '\0';
    }
    printf("%lu: %s\n", k, buf);
  }
  free(buf);
}
#endif

biclique *biclique_alloc(uint64_t M,uint64_t N)
{
  biclique *B;
  B = malloc(sizeof(biclique));
  if(B == NULL) goto fail;
  B->M = M;
  B->N = N;
  B->U = NULL;
  B->V = NULL;

  B->U = bitset_alloc(bitset_required_capacity(M));
  if(B->U == NULL) goto fail;
  B->V = bitset_alloc(bitset_required_capacity(N));
  if(B->V == NULL) goto fail;

  return B;
 fail:
  biclique_free(B);
  return NULL;
}
void biclique_clear(biclique *B)
{
  bitset_clear(B->U, bitset_required_capacity(B->M));
  bitset_clear(B->V, bitset_required_capacity(B->N));
}
void biclique_free(biclique *B)
{
  if(B) {
    if(B->U) bitset_free(B->U);
    if(B->V) bitset_free(B->V);
    free(B);
  }
}
uint64_t biclique_U_length(biclique *B)
{
  return bitset_length(B->U, bitset_required_capacity(B->M));
}
uint64_t biclique_V_length(biclique *B)
{
  return bitset_length(B->V, bitset_required_capacity(B->N));
}
uint64_t biclique_size(biclique *B)
{
  return bitset_length(B->U, bitset_required_capacity(B->M)) * bitset_length(B->V, bitset_required_capacity(B->N));
}
int biclique_copy(biclique *D, biclique *S)
{
  if(D->M != S->M || D->N != S->N) return -1;
  bitset_copy(D->U, S->U, D->M);
  bitset_copy(D->V, S->V, D->N);
  return 0;
}
void biclique_transpose(biclique *B)
{
  uint64_t tmpI;
  bitset tmpS;
  
  if(B == NULL) return;

  tmpI = B->M;
  B->M = B->N;
  B->N = tmpI;

  tmpS = B->U;
  B->U = B->V;
  B->V = tmpS;
}

void biclique_add_U(biclique *B, uint64_t u)
{
  bitset_add(B->U, u, bitset_required_capacity(B->M));
}
void biclique_add_V(biclique *B, uint64_t v)
{
  bitset_add(B->V, v, bitset_required_capacity(B->N));
}
void biclique_set_U(biclique *B, bitset U)
{
  bitset_copy(B->U, U, bitset_required_capacity(B->M));
}
void biclique_set_V(biclique *B, bitset V)
{
  bitset_copy(B->V, V, bitset_required_capacity(B->N));
}
#ifndef NDEBUG
int biclique_is_biclique(bipartite_graph *G, biclique *B)
{
  bitset_iterator it;
  uint64_t M,N;
  size_t M_ndata,N_ndata;

  if(G == NULL || B == NULL) return 0;

  if(G->M != B->M || G->N != B->N) return 0;

  M = G->M;
  N = G->N;
  M_ndata = bitset_required_capacity(M);
  N_ndata = bitset_required_capacity(N);
  bitset_iterator_begin(&it, B->U, M_ndata);
  while(bitset_iterator_next(&it)) {
    if(!bitset_is_subset(B->V, G->U[it.current], N_ndata)) return 0;
  }

  bitset_iterator_begin(&it, B->V, N_ndata);
  while(bitset_iterator_next(&it)) {
    if(!bitset_is_subset(B->U, G->V[it.current], M_ndata)) return 0;
  }
  
  return 1;
}
#endif

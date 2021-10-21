#include <stdlib.h>
#include <assert.h>
#include <string.h>

#ifndef NDEBUG
#include <stdio.h>
#endif

#include "bipartite_graph.h"
#include "lyu.h"

#define FREE_BITSET(x) bitset_free(x); x = NULL

#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define TEST_LYU_STOP if(lyu_stop) goto fail

volatile int lyu_stop = 0;

static struct random_data _internal_random_state;
static char _internal_random_buf[8];

int BranchBound(bipartite_graph *G,
		bitset U, bitset V, bitset C_V, bitset X_V,
		uint64_t tau_u,
		uint64_t tau_v, biclique *C_star)
{
  uint64_t lenU, lenV;
  uint64_t lenU1, lenV1, lenC_V1, C_star_size;
  uint64_t v_star,v;
  size_t M_ndata, N_ndata;
  bitset U1 = NULL, V1 = NULL, C_V1 = NULL, X_V1 = NULL;
  bitset test = NULL;
  bitset N;
  bitset_iterator it;
  int go;

  assert(G->M == C_star->M);
  assert(G->N == C_star->N);

  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);

  lenU = bitset_length(U, M_ndata);
  lenV = bitset_length(V, N_ndata);
  C_star_size = biclique_size(C_star);
  if(lenV >= tau_v && lenU*lenV > C_star_size) {
    biclique_set_U(C_star, U);
    biclique_set_V(C_star, V);
    C_star_size = biclique_size(C_star);
  }
  
  test = bitset_alloc(M_ndata);
  if(test == NULL) goto fail;
  
  U1 = bitset_alloc(M_ndata);
  if(U1 == NULL) goto fail;

  V1 = bitset_alloc(N_ndata);
  if(V1 == NULL) goto fail;
  bitset_clear(V1, N_ndata);
  
  C_V1 = bitset_alloc(N_ndata);
  if(C_V1 == NULL) goto fail;
  bitset_clear(C_V1, N_ndata);
  
  X_V1 = bitset_alloc(N_ndata);
  if(X_V1 == NULL) goto fail;
  bitset_clear(X_V1, N_ndata);

  while(bitset_pop(C_V, &v_star, N_ndata, 0, 0)) {
    TEST_LYU_STOP;
    
    bitset_intersect(U1,U,G->V[v_star],M_ndata);
    if(bitset_length(U1,M_ndata) < tau_u) {
      bitset_free(U1);
      continue;
    }
    bitset_copy(V1, V, N_ndata);
    bitset_add(V1, v_star, N_ndata);
    bitset_clear(C_V1, N_ndata);
    
    bitset_iterator_begin(&it, C_V, N_ndata);
    while(bitset_iterator_next(&it)) {
      if(lyu_stop) goto fail;
      
      v = it.current; // better symbol      
      N = G->V[v];
      if(bitset_is_subset(U1,N,M_ndata)) {
	bitset_add(V1, v, N_ndata);
      }
      else {
	bitset_intersect(test,N,U1,M_ndata);
	if(bitset_length(test,M_ndata) >= tau_u) {
	  bitset_add(C_V1, v, N_ndata);
	}
      }
    }
    
    lenV1 = bitset_length(V1,N_ndata);
    lenC_V1 = bitset_length(C_V1,N_ndata);

    if((lenV1 + lenC_V1) < tau_v) {
      continue;
    }
    lenU1 = bitset_length(U1,M_ndata);
    if(lenU1*(lenV1 + lenC_V1) <= C_star_size) {
      continue;
    }
    bitset_clear(X_V1, N_ndata);
    bitset_iterator_begin(&it, X_V, N_ndata);
    go = 1;
    while(bitset_iterator_next(&it)) {
      if(lyu_stop) goto fail;

      v = it.current;
      N = G->V[v];
      if(bitset_is_subset(U1,N,M_ndata)) {
	go = 0;
	break;
      }
      bitset_intersect(test,N,U1,M_ndata);
      if(bitset_length(test,M_ndata) >= tau_u) {
	bitset_add(X_V1, v, N_ndata);
      }
    }
    if(go) {
      if(BranchBound(G,U1,V1,C_V1,X_V1,tau_u,tau_v,C_star)<0) return -1;
    }
    bitset_add(X_V, v_star, N_ndata);
  }
  
  FREE_BITSET(U1);
  FREE_BITSET(V1);
  FREE_BITSET(C_V1);
  FREE_BITSET(X_V1);
  FREE_BITSET(test);
  return 0;
  
 fail:
  FREE_BITSET(U1);
  FREE_BITSET(V1);
  FREE_BITSET(C_V1);
  FREE_BITSET(X_V1);
  FREE_BITSET(test);
  return -1;
}
int Reduce1Hop(bipartite_graph *G, uint64_t tau_u, uint64_t tau_v)
{
  int c,tc;
  uint64_t u,v;
  uint64_t deg;
  
  c = 1;
  tc = 0;
  while(c > 0) {
    TEST_LYU_STOP;
    
    c = 0;
    // try u
    for(u=0;u<G->M;u++) {
      if(G->U[u] == NULL) continue;
      deg = bipartite_graph_deg_u(G,u);
      if(deg < tau_v) {
	bipartite_graph_remove_u(G,u);
	c++;
	break;
      }
    }

    // try v
    for(v=0;v<G->N;v++) {
      if(G->V[v] == NULL) continue;
      deg = bipartite_graph_deg_v(G,v);
      if(deg < tau_u) {
	bipartite_graph_remove_v(G,v);
	c++;
	break;
      }
    }
    tc += c;
  }
  return tc;

 fail:
  return -1;
}
static int _Reduce2H_optimize1_compare(const void *a, const void *b, void *c)
{
  const uint64_t i = *(const uint64_t *)a;
  const uint64_t j = *(const uint64_t *)b;
  uint64_t *score = (uint64_t *)c;
  if(score[i] < score[j]) return -1;
  else if(score[i] > score[j]) return 1;
  else return 0;
}
void Reduce2H_optimize1(bipartite_graph *G, uint64_t *idx, uint64_t *score)
{
  bitset_iterator it;
  size_t M_ndata,N_ndata;
  uint64_t u,c;
  
  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);
  for(u=0;u<G->M;u++) {
    TEST_LYU_STOP;
    
    idx[u] = u;
    if(G->U[u] == NULL) {
      score[u] = 0;
    }
    else {
      bitset_iterator_begin(&it, G->U[u], N_ndata);
      c = 0;
      while(bitset_iterator_next(&it)) {
	c += bitset_length(G->V[it.current], M_ndata);
      }
      score[u] = c;
    }
  }
  qsort_r(idx, G->M, sizeof(uint64_t), _Reduce2H_optimize1_compare, score);

 fail:
  return;
}
int Reduce2H(bipartite_graph *G, uint64_t tau_u, uint64_t tau_v, unsigned int optimize)
{
  uint64_t u,v,u1;
  bitset_iterator it1,it2;
  size_t M_ndata, N_ndata;
  uint64_t *S = NULL;
  uint64_t i,k,c,tc;

  // optimize1 (early pruning)
  uint64_t *score = NULL, *idx = NULL;
  // optimize2 (early skipping)
  uint64_t *u_c = NULL;
  bitset test = NULL;

  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);

  S = malloc(sizeof(uint64_t)*G->M);
  if(S == NULL) goto fail;
  tc = 0;

  if(optimize & 1U) {
    score = malloc(sizeof(uint64_t)*G->M);
    if(score == NULL) goto fail;
    idx = malloc(sizeof(uint64_t)*G->M);
    if(idx == NULL) goto fail;
    Reduce2H_optimize1(G,idx,score);
  }
  if(optimize & 2U) {
    u_c = malloc(sizeof(uint64_t)*G->M);
    if(u_c == NULL) goto fail;
    memset(u_c, 0 ,sizeof(uint64_t)*G->M);
    test = bitset_alloc(N_ndata);
    if(test == NULL) goto fail;
  }

  for(k=0;k<G->M;k++) {
    TEST_LYU_STOP;

    u = (optimize&1U) ? idx[k] : k;
    if(G->U[u] == NULL) continue;
    if((optimize&2U) && (u_c[u]+1) >= tau_u) {
      /* printf("Skipped %lu by optimize2, u_c[u]=%lu, tau_u=%lu\n", u, u_c[u], tau_u); */
      continue;
    }

    bitset_iterator_begin(&it1, G->U[u], N_ndata);
    memset(S, 0, sizeof(uint64_t)*G->M);
    
    while(bitset_iterator_next(&it1)) {
      TEST_LYU_STOP;

      v = it1.current;
      bitset_iterator_begin(&it2, G->V[v], M_ndata);
      while(bitset_iterator_next(&it2)) {
	TEST_LYU_STOP;

	u1 = it2.current;
	S[u1]++;

	if(optimize&2U) { // optimize 2
	  // page 1366 N_{tau_v} = {u': | N(u',G) \intersects N(u,G) | >= tau_v }
	  bitset_intersect(test, G->U[u], G->U[u1], N_ndata);
	  if(bitset_length(test, N_ndata) >= tau_v) u_c[u1]++;
	}
      }
    }

    c = 0;
    for(i=0;i<G->M;i++) {
      TEST_LYU_STOP;
      if(S[i] >= tau_v) c++;
    }
    if(c < tau_u) {
      bipartite_graph_remove_u(G,u);
      tc++;
    }
  }
  
  if(S) free(S);
  if(score) free(score);
  if(idx) free(idx);
  if(u_c) free(u_c);
  if(test) free(test);
  return tc;

 fail:
  if(S) free(S);
  if(score) free(score);
  if(idx) free(idx);
  if(u_c) free(u_c);
  if(test) free(test);
  return -1;
}
int Reduce2Hop(bipartite_graph *G, uint64_t tau_u, uint64_t tau_v, unsigned int optimize)
{
  int c1,c2;

  c1 = Reduce2H(G,tau_u,tau_v,optimize);
  if(c1 < 0) return -1;

  bipartite_graph_transpose(G);
  c2 = Reduce2H(G,tau_v,tau_u,optimize);
  if(c2 < 0) return -1;
  
  bipartite_graph_transpose(G);
  return c1+c2;
}

int MBC(bipartite_graph *G, biclique *C, uint64_t tau_u, uint64_t tau_v)
{
  bitset U = NULL,V = NULL,e1 = NULL,e2 = NULL;
  size_t M_ndata, N_ndata;
  uint64_t u;
  int ret;
  
  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);
  U = bitset_alloc(M_ndata);
  V = bitset_alloc(N_ndata);
  e1 = bitset_alloc(N_ndata);
  e2 = bitset_alloc(N_ndata);
  if(U == NULL || V == NULL || e1 == NULL || e2 == NULL) goto fail;
  bitset_clear(U, M_ndata);
  bitset_clear(V, N_ndata);
  bitset_clear(e1, N_ndata);
  bitset_clear(e2, N_ndata);

  for(u=0;u<G->M;u++) {
    if(G->U[u]) bitset_add(U, u, M_ndata);
  }
  for(u=0;u<G->N;u++) { // u->v
    if(G->V[u]) bitset_add(V, u, N_ndata);
  }

  ret = BranchBound(G, U, e1, V, e2, tau_u, tau_v, C);
  
  FREE_BITSET(U);
  FREE_BITSET(V);
  FREE_BITSET(e1);
  FREE_BITSET(e2);

  return ret;

 fail:
  FREE_BITSET(U);
  FREE_BITSET(V);
  FREE_BITSET(e1);
  FREE_BITSET(e2);
  return -1;
}
int MBC_star(bipartite_graph *G, biclique *C, uint64_t tau_u, uint64_t tau_v, int max_iter, unsigned int optimize)
{
  uint64_t tau_u1, tau_v1, tau_v0;
  uint64_t C_star_size;
  bipartite_graph *G1;
  int i,c1,c2;
  
  tau_v0 = bipartite_graph_max_deg_u(G);
  G1 = NULL;
  /* printf("tau_v0=%lu, tau_u=%lu, tau_v=%lu\n",tau_v0,tau_u,tau_v); */
  while(tau_v0 > tau_v) {
    C_star_size = biclique_size(C);
    tau_u1 = MAX(C_star_size / tau_v0, tau_u);
    tau_v1 = MAX(tau_v0 / 2UL, tau_v);

    if(G1 == NULL) G1 = bipartite_graph_dup(G);
    else {
      if(bipartite_graph_copy(G1,G)<0) {
	bipartite_graph_free(G1);
	return -1;
      }
    }

    /* printf("Start reduction: |C|=%lu, tau_u1=%lu, tau_v1=%lu\n", C_star_size, tau_u1, tau_v1); */

    for(i=0;i<max_iter;i++) {
      c1 = Reduce1Hop(G1,tau_u1,tau_v1);
      /* printf("%d/%d: Reduce1Hop %d\n", i, max_iter, c1); */
      if(c1 < 0) return -1;
      c2 = Reduce2Hop(G1,tau_u1,tau_v1,optimize);
      /* printf("%d/%d: Reduce2Hop %d\n", i, max_iter, c2); */
      if(c2 < 0) return -1;
      if(c1 == 0 && c2 == 0) break;
    }

    /* printf("Call MBC: |C|=%lu, tau_u1=%lu, tau_v1=%lu\n", C_star_size, tau_u1, tau_v1); */
    if(MBC(G1, C, tau_u1, tau_v1) < 0) {
      if(G1) bipartite_graph_free(G1);
      return -1;
    }

    tau_v0 = tau_v1;
  }
  return 0;
}
int lyu_random_init(unsigned int seed)
{
  return initstate_r(seed, _internal_random_buf, sizeof(_internal_random_buf), &_internal_random_state);
}

biclique *InitMBC_Star(bipartite_graph *G, int byV)
{
  biclique *B = NULL;
  uint64_t u;
  uint64_t best_length, length;
  size_t ndata;
  
  B = biclique_alloc(G->M, G->N);
  if(B == NULL) goto fail;

  biclique_clear(B);

  best_length = 0;
  if(byV) {
    ndata = bitset_required_capacity(G->M);
    
    for(u=0;u<G->N;u++) { // u->v
      TEST_LYU_STOP;
      
      length = bitset_length(G->V[u], ndata);
      if(length == 0) continue; // u->v
      
      if(length > best_length) {
	biclique_clear(B);
	biclique_add_V(B,u); // u->v
	biclique_set_U(B,G->V[u]); // u->v
	best_length = length;
      }
    }
  }
  else {
    ndata = bitset_required_capacity(G->N);

    for(u=0;u<G->M;u++) {
      TEST_LYU_STOP;

      length = bitset_length(G->U[u], ndata);
      if(length == 0) continue;

      if(length > best_length) {
	biclique_clear(B);
	biclique_add_U(B,u);
	biclique_set_V(B,G->U[u]);
	best_length = length;
      }
    }
  }

  if(best_length == 0) goto fail;
  return B;

 fail:
  if(B) biclique_free(B);
  return NULL;
}
static int _FisherYates(uint64_t *a, size_t n)
{
  size_t i,j;
  uint64_t tmp;
  int32_t randval;

  if(n >= RAND_MAX) return -1;

  for(i=n-1;i > 0;i--) {
    if(random_r(&_internal_random_state, &randval)<0) return -1;
    j = (uint64_t)randval % (i+1);
    tmp = a[j];
    a[j] = a[i];
    a[i] = tmp;
  }
  return 0;
}
static biclique *_InitMBC_Greedy(bipartite_graph *G)
{
  biclique *B = NULL;
  bitset V_cand = NULL;
  bitset test = NULL;
  uint64_t *U_cand = NULL;
  uint64_t i, u, v, U_cand_length, V_length;
  size_t M_ndata,N_ndata;
  int a,b;
  
  M_ndata = bitset_required_capacity(G->M);
  N_ndata = bitset_required_capacity(G->N);
  
  B = biclique_alloc(G->M, G->N);
  biclique_clear(B);

  U_cand_length = bipartite_graph_U_length(G);
  if(U_cand_length == 0) goto fail;
  
  V_length = bipartite_graph_V_length(G);
  if(V_length == 0) goto fail;

  U_cand = malloc(sizeof(uint64_t)*U_cand_length);
  if(U_cand == NULL) goto fail;

  V_cand = bitset_alloc(N_ndata);
  if(V_cand == NULL) goto fail;
  
  bipartite_graph_get_U_array(G, U_cand);
  if(_FisherYates(U_cand, U_cand_length)<0) goto fail;
  
  bipartite_graph_get_V_bitset(G, V_cand);

  test = bitset_alloc(N_ndata);
  if(test == NULL) goto fail;
  
  for(i=0;i<U_cand_length;i++) {
    TEST_LYU_STOP;
    
    u = U_cand[i];
    a = bitset_is_subset(B->V,G->U[u],N_ndata);
    // first test if u can be simply added without problem, if not skip
    if(!a) continue;
    
    biclique_add_U(B,u);

    bitset_intersect(test,V_cand,G->U[u],N_ndata);
    // attempt to add v
    while((a=bitset_random_pop(test, &v, N_ndata, &_internal_random_state, 0))>0) {
      a = bitset_is_member(G->V[v],u,M_ndata);
      b = bitset_is_subset(B->U,G->V[v],M_ndata);
      if(a && b) {
	// ok to add
	biclique_add_U(B,u);
	biclique_add_V(B,v);
	bitset_remove(V_cand,v,N_ndata);
	break;
      }
    }
    if(a < 0) goto fail;
  }

  free(U_cand);
  bitset_free(V_cand);
  bitset_free(test);
  
  return B;

 fail:
  if(U_cand) free(U_cand);
  if(V_cand) bitset_free(V_cand);
  if(test) bitset_free(test);
  if(B) biclique_free(B);
  return NULL;
}

biclique *InitMBC_Greedy(bipartite_graph *G, int init_iter)
{
  biclique *B = NULL, *NB = NULL;
  uint64_t current_size = 0, new_size;
  int i;

  for(i=0;i<init_iter;i++) {
    TEST_LYU_STOP;

    NB = _InitMBC_Greedy(G);
    if(NB == NULL) goto fail;
    new_size = biclique_size(NB);
    if(new_size > current_size) {
      if(B) biclique_free(B);
      B = NB;
      NB = NULL;
      current_size = new_size;
    }
  }

  return B;
  
 fail:
  if(B) biclique_free(B);
  if(NB) biclique_free(NB);
  return NULL;
 
}

biclique *InitMBC_Prune(bipartite_graph *G)
{
  uint64_t u,v,u_min,v_min;
  uint64_t deg, min_deg;
  int found;
  biclique *B = NULL;
  uint64_t length;

  G = bipartite_graph_dup(G);
  if(G == NULL) return NULL;
  
  for(;;) {
    if(bipartite_graph_is_complete(G)) break;
    
    // look for minimum degree U, and remove it
    min_deg = G->N;
    found = 0;
    for(u=0;u<G->M;u++) {
      if(G->U[u] == NULL) continue;
      deg = bipartite_graph_deg_u(G,u);
      if(deg < min_deg) {
	min_deg = deg;
	u_min = u;
	found = 1;
      }
    }
    
    if(found) {
      bipartite_graph_remove_u(G,u_min);
      if(bipartite_graph_is_complete(G)) break;
    }
    
    // look for minimum degree V, and remove it
    min_deg = G->M;
    found = 0;
    for(v=0;v<G->N;v++) {
      if(G->V[v] == NULL) continue;
      deg = bipartite_graph_deg_v(G,v);
      if(deg < min_deg) {
	min_deg = deg;
	v_min = v;
	found = 1;
      }
    }
    if(found) {
      bipartite_graph_remove_v(G,v_min);
      if(bipartite_graph_is_complete(G)) break;
    }
  }

  B = biclique_alloc(G->M, G->N);
  length = bipartite_graph_get_U_bitset(G, B->U);
  if(length == 0) goto fail;
  length = bipartite_graph_get_V_bitset(G, B->V);
  if(length == 0) goto fail;

  bipartite_graph_free(G);
  return B;

 fail:
  bipartite_graph_free(G);
  if(B) biclique_free(B);
  return NULL;  
}
biclique *InitMBC_Best(bipartite_graph *G, int init_iter)
{
  uint64_t C_size, NC_size;
  biclique *C = NULL, *NC = NULL;

#define TEST_AND_UPDATE				\
  do {						\
    if(NC != NULL) {				\
      NC_size = biclique_size(NC);		\
      if(NC_size > C_size) {			\
	if(C) biclique_free(C);			\
	C = NC;					\
	NC = NULL;				\
	C_size = NC_size;			\
      }						\
    }						\
  } while(0)				

  C_size = 0;

  NC = InitMBC_Star(G, 0);
  TEST_AND_UPDATE;

  NC = InitMBC_Star(G, 1);
  TEST_AND_UPDATE;

  // no need to try if even no star was found
  if(C == NULL) return NULL; 
  
  NC = InitMBC_Prune(G);
  TEST_AND_UPDATE;

  NC = InitMBC_Greedy(G, init_iter);
  TEST_AND_UPDATE;

#undef TEST_AND_UPDATE
  return C;
}

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bitset.h"

#ifndef NDEBUG
#include <stdio.h>
#endif

#define divmod(d,r,x,y) d = (x)/(y); r = (x)%(y)

size_t bitset_required_capacity(uint64_t maxlen)
{
  size_t ndata;
    
  ndata = (size_t)((maxlen / BITWISE_SET_UNIT_BITS) + !!(maxlen % BITWISE_SET_UNIT_BITS));
  if(ndata == 0) ndata = 1;
  assert(ndata > 0);
  assert(ndata*BITWISE_SET_UNIT_BITS >= maxlen);
  return ndata;
}

bitset bitset_alloc(size_t ndata)
{  
  return malloc(sizeof(uint64_t)*ndata);
}
void bitset_free(bitset S)
{
  if(S) free(S);
}
void bitset_clear(bitset S, size_t ndata)
{
  memset(S,0,sizeof(uint64_t)*ndata);
}
void bitset_setall(bitset S, size_t ndata)
{
  memset(S,0xFF,sizeof(uint64_t)*ndata);
}
void bitset_copy(bitset D, bitset S, size_t ndata)
{
  if(S == NULL) {
    bitset_clear(D,ndata);
  }
  else {
    memcpy(D,S,BITWISE_SET_UNIT_SIZE*ndata);
  }
}
uint64_t bitset_length(bitset S, size_t ndata)
{
  size_t i;
  uint64_t c = 0ULL;

  if(S == NULL) return 0;
  else {
    for(i=0;i<ndata;i++) c += __builtin_popcountll(S[i]);
    return c;
  }
}

int bitset_is_member(bitset S, uint64_t n, size_t ndata)
{
  size_t d,r;
  if(S == NULL) return 0;
  else {
    divmod(d,r,n,BITWISE_SET_UNIT_BITS);
    assert(d < ndata);
    return !!(S[d] & (1UL << r));
  }
}
int bitset_is_subset(bitset A, bitset B, size_t ndata)
{
  size_t i;
  uint64_t test;

  if(A == NULL) return 1;
  else if(B == NULL) return 0;
  else {
    for(i=0;i<ndata;i++) {
      test = A[i] & B[i];
      if(test != A[i]) return 0; // false
    }
    return 1; // true
  }
}
void bitset_add(bitset S, uint64_t n, size_t ndata)
{
  size_t d,r;

  divmod(d,r,n,BITWISE_SET_UNIT_BITS);
  assert(d < ndata);
  S[d] |= (1UL << r);
}
void bitset_remove(bitset S, uint64_t n, size_t ndata)
{
  size_t d,r;

  divmod(d,r,n,BITWISE_SET_UNIT_BITS);
  assert(d < ndata);
  S[d] &= ~(1UL << r);
}
void bitset_union(bitset O, bitset A, bitset B, size_t ndata)
{
  size_t i;

  if(A == NULL && B == NULL) bitset_clear(O,ndata);
  else if(A != NULL && B == NULL) bitset_copy(O,A,ndata);
  else if(A == NULL && B != NULL) bitset_copy(O,B,ndata);
  else {
    for(i=0;i<ndata;i++) O[i] = A[i] | B[i];
  }
}
void bitset_intersect(bitset O, bitset A, bitset B, size_t ndata)
{
  size_t i;

  if(A == NULL || B == NULL) bitset_clear(O,ndata);
  else {
    for(i=0;i<ndata;i++) O[i] = A[i] & B[i];
  }
}
void bitset_difference(bitset O, bitset A, bitset B, size_t ndata)
{
  size_t i;

  if(A == NULL) bitset_clear(O,ndata);
  else if(B == NULL) {} // nothing to do
  else {
    for(i=0;i<ndata;i++) O[i] = A[i] & (~(B[i]));
  }
}

int bitset_from_C_uint64_array(uint64_t *values,size_t N, bitset *S, size_t ndata, size_t *ndata_out)
{
  size_t i;
  uint64_t max_val = 0;
  uint64_t val;
  int c;

  if(ndata == 0) {
    for(i=0;i<N;i++) {
      val = values[i];
      if(val > max_val) max_val = val;
    }
    max_val += 1;
    ndata = bitset_required_capacity(max_val);
  }
  if(ndata_out != NULL) {
    *ndata_out = ndata;
  }
  
  *S = bitset_alloc(ndata);
  if(*S == NULL) return -1;
  
  bitset_clear(*S,ndata);

  c = 0;
  for(i=0;i<N;i++) {
    bitset_add(*S, values[i], ndata);
    c += 1;
  }
  return c;
}

void bitset_iterator_begin(bitset_iterator *it, bitset S, size_t ndata)
{
  it->S = S;
  it->ndata = ndata;
  it->i = 0;
  it->j = 0ULL;
  it->offset = 0ULL;
  it->val = 0ULL;
  it->current = 0ULL;
}

int bitset_iterator_next(bitset_iterator *it)
{
  if(it->S == NULL) return 0; // treat as empty set

  for(;;) {
    while(it->val == 0ULL) {
      if(it->i >= it->ndata) return 0;
      it->val = it->S[it->i];
      it->offset = BITWISE_SET_UNIT_BITS * it->i;
      it->i++;
      it->j=0;
    }
    
    while(it->val != 0UL && it->j < (uint64_t)BITWISE_SET_UNIT_BITS) {
      if(it->val & 1) {
	it->current = it->j + it->offset;
	it->val >>= 1;
	it->j++;
	return 1;
      }
      else {
	it->val >>= 1;
	it->j++;
      }
    }
  }
  assert(0); // unreachable
  return 0;
}

int bitset_pop(bitset S, uint64_t *n, size_t ndata, uint64_t pos, int keep)
{
  uint64_t i,j;
  uint64_t val,test;
  uint64_t c;
  
  if(S == NULL) return 0; // treat as empty set

  c = 0ULL;

  for(i=0;i<ndata;i++) {
    val = S[i];
    test = 1ULL;
    for(j=0;j<(uint64_t)BITWISE_SET_UNIT_BITS;j++, test <<= 1) {
      if(test & val) {
	if(c == pos) {
	  *n = j + i*BITWISE_SET_UNIT_BITS;
	  if(!keep) S[i] &= ~test;
	  return 1;
	}
	else c++;
      }
    }
  }
  return 0;
}

int bitset_random_pop(bitset S, uint64_t *n, size_t ndata, struct random_data *buf, int keep)
{
  uint64_t length;
  int32_t randval;
  uint64_t pos;
  
  if(S == NULL) return 0; // treat as empty set
  
  length = bitset_length(S,ndata);
  if(length == 0) return 0;

  if(length > RAND_MAX) return -1;

  if(random_r(buf, &randval)<0) return -1;

  pos = ((uint64_t)randval)%length;
  return bitset_pop(S,n,ndata,pos,keep);
}

void bitset_to_string(bitset S, char *buf, uint64_t N)
{
  size_t i;
  uint64_t j;
  uint64_t val;
  uint64_t offset;
  size_t ndata;
  uint64_t c;

  if(S == NULL) { // treat as empty set
    memset(buf, '0', N);
    buf[N+1] = '\0';
    return;
  }
  c = 0;
  ndata = bitset_required_capacity(N);
  for(i=0, offset=0;i<ndata && c<N;i++, offset += BITWISE_SET_UNIT_BITS) {
    val = S[i];
    for(j=0;j<(uint64_t)BITWISE_SET_UNIT_BITS && c<N;j++, val >>= 1UL) {
      *buf++ = (val & 1) ? '1' : '0';
      c++;
    }
  }
  *buf = '\0';
}

#ifndef NDEBUG
static void _bitset_print(bitset S, size_t ndata)
{
  int c;
  bitset_iterator it;
  
  bitset_iterator_begin(&it, S, ndata);
  printf("ndata=%lu:{", ndata);
  c = 0;
  while(bitset_iterator_next(&it)) {
    printf("%lu,",it.current);
    c++;
  }
  printf("} => %d\n", c);
}
void bitset_test()
{
  bitset S;
  uint64_t maxlen = 100;
  size_t ndata;
  uint64_t array[3] = {1,6,100};
  uint64_t n;
  
  ndata = bitset_required_capacity(maxlen);
  S = bitset_alloc(ndata);
  bitset_clear(S, ndata);
  
  bitset_add(S,11,ndata);
  bitset_add(S,63,ndata);
  bitset_add(S,13,ndata);
  bitset_add(S,99,ndata);
  bitset_add(S,5,ndata);
  bitset_add(S,7,ndata);
  bitset_add(S,127,ndata);
  _bitset_print(S,ndata);

  while(bitset_pop(S,&n,ndata,0,0)) {
    printf("pop %lu\n", n);
    _bitset_print(S,ndata);
  }
  printf("length=%lu\n", bitset_length(S,ndata));
  bitset_remove(S,10,ndata);
  bitset_remove(S,11,ndata);
  bitset_remove(S,64,ndata);
  _bitset_print(S,ndata);
  printf("length=%lu\n", bitset_length(S,ndata));
  bitset_free(S);
  
  bitset_from_C_uint64_array(array,3,&S,0,&ndata);
  _bitset_print(S,ndata);
  bitset_free(S);
}
#endif

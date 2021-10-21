#pragma once

#include <stdint.h>
#include <stdlib.h>

// assume uint64_t
typedef uint64_t *bitset;

typedef struct {
  bitset S;
  size_t ndata;
  size_t i;
  uint64_t j;
  uint64_t offset;
  uint64_t val;
  uint64_t current;
} bitset_iterator;

// number of bits in each unit of data
#define BITWISE_SET_UNIT_SIZE 8 // sizeof(uint64_t)
#define BITWISE_SET_UNIT_BITS 64

size_t bitset_required_capacity(uint64_t maxlen);

bitset bitset_alloc(size_t ndata);
void bitset_free(bitset S);
void bitset_setall(bitset S, size_t ndata);
void bitset_clear(bitset S, size_t ndata);

// convention: NULL treats as empty set for all readonly operation
void bitset_copy(bitset D, bitset S, size_t ndata);
uint64_t bitset_length(bitset S, size_t ndata);

int bitset_is_member(bitset S, uint64_t n, size_t ndata);
void bitset_add(bitset S, uint64_t n, size_t ndata);
void bitset_remove(bitset S, uint64_t n, size_t ndata);

int bitset_is_subset(bitset A, bitset B, size_t ndata);
// O = A \union B
void bitset_union(bitset O, bitset A, bitset B, size_t ndata);
// O = A \intersects B
void bitset_intersect(bitset O, bitset A, bitset B, size_t ndata);
// O = A \ B
void bitset_difference(bitset O, bitset A, bitset B, size_t ndata);

int bitset_from_C_uint64_array(uint64_t *values,size_t N, bitset *S, size_t ndata, size_t *ndata_out);

void bitset_iterator_begin(bitset_iterator *it, bitset S, size_t ndata);
int bitset_iterator_next(bitset_iterator *it); // return 0 if stop
int bitset_pop(bitset S, uint64_t *n, size_t ndata, uint64_t pos, int keep);
int bitset_random_pop(bitset S, uint64_t *n, size_t ndata, struct random_data *buf, int keep);

void bitset_to_string(bitset S, char *buf, uint64_t N);

#ifndef NDEBUG
void bitset_test();
#endif

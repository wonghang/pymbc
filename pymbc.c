#include <Python.h>
#include <numpy/arrayobject.h>
#include "structmember.h"
#include <stdbool.h>
#include <signal.h>

#ifndef NDEBUG
#include <stdio.h>
#endif

#include "bitset.h"
#include "bipartite_graph.h"
#include "lyu.h"

#if (PY_VERSION_HEX < 0x03000000)
#error "require Python 3.0"
#endif

extern int lyu_stop;

static void my_sigint_handler(int signum)
{
  lyu_stop = 1;
}

#define SIGINT_TRAP sighandler_t old_sigint_handler = NULL; old_sigint_handler = signal(SIGINT, my_sigint_handler)
#define SIGINT_RECOVER if(old_sigint_handler) signal(SIGINT, old_sigint_handler)

static bitset python_numpy_uint64_to_bitset(PyObject *obj, size_t ndata)
{
  bitset S = NULL;
  void *ptr = NULL;
  npy_intp shape[1];
  
  if(obj == Py_None) {
    S = bitset_alloc(ndata);
    if(S == NULL) return NULL;
    bitset_setall(S,ndata);
  }
  else {
    if(PyArray_AsCArray(&obj,
			(void *)&ptr,
			shape,
			1,
			PyArray_DescrFromType(NPY_UINT64)) < 0) {
      goto fail;
    }
  
    if(bitset_from_C_uint64_array(ptr, shape[0], &S, ndata, NULL) < 0) {
      goto fail;
    }
  }

  return S;

 fail:
  if(S) bitset_free(S);
  if(ptr) PyArray_Free(obj, ptr);
  return NULL;
}

static PyObject *bitset_to_python_list(bitset S, size_t ndata)
{
  PyObject *ret = NULL;
  bitset_iterator it;
  Py_ssize_t idx;

  ret = PyList_New(bitset_length(S, ndata));
  if(ret == NULL) goto fail;

  bitset_iterator_begin(&it, S, ndata);
  idx = 0;
  while(bitset_iterator_next(&it)) {
    PyList_SetItem(ret, idx++, PyLong_FromUnsignedLongLong(it.current));
  }
  return ret;

 fail:
  Py_XDECREF(ret);
  return NULL;
}
/*
  (U,V,D, tau_u, tau_v, init_type, init_iter, seed, use_star, star_max_iter, optimize)
 */
static PyObject *maximum_biclique_search(PyObject *self, PyObject *args)
{
  PyObject *objU, *objV, *objD;
  uint64_t tau_u, tau_v;
  unsigned int seed;
  int init_type, init_iter;
  int star, star_max_iter;
  unsigned int optimize;
  
  npy_intp shape[2];
  uint64_t M, N;
  size_t M_ndata, N_ndata;
  int8_t **Dptr = NULL;
  bitset U = NULL, V = NULL;
  biclique *C = NULL;
  bipartite_graph *G = NULL;
  int ret_c;

  PyObject *ret, *retU = NULL, *retV = NULL;

  SIGINT_TRAP;

  if(!PyArg_ParseTuple(args, "OOOKKiiIpiI",
		       &objU, &objV, &objD,
		       &tau_u, &tau_v,
		       &init_type, &init_iter,
		       &seed,
		       &star, &star_max_iter,
		       &optimize)) {
    PyErr_SetString(PyExc_TypeError, "Unable to parse arguments");
    goto fail;
  }

  if(init_iter <= 0) {
    PyErr_SetString(PyExc_ValueError, "init_iter should be > 0");
    goto fail;
  }

  if(star_max_iter <= 0) {
    PyErr_SetString(PyExc_ValueError, "star_max_iter should be > 0");
    goto fail;
  }

  if(PyArray_AsCArray(&objD,
		      (void *)&Dptr,
		      shape,
		      2,
		      PyArray_DescrFromType(NPY_BOOL)) < 0) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert D");
    goto fail;
  }
  M = (uint64_t)shape[0];
  M_ndata = bitset_required_capacity(M);
  N = (uint64_t)shape[1];
  N_ndata = bitset_required_capacity(N);

  U = python_numpy_uint64_to_bitset(objU,M_ndata);
  if(U == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert U");
    goto fail;
  }

  V = python_numpy_uint64_to_bitset(objV,N_ndata);
  if(V == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert V");
    goto fail;
  }

  G = bipartite_graph_from_C_bool_array(M,N,U,V,Dptr);
  if(!G) {
    PyErr_SetString(PyExc_MemoryError, "Failed to allocate G");
    goto fail;
  }

  if(lyu_random_init(seed)<0) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to setup random generator");
    goto fail;
  }

  Py_BEGIN_ALLOW_THREADS;
  lyu_stop = 0;

  switch(init_type) {
  case 0: 
    C = InitMBC_Best(G, init_iter);
    break;
  case 1:
    C = InitMBC_Greedy(G, init_iter);
    break;
  case 2:
    C = InitMBC_Prune(G);
    break;
  case 3:
    C = InitMBC_Star(G, 0);
    break;
  case 4:
    C = InitMBC_Star(G, 1);
    break;
  default:
    C = NULL;
    break;
  }

  Py_END_ALLOW_THREADS;

  if(lyu_stop) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "SIGINT received inside pymbc");
    goto fail;
  }
  
  if(!C) {
    PyErr_SetString(PyExc_ValueError, "Failed to InitMBC, try another approach");
    goto fail;
  }

#ifndef NDEBUG
  //  bipartite_graph_dump(G);
  /* printf("initial size (%lu,%lu) => %lu\n", */
  /* 	 biclique_U_length(C), */
  /* 	 biclique_V_length(C), */
  /* 	 biclique_size(C)); */
  if(!biclique_is_biclique(G,C)) {
    PyErr_SetString(PyExc_RuntimeError, "BUG #1");
    goto fail;
  }
#endif
  
  Py_BEGIN_ALLOW_THREADS;
  lyu_stop = 0;

  if(star) {
    ret_c = MBC_star(G, C, tau_u, tau_v, star_max_iter, optimize);
  }
  else {
    ret_c = MBC(G, C, tau_u, tau_v);
  }

  Py_END_ALLOW_THREADS;

  if(lyu_stop) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "SIGINT received inside pymbc");
    goto fail;
  }

  if(ret_c < 0) {
    PyErr_SetString(PyExc_RuntimeError, "MBC return error");
    goto fail;
  }

#ifndef NDEBUG
  if(!biclique_is_biclique(G,C)) {
    PyErr_SetString(PyExc_RuntimeError, "BUG #2");
    goto fail;
  }
#endif

  retU = bitset_to_python_list(C->U, M_ndata);
  if(retU == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to get back U");
    goto fail;
  }
  retV = bitset_to_python_list(C->V, N_ndata);
  if(retV == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to get back V");
    goto fail;
  }

  bitset_free(U);
  bitset_free(V);
  PyArray_Free(objD, Dptr);
  bipartite_graph_free(G);
  biclique_free(C);
  
  ret = Py_BuildValue("(OO)", retU, retV);
  Py_DECREF(retU);
  Py_DECREF(retV);
  SIGINT_RECOVER;
  return ret;

 fail:
  Py_XDECREF(retU);
  Py_XDECREF(retV);
  if(U) bitset_free(U);
  if(V) bitset_free(V);
  if(Dptr) PyArray_Free(objD, Dptr);
  if(G) bipartite_graph_free(G);
  if(C) biclique_free(C);
  SIGINT_RECOVER;
  return NULL;
}

static const char * const maximum_biclique_search_docstring =
  "maximum_biclique_search(U: Optional[numpy.typing.NDArray[np.uint64]],\n"
  "                        V: Optional[numpy.typing.NDArray[np.uint64]],\n"
  "                        D: numpy.typingNDArray[np.bool],\n"
  "                        tau_u: int,\n"
  "                        tau_v: int,\n"
  "                        init_type: int,\n"
  "                        init_iter: int,\n"
  "                        seed: int,\n"
  "                        ues_star: bool,\n"
  "                        star_max_iter: int,\n"
  "                        optimize: int) -> Tuple[List,List]\n"
  "\n"
  "Perform a maximum biclique search by using the algorithm MBC or MBC* on the paper.\n"
  "Algorithm of InitMBC can be selected from 5 different options.\n"
  "\n"
  "Parameters:\n"
  "U             -- array of integers (starts from 0) indicate the vertices of one side of a bipartite graph (e.g. [0,1,2,4])) or None\n"
  "V             -- array of integers (starts from 0) indicate the vertices of another side of a bipartite graph (e.g. [0,3,5,6]) or None\n"
  "D             -- Boolean matrix with entries (i,j)=True if (U_i,V_j) is connected with an edge, False otherwise\n"
  "                 **The code determines the size of the graph based on the shape of D. U and V should not contain entries greater than the size**\n"
  "tau_u         -- tau_u as on the paper (minimum number of U to search)\n"
  "tau_v         -- tau_v as on the paper (minimum number of V to search)\n"
  "init_type     -- Algorithm for InitMBC:\n"
  "                 0 - Best (try all below and use the one with maximum size)\n"
  "                 1 - Greedy (start from an empty graph and add vertex until it is not possible)\n"
  "                 2 - Prune (start from everything and remove vertex until it finds a biclique)\n"
  "                 3 - Star (start from a star with maximum size, check from U, i.e. U has one vertex)\n"
  "                 4 - Star (start from a star with maximum size, check from V, i.e. V has one vertex)\n"
  "init_iter     -- Greedy algorithm is randomized and indicate how many iterations to run\n"
  "seed          -- Random seed for any randomness in the algorithm\n"
  "use_star      -- Whether to use MBC* on the paper, if use_star==False, the program will only use MBC\n"
  "star_max_iter -- the MAX_ITER variable on the paper, i.e. maximum of pruning step in MBC*\n"
  "optimize      -- Flags to indicate whether to use further optimization of Reduce2Hop written on the paper\n"
  "                 0 - Disable\n"
  "                 1 - Enable #1 (Early Pruning)\n"
  "                 2 - Enable #2 (Early Skipping)\n"
  "                 3 - Enable Both\n"
  "\n"
  "Return:\n"
  "U             -- List of the vertices of the found biclique on one side\n"
  "V             -- List of the vertices of the found biclique on another side\n"
  "\n"
  "Remarks:\n"
  "(1) len(U) >= tau_u and len(V) >= tau_v are not guaranteed. The code will return the one found by InitMBC if MBC/MBC* is unable to find a better biclique\n"
  "(2) The code raises a ValueError if no biclique is found during InitMBC. User should try different value of init_type. If init_type = 0 is already indicated, then there are no bicliques in the bipartite graph.\n"
  "(3) If you get a TypeError exception, try to recast all arguments according to the type annotations";

static PyObject *greedy_biclique(PyObject *self, PyObject *args)
{
  PyObject *objU, *objV, *objD;
  unsigned int seed;
  int iter;
  
  npy_intp shape[2];
  uint64_t M, N;
  size_t M_ndata, N_ndata;
  int8_t **Dptr = NULL;
  bitset U = NULL, V = NULL;
  biclique *C = NULL;
  bipartite_graph *G = NULL;

  PyObject *ret, *retU = NULL, *retV = NULL;

  SIGINT_TRAP;

  if(!PyArg_ParseTuple(args, "OOOIi",
		       &objU, &objV, &objD,
		       &seed, &iter)) {
    PyErr_SetString(PyExc_TypeError, "Unable to parse arguments");
    goto fail;
  }

  if(PyArray_AsCArray(&objD,
		      (void *)&Dptr,
		      shape,
		      2,
		      PyArray_DescrFromType(NPY_BOOL)) < 0) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert D");
    goto fail;
  }
  M = (uint64_t)shape[0];
  M_ndata = bitset_required_capacity(M);
  N = (uint64_t)shape[1];
  N_ndata = bitset_required_capacity(N);

  U = python_numpy_uint64_to_bitset(objU,M_ndata);
  if(U == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert U");
    goto fail;
  }

  V = python_numpy_uint64_to_bitset(objV,N_ndata);
  if(V == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to convert V");
    goto fail;
  }

  G = bipartite_graph_from_C_bool_array(M,N,U,V,Dptr);
  if(!G) {
    PyErr_SetString(PyExc_MemoryError, "Failed to allocate G");
    goto fail;
  }

  if(lyu_random_init(seed)<0) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to setup random generator");
    goto fail;
  }

  Py_BEGIN_ALLOW_THREADS;
  lyu_stop = 0;

  C = InitMBC_Greedy(G, iter);

  Py_END_ALLOW_THREADS;

  if(lyu_stop) {
    PyErr_SetString(PyExc_KeyboardInterrupt, "SIGINT received inside pymbc");
    goto fail;
  }
  
  if(!C) {
    PyErr_SetString(PyExc_ValueError, "Failed to run InitMBC_Greedy");
    goto fail;
  }
  
#ifndef NDEBUG
  if(!biclique_is_biclique(G,C)) {
    PyErr_SetString(PyExc_RuntimeError, "BUG #3");
    goto fail;
  }
#endif

  retU = bitset_to_python_list(C->U, M_ndata);
  if(retU == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to get back U");
    goto fail;
  }
  retV = bitset_to_python_list(C->V, N_ndata);
  if(retV == NULL) {
    PyErr_SetString(PyExc_TypeError, "Unable to get back V");
    goto fail;
  }

  bitset_free(U);
  bitset_free(V);
  PyArray_Free(objD, Dptr);
  bipartite_graph_free(G);
  biclique_free(C);

  ret = Py_BuildValue("(OO)", retU, retV);
  Py_DECREF(retU);
  Py_DECREF(retV);
  SIGINT_RECOVER;
  return ret;
  
 fail:
  Py_XDECREF(retU);
  Py_XDECREF(retV);
  if(U) bitset_free(U);
  if(V) bitset_free(V);
  if(Dptr) PyArray_Free(objD, Dptr);
  if(G) bipartite_graph_free(G);
  if(C) biclique_free(C);
  SIGINT_RECOVER;
  return NULL;
}

/*
  if(!PyArg_ParseTuple(args, "OOOIi",
		       &objU, &objV, &objD,
		       &seed, &iter)) {
*/
static const char * const greedy_biclique_docstring =
  "greedy_biclique(U: Optional[numpy.typing.NDArray[np.uint64]],\n"
  "                V: Optional[numpy.typing.NDArray[np.uint64]],\n"
  "                D: numpy.typingNDArray[np.bool],\n"
  "                seed: int,\n"
  "                iter: bool) -> Tuple[List,List]\n"
  "\n"
  "This method exposes an interface for the greedy algorithm used in InitMBC.\n"
  "\n"
  "Parameters:\n"
  "U             -- array of integers (starts from 0) indicate the vertices of one side of a bipartite graph (e.g. [0,1,2,4])) or None\n"
  "V             -- array of integers (starts from 0) indicate the vertices of another side of a bipartite graph (e.g. [0,3,5,6]) or None\n"
  "D             -- Boolean matrix with entries (i,j)=True if (U_i,V_j) is connected with an edge, False otherwise\n"
  "                 **The code determines the size of the graph based on the shape of D. U and V should not contain entries greater than the size**\n"
  "seed          -- Random seed for any randomness in the algorithm\n"
  "iter          -- This algorithm is randomized and indicate how many iterations to run\n"
  "\n"
  "Return:\n"
  "U             -- List of the vertices of the found biclique on one side\n"
  "V             -- List of the vertices of the found biclique on another side\n"
  "\n"
  "Remarks:\n"
  "(1) If you get a TypeError exception, try to recast all arguments according to the type annotations";

static PyMethodDef pymbcmethods[] =
  {
   { "maximum_biclique_search", (PyCFunction)maximum_biclique_search, METH_VARARGS, maximum_biclique_search_docstring },
   { "greedy_biclique", (PyCFunction)greedy_biclique, METH_VARARGS, greedy_biclique_docstring },
   { NULL, NULL, 0, NULL }
  };

static const char * const pymbc_docstring = 
  "Maximum Biclique Search\n"
  "\n"
  "An implementation of\n"
  "Lyu, Bingqing, et al. \"Maximum biclique search at billion scale.\" Proceedings of the VLDB Endowment (2020).\n"
  "\n";

static PyModuleDef pymbcmodule =
  {
   PyModuleDef_HEAD_INIT,
   "pymbc",
   pymbc_docstring,
   -1,
   pymbcmethods, NULL, NULL, NULL, NULL
  };

PyMODINIT_FUNC PyInit_pymbc(void)
{
  import_array();
  return PyModule_Create(&pymbcmodule);
}

#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <Python.h>
#include "spkmeans.c"

/* algorithm */
double** calc_means_impl(double **vectors, double **cntrd, int v_len);
double*** calc_means(double ***clusters, double **cntrd, int *clstr_lens); 
int update_centroids(double ***clusters, double **cntrd, int *clstr_lens);
int calc_new_centroid(double **cluster, double **cntrds, int len, int c);
double*** classify_vectors(double ***clusters, double **cntrd, int *clstr_lens);
double*** set_vectors_in_clusters(double **vectors, double **cntrd,int v_len, int *clstr_lens);
/* helpers */
void set_a_vector_in(double ***clusters, double *vector, int index, int* clstr_lens);
/* calculators */
int find_closest_centroid(double **cntrds, double *vector);
/* c- api*/
PyObject* to_py_lst(double **centroids, int kn, int dimn);
/* contructors*/
int* init_clusters_lens();
double* init_double_array(int len);
double*** init_clusters_array();
double** init_2d_double_array(int len, int inner_size);
/* free memory methods*/
void free_clusters(double ***clusters, int *clstr_lens);


int iter_limit=300,dim,k;
double epsilon=0;



/* algorithm */
double** calc_means_impl(double **vectors, double **cntrd, int v_len){
    double ***clusters;
    int *clusters_lens;
    clusters_lens = init_clusters_lens(); 
    clusters = set_vectors_in_clusters(vectors, cntrd, v_len, clusters_lens);
    clusters = calc_means(clusters, cntrd, clusters_lens);
    free_clusters(clusters, clusters_lens);
    free(clusters_lens);
    return cntrd;
}

double*** set_vectors_in_clusters(double **vectors, double **cntrds, int v_len, int *clstr_lens){
    double ***clusters;
    int i;
    int index = 0;
    clusters = init_clusters_array();
    for (i=0; i < v_len; i++){
        index = find_closest_centroid(cntrds, vectors[i]);
        set_a_vector_in(clusters, vectors[i], index, clstr_lens);
    }
    free(vectors);
    return clusters;
}

 double*** calc_means(double ***clusters, double **cntrd, int *clstr_lens){
    int iter = 1;
    int converg = 0;
    while (iter < iter_limit && converg < k  ){
        converg ++;
        clusters = classify_vectors(clusters, cntrd, clstr_lens);
        converg = update_centroids(clusters, cntrd, clstr_lens);
        iter++;
    } 
    return clusters;
 }
 
/* calculate & update the new centroids after clustering: return number of (delta < e) */
 int update_centroids(double ***clusters, double **cntrds, int *clstr_lens){
    int i;
    int converg = 0;
    for (i=0; i <k ; i++){
        converg += calc_new_centroid(clusters[i], cntrds, clstr_lens[i], i);
    }
    return converg;

 }


int calc_new_centroid(double **cluster, double **cntrds, int len, int c){
    double *vector = (double*) calloc(dim, sizeof(double)); 
    if (vector==NULL){ error_message(); }
    int i, j;
    double dist;
    int res = 0;
    for (i=0; i< len; i++){
        for (j=0; j< dim; j++){ 
            vector[j] += cluster[i][j];
        }
    }
    for (i=0; i<dim; i++){ vector[i] = vector[i]/len; }
    dist = calc_distance(vector, cntrds[c], dim); 
    if (dist < epsilon){ res++; }
    free(cntrds[c]);
    cntrds[c]=vector;

    return res;
}   

/* classify_vectors to closest cluster: return new clusters array */    
double*** classify_vectors(double ***clusters, double **cntrds, int *clstr_lens){ 
    int *new_clstr_lens;
    double ***new_clusters;
    int i, j; 
    int index=0;
    new_clstr_lens = init_clusters_lens();  
    new_clusters = init_clusters_array(); 
    for (i=0; i<k; i++){
        for (j=0; j < clstr_lens[i]; j++){
            index = find_closest_centroid(cntrds, clusters[i][j]);
            set_a_vector_in(new_clusters, clusters[i][j], index, new_clstr_lens);
        }
        free(clusters[i]); 
    }
    free(clusters);
    for (i=0;i<k;i++){ clstr_lens[i] = new_clstr_lens[i]; }
    free(new_clstr_lens);

    return new_clusters;    
} 




 /*helpers*/ 
 void set_a_vector_in(double ***clusters, double *vector, int index, int* clstr_lens){
    clusters[index] = (double**) realloc(clusters[index],((clstr_lens[index]+1)* sizeof(double*)));
    clusters[index][clstr_lens[index]] = vector; 
    clstr_lens[index]++;
}

int find_closest_centroid(double **cntrds, double *vector){
    double min, dist;
    int c;
    int index=0;
    min = INFINITY;
    for (c=0; c < k; c++){
        dist = calc_distance(vector, cntrds[c], dim);
        if (dist < min){
            min = dist;
            index = c;
        }
    }
    return index; 
}


/* arrays contructers */ 
int* init_clusters_lens(){
    int* arr;
    arr = (int*) calloc(k, sizeof(int)); 
    if (arr == NULL){ error_message();}
    return arr;
}

double*** init_clusters_array(){
    int i;
    double*** clusters;
    clusters = (double***) malloc( k* sizeof(double**)); 
    if (clusters == NULL){ error_message(); }
    for (i=0; i<k; i++){
        clusters[i] = (double**) malloc(sizeof(double*));
        if (clusters[i] == NULL) { error_message(); }
    }
    return clusters;
}

double* init_double_array(int len){
    double* array = (double*) malloc(len*sizeof(double));
    if (array == NULL){ error_message(); } 
    return array;
}

double** init_2d_double_array(int len, int inner_size){
    int i; double* inner_arr; double** array;
    array = (double**) malloc(len*sizeof(double**));
    if (array == NULL){ error_message(); } 
    for (i=0; i < len; i++){  
        inner_arr = init_double_array(inner_size);
        array[i] = inner_arr; 
    }
    return array;
}



void free_clusters(double ***clusters, int *clstr_lens){
    int i;
    for (i=0; i < k; i++){
        free_2d_double_array(clusters[i],clstr_lens[i]);
    }
    free(clusters);
}






PyObject* to_py_lst(double **centroids, int kn, int dimn){
    int i,c;
    PyObject* python_lst;
    PyObject* python_vector;
    PyObject* python_float;
    python_lst = PyList_New(kn);
    for (i = 0; i < kn; ++i){
        python_vector = PyList_New(dimn);
        for (c=0; c < dimn; c++){
            python_float = Py_BuildValue("d", centroids[i][c]);
            PyList_SetItem(python_vector, c, python_float);
        } 
        PyList_SetItem(python_lst, i, python_vector);
    }
    free_2d_double_array(centroids, kn);
    return python_lst;
}

static PyObject* spk(PyObject *self, PyObject *args){
    PyObject *vectors_arr;
    PyObject *centroids_indexes;
    PyObject *vector;
    int v_len;
/* This parses the Python arguments into a double** (O) variable named vectors and double* (O) variable named centroids */
    if(!PyArg_ParseTuple(args, "OOi", &vectors_arr, &centroids_indexes, &dim)) {  /* change to k !!!!!!! */
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    v_len =  PyObject_Length(vectors_arr);
    k = PyObject_Length(centroids_indexes);
    if (k < -1 || v_len < 0 || dim < 0 ) {
        return NULL;
    }
    double **cntrds;
    cntrds = init_2d_double_array(k, dim);
    int i, c;
    for (i=0; i < k; i++){
        vector = PyList_GetItem(vectors_arr,PyLong_AsLong(PyList_GetItem(centroids_indexes,i)));
        for (c=0; c < dim; c++){
            cntrds[i][c]= PyFloat_AsDouble(PyList_GetItem(vector,c));
        }
    }
    double **vectors;
    vectors = init_2d_double_array(v_len, dim);
    for (i=0; i < v_len; i++){
        for (c=0; c < dim; c++){
            vectors[i][c] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(vectors_arr,i),c));
        }
    }
    return to_py_lst(calc_means_impl(vectors, cntrds, v_len), k, dim);
}


static PyObject* wam(PyObject *self, PyObject *args){
    PyObject *vectors_arr;
    int v_len;
/* This parses the Python arguments into a double** (O) variable named vectors_arr */
    if(!PyArg_ParseTuple(args, "O", &vectors_arr)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    v_len =  PyObject_Length(vectors_arr);
    dim = PyObject_Length(PyList_GetItem(vectors_arr, 0));
    if (v_len < 0 || dim < 0 ) {
        return NULL;
    }   
    double **vectors;
    vectors = init_2d_double_array(v_len, dim);
    for (int i=0; i < v_len; i++){
        for (int c=0; c < dim; c++){
            vectors[i][c] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(vectors_arr,i),c));
        }
    }
    double **w = build_adjacency_matrix(vectors, v_len, dim);   
    free_2d_double_array(vectors, v_len);
    return to_py_lst(w, v_len, v_len);
}

static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *vectors_arr;
    int v_len;
/* This parses the Python arguments into a double** (O) variable named vectors_arr */
    if(!PyArg_ParseTuple(args, "O", &vectors_arr)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    v_len =  PyObject_Length(vectors_arr);
    dim = PyObject_Length(PyList_GetItem(vectors_arr, 0));
    if (v_len < 0 || dim < 0 ) {
        return NULL;
    }   
    double **vectors;
    vectors = init_2d_double_array(v_len, dim);
    for (int i=0; i < v_len; i++){
        for (int c=0; c < dim; c++){
            vectors[i][c] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(vectors_arr,i),c));
        }
    }
    double **w = build_adjacency_matrix(vectors, v_len, dim);
    double **d = build_dd_matrix(v_len, w);  
    free_2d_double_array(w, v_len);
    free_2d_double_array(vectors, v_len);
    return to_py_lst(d, v_len, v_len);
}


static PyObject* gl(PyObject *self, PyObject *args){
    PyObject *vectors_arr;
    int v_len;
/* This parses the Python arguments into a double** (O) variable named vectors_arr */
    if(!PyArg_ParseTuple(args, "O", &vectors_arr)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    v_len =  PyObject_Length(vectors_arr);
    dim = PyObject_Length(PyList_GetItem(vectors_arr, 0));
    if (v_len < 0 || dim < 0 ) {
        return NULL;
    }   
    double **vectors;
    vectors = init_2d_double_array(v_len, dim);
    for (int i=0; i < v_len; i++){
        for (int c=0; c < dim; c++){
            vectors[i][c] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(vectors_arr,i),c));
        }
    }
    double **w = build_adjacency_matrix(vectors, v_len, dim);
    double **d = build_dd_matrix(v_len, w);
    build_lm_matrix(v_len, w, d);  
    free_2d_double_array(w, v_len);
    free_2d_double_array(vectors, v_len);
    return to_py_lst(d, v_len, v_len);
}

static PyObject* jacobi(PyObject *self, PyObject *args){/*not completed*/
    PyObject *vectors_arr;
    int v_len;
/* This parses the Python arguments into a double** (O) variable named vectors_arr */
    if(!PyArg_ParseTuple(args, "O", &vectors_arr)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    v_len =  PyObject_Length(vectors_arr);
    dim = PyObject_Length(PyList_GetItem(vectors_arr, 0));
    if (v_len < 0) {
        return NULL;
    }   
    double **vectors;
    vectors = init_2d_double_array(v_len, dim);
    for (int i=0; i < v_len; i++){
        for (int c=0; c < dim; c++){
            vectors[i][c] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(vectors_arr,i),c));
        }
    }
    double *eigenvals = jacobi_algo(vectors, v_len);  
    double **eigen = init_2d_double_array(v_len+1,v_len);
    for(int i=0;i<v_len;i++){
        for(int j=0;j<v_len;j++){
            eigen[i][j] = vectors[i][j];
        }
        eigen[v_len][i] = eigenvals[i];
    }    
    free(eigenvals);
    free_2d_double_array(vectors, v_len);
    return to_py_lst(eigen, v_len+1, v_len);
}


static PyMethodDef kmeansMethods[] = {
    {"gl",                   /* the Python method name that will be used */
      (PyCFunction) gl, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("Calculate the laplacian graph accordingly for the arguments") /*  The docstring for the function */
    },
    {"ddg",                   /* the Python method name that will be used */
      (PyCFunction) ddg, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("Calculate the degree graph accordingly for the arguments") /*  The docstring for the function */
    },
    {"wam",                   /* the Python method name that will be used */
      (PyCFunction) wam, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("Calculate the weight graph accordingly for the arguments") /*  The docstring for the function */
    },
    {"jacobi",                   /* the Python method name that will be used */
      (PyCFunction) jacobi, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("Calculate the eigenvalues and eigenvectors accordingly for the arguments") /*  The docstring for the function */
    },
    {"spk",                   /* the Python method name that will be used */
      (PyCFunction) spk, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("Run spkmean clustering algorithem accordingly for the arguments") /*  The docstring for the function */
    },
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    kmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};


PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}






#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"


double** init_vectors_array( FILE *f);
int calc_vector_dim(char *line1);
double* convert_vector_string_to_double(char *line);
void print_centroids_array(double **cntrd, int k, int dim);
void print_ev(double *ev_arr);
void free_2d_double_array(double **arr, int len);
const int Max_iter = 100;
const double Epsilon = 0.00001;
int dim, n;

double** build_adjacency_matrix(double **vectors, int n, int dim){
    int i,j;
    double dist, weight;
    double** res = (double**) malloc(sizeof(double*)*n);
    if (res==NULL) { error_message(); }
    for (i=0; i<n;i++){
        res[i]= (double*) malloc(sizeof(double)*n);
        if (res[i]==NULL) { error_message(); }
    }
    for (i=0; i<n;i++){
        res[i][i]=0;
        for (j=i+1; j<n;j++) {
            dist = calc_distance(vectors[i], vectors[j], dim);
            weight = exp(-(pow(dist, 2) / 2));
            res[i][j] = weight;
            res[j][i] = weight;
        }
    }
    return res;
}

double** build_dd_matrix(int n, double **am){
    int i,j;
    double d;
    double** res = (double**) malloc(sizeof(double*)*n);
    if (res==NULL) { error_message(); }
    for (i=0; i<n;i++){
        res[i] = (double*) calloc(n,sizeof(double));
        if (res[i]==NULL) { error_message(); }
        d=0;
        for (j=0; j<n;j++){
            d+=am[i][j];
        }
        res[i][i]=d;
    }
    return res;
}

void build_lm_matrix(int n, double **am, double **ddm){
    int i,j;
    for (i=0; i<n;i++){
        for (j=i+1; j<n;j++){
            ddm[i][j] = am[i][j]*-1;
            ddm[j][i] = am[i][j]*-1;
        }
    }
}

double* jacobi_algo(double **lm, int n){
    double *res; int *max_entries;
    double **pm, **p1;
    int i,j;
    double of_delta,of_A,of_B;
    int iter = 0;
    p1 = (double**) malloc(sizeof(double*)*n);
    if (p1==NULL) { error_message(); }
    for (i=0;i<n;i++){
        p1[i] = (double*) calloc(n,sizeof(double));
        if (p1[i]==NULL) { error_message(); }
        p1[i][i] = 1.0;
    }
    of_delta = 1;
    while (iter < Max_iter && of_delta > Epsilon){
        of_A = of_f(lm, n);
        max_entries = find_maximal_entry(lm, n);
        pm = build_p_matrix(lm, n, max_entries);
        if (pm==NULL) { error_message(); }
        mult_ptap(lm, pm,n,max_entries);
        of_B = of_f(lm, n);
        of_delta = of_A-of_B;
        of_B = of_A;
        p1 = matrix_mult(p1, pm, n);
        free(max_entries);
        iter++;
        
    }
    res = (double*) malloc(sizeof(double)*n);
    if (res==NULL) { error_message(); }
    for (i=0; i<n; i++){
        res[i] = lm[i][i];
    }
    for (i=0; i<n;i++){
        for (j=0; j<n; j++){
            lm[i][j] = p1[i][j];
        }
    }
    free_2d_double_array(p1,n);
    return res;
}

int* find_maximal_entry(double **lm, int n){
    int i, j;
    double max = 0;
    int *res = (int*) calloc(2,sizeof(int));
    if (res==NULL){ error_message();}
    for (i=0; i<n;i++){
        for (j=i+1; j<n;j++){
            if (fabs(lm[i][j])> max){
                res[0]=j;
                res[1]=i;
                max=fabs(lm[i][j]);
            }
        }
    }
    return res;
}

double** build_p_matrix(double **lm,int n, int *pos){
    int i,j; double **pm;
    double theta,c,t,s,sign;
    pm = (double**) malloc(sizeof(double*)*n);
    if (pm==NULL) { error_message(); }

    for (i=0;i<n;i++){
        pm[i] = (double*) calloc(n,sizeof(double));
        if (pm[i]==NULL) { error_message(); }
        pm[i][i]=1;
    }
    theta = (lm[pos[1]][pos[1]]-lm[pos[0]][pos[0]]) / (2*lm[pos[0]][pos[1]]);
    sign = theta < 0 ? -1.0 : 1.0;
    t = sign/(fabs(theta)+ sqrt(pow(theta,2)+1));
    c = 1/sqrt(pow(t,2)+1);
    s = t*c;
    i= pos[0]<pos[1] ? pos[0]:pos[1] ;
    j = i==pos[0] ? pos[1]:pos[0];
    pm[i][i]=c;
    pm[j][j]=c;
    pm[i][j]=-s;
    pm[j][i]=s;

    return pm;
}

void mult_ptap(double **am, double **pm,int n, int *pos){
    int r,i,j;
    double c,s;
    double **tm = copy_matrix(am,n);
    if (tm==NULL) { error_message(); }
    i = pos[0] < pos[1] ? pos[0] : pos[1];
    j = i==pos[0] ? pos[1] : pos[0];
    c = pm[i][i];
    s = pm[i][j];

    for (r=0; r<n; r++){
        if (r!=i && r!=j){
            am[r][i] = c*tm[r][i]-s*tm[r][j];
            am[r][j] = c*tm[r][j]+s*tm[r][i];
            am[i][r] = am[r][i];
            am[j][r] = am[r][j];
        }
    }
    am[i][i]= pow(c,2)*(tm[i][i])+pow(s,2)*(tm[j][j])-2*s*c*tm[i][j];
    am[j][j]= pow(s,2)*(tm[i][i])+pow(c,2)*(tm[j][j])+2*s*c*tm[i][j];
    am[i][j]=0;
    am[j][i]=0;
    free_2d_double_array(tm,n);
}

double** matrix_mult(double **p1, double **p2, int n){
    int i,j,k;
    double **res = (double**) malloc(sizeof(double*)*n);
    if (res==NULL){ error_message(); }
    for (i=0;i<n;i++){
        res[i] = (double*) calloc(n,sizeof(double));
        if (res[i]==NULL){ error_message(); }
        for (j=0; j<n; j++){
            for (k=0; k<n; k++){
                res[i][j]+=p1[i][k]*p2[k][j];
            }
        }
    }
    free_2d_double_array(p1,n);
    free_2d_double_array(p2,n);
    return res;
}

double **copy_matrix(double **sm, int n){
    int i,j;
    double **tm = (double**) malloc(sizeof(double*)*n);
    if (tm==NULL) { error_message(); }
    for (i=0; i<n; i++){
        tm[i] = (double*) malloc(sizeof(double)*n);
        if (tm[i]==NULL) { error_message(); }
        for (j=0; j<n; j++){
            tm[i][j] = sm[i][j];
        }
    }
    return tm;
}
double of_f(double **m, int n){
    int i,j;
    double res = 0;
    for (i=0; i<n; i++){
        for (j=i+1; j<n; j++){
            if (j!=i){
                res+= 2*pow(m[i][j],2);
            }
        }
    }
    return res;
}

int calc_vector_dim(char *line1){
    int res = 1;
    char comma = ',';
    int i = 0;
    while (line1[i] != '\n'){
        if (line1[i] == comma){
            res++;
        }
        i++;
    }
    return res;
}

double calc_distance(double *v1 , double *v2, int dim){
    int i;
    double res = 0; double power = 2;
    for (i=0; i<dim; i++){
        res+=pow((v1[i]-v2[i]),power); }
    res = sqrt(res);
    return res;
}

double* convert_vector_string_to_double(char *line){
    int j = 0;int c=0; int offset=0;
    double *vector;
    char* num_str;
    double coordinate;

    vector= (double*) malloc(dim* sizeof(double));
    if (vector == NULL){ error_message();}

    while (line[j] != '\n'){
        while (line[j+offset] != ',' && line[j+offset] != '\n'){
            offset++;
        }
        num_str = (char*) malloc(offset * sizeof(char));
        if (num_str == NULL){ error_message();}

        offset = 0;
        while (line[j] != ',' && line[j] != '\n'){
            num_str[offset] = line[j];
            j++;
            offset++;
        }
        offset=0;
        if (line[j] != '\n'){
            j++;
        }
        coordinate = strtod(num_str, NULL);
        free(num_str);
        vector[c++] = coordinate;
    }
    return vector;
}

double** init_vectors_array(FILE *f){
    int read; char *line;
    size_t len = 0;
    double *vector;
    double **v_array;
    read = getline(&line, &len, f);
    if (read <=0 ){ error_message(); }
    dim = calc_vector_dim(line);
    n=1;
    v_array = (double**) malloc(sizeof(double*)*n);
    if (v_array==NULL) { error_message(); }
    v_array[0] = convert_vector_string_to_double(line);
    while ( (read=getline(&line, &len, f)) > 0 ) {
        vector = convert_vector_string_to_double(line);
        v_array = (double**) realloc(v_array, (sizeof(double*))*(n+1));
        if (v_array==NULL){ error_message(); }
        v_array[n] = vector;
        /* TODO- find a better condition to not take the last line */
        if (strcmp("\r\n",line)!=0){
            n++;
        }
    }
    free(line);
    return v_array;
}

void error_message(){
    printf("An Error Has Occurred\n");
    exit(1);
}


int main(int argc, char *argv[]) {
    FILE *f; double *ev; double **vectors, **am, **lm;
    char *command = argv[1];
    if (argc!=3){ error_message(); }
    f = fopen(argv[2],"r");
    if (f==NULL){ error_message(); }
    vectors = init_vectors_array(f);
    if (strcmp("jacobi",command)==0){
        ev = jacobi_algo(vectors, n);
        print_ev(ev);
        print_centroids_array(vectors, n, n);
        free(ev);
    }
    else{
        am = build_adjacency_matrix(vectors,n,dim);
        if (strcmp("wam",command)==0) {
            print_centroids_array(am, n, n);
        }
        else{
            lm = build_dd_matrix(n, am);
            if (strcmp("ddg",command)==0){
                print_centroids_array(lm, n, n);
                
            }else{
                if (strcmp("gl", command)==0){
                    build_lm_matrix(n,am,lm);
                    print_centroids_array(lm, n, n);
                }
                else{ free_2d_double_array(lm,n);free_2d_double_array(am,n);free_2d_double_array(vectors,n);error_message(); }
            }
            free_2d_double_array(lm,n);

        }
        free_2d_double_array(am,n);
    }
    free_2d_double_array(vectors,n);

    return 0;
}

void free_2d_double_array(double **arr, int len){
    int i;
    for (i=0; i < len; i++){
        free(arr[i]);
    }
    free(arr);
}

void print_centroids_array(double **cntrd, int k, int dim){
    int i,j;
    double vector_coordinate;
    double *vector_array;
    for (i=0; i<k ; i++){
        for (j=0; j<dim; j++){
            vector_array = cntrd[i];
            vector_coordinate = vector_array[j];
            if (j==(dim-1)){
                printf("%0.4f\n", vector_coordinate);
            }
            else {
                printf("%0.4f,", vector_coordinate);
            }
        }
    }
}

void print_ev(double *ev_arr){
    int i;
    for (i=0;i<n;i++){
        printf("%0.4f",ev_arr[i]);
        if (i<n-1){
            printf(",");
        }
    }
    printf("\n");
}

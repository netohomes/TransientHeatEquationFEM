//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
/** FUNCIONES Y LIBRERIAS **/
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#ifndef Heat_equation_2D_funciones_h
#define Heat_equation_2D_funciones_h
#endif
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
// Librerías
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
// Estructuras
struct elem{
    int num;
    int *conec;
    double *prop;
};

struct tipo_2{
    int num;
    double *val;
    int *ind;
    double f;
    int cols;
};
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
// Prototipo de funciones

void multiplica ( double **C, double **A, double **B, int rengA, int col, int colB );
void mult_mat_vector ( double *C, double **A, double *B, int rengA, int colA, int rengB );
void transponer_mat ( double **A, double **B, int reng, int col );
void mult_mat_escalar ( double **A, double b, int reng, int col);
void suma_mat ( double **A, double **B, int reng, int col);
double **mem_matriz ( int n, int  m, char *nombre );
double *mem_vector ( int n, char  *nombre );
void inicializa_matriz ( double **A, int reng, int col );
void inicializa_vector ( double *A, int reng );
void inicializa_mat_identidad ( double **A, int n );
void lectura ( char *r );
int *mem_vector_int ( int n, char  *nombre );
struct elem *mem_vector_elem ( int n, char  *nombre );
void inicializa_vector_int ( int *A, int reng );
int **mem_matriz_int ( int n, int  m, char *nombre );
void ady ( struct elem *elems, int n_elem, int nodos, int nodos_elem );
struct tipo_2 *mem_vector_tipo_2 ( int n, char  *nombre );
void redimensiona_vector_tipo_2 ( struct tipo_2 *A, int m,int n, char *nombre );
void redimensiona_vector ( int **A, int m, int n, char *nombre );
//void redim_vector_elem ( struct elem *a, int n, char  *nombre );
//void redim_matriz ( double **A, int n, int  m, char *nombre );
void fk_elem ( double **K_elem, double *f_elem, double **coord, double *prop, int tipo_elem,
               int dimension, int nodos_elem, int pg_elem, double **pg, double **wpg,
               double **jacobiano, double *vect_aux, double *vect_aux2, double **B,
               double **BT, double **Baux, double **K_elemaux, double **derin,
               double **difusion, double *det, int m );
void puntos_gauss ( double **pg, double **wpg, int tipo_elem );
void fjacobiano ( double **jacobiano, double **coord, double p, double n, double c, int nodos_elem,
                  int tipo_elem, int dimension, double **derin );
void fderin ( double **derin, double p, double n, double c, int tipo_elem );
double determinante ( double **A, int dimension );
void inversa ( double **A, double **inversa_A, int n, double det );
void rigid ( struct tipo_2 *tabla_ady, int tipo_elem, int nodos_elem, int pg_elem, int dimension,
             double *det );
void cond_cont( struct tipo_2 *a, int TP, int TL, int FL, int TS, int FS, int FP, int *t_points_ind,
                double *t_points_val, int *f_points_ind, double *f_points_val, int *t_line_ind,
                double *t_line_val, int *f_line_ind, double **f_line_val, int **t_surf_ind,
                double *t_surf_val, int **f_surf_ind, double **f_surf_val, int nodos, double **coord );
void grad_conj ( struct tipo_2 *A, double *x, int n, int n_iter, double tol );
void libera_matriz ( double **a, int col );
void libera_matriz_int ( int **a, int col );
void libera_vector ( double *v );
void libera_vector_int ( int *v );
void libera_elem ( struct elem *a, int n_elem );
void libera_tipo_2 ( struct tipo_2 *a, int nodos );
void libera_mem1 ( double *T, double **q_pg, double **q_pg2, double **q_nodos, double **error_PG,
                   int dimension, double *det );
void libera_mem2 ( void );
void suavizado_flujos ( struct elem *a, double *T, int dimension, int tipo_elem, int n_elem,
                        int pg_elem, int nodos_elem, double **coord, int nodos, double **q_pg,
                        double **q_nodos, double **q_pg2, double **error_PG );
void flujos_1D ( struct elem *elems, double *T, double **q_nodos, int nodos, int n_elem );
void escribe_result ( char *r, int pg_elem, int tipo_elem,  int nodos, double *T, double **q_pg,
                      double **q_pg2, double **q_nodos, double **error_PG );
void deriv_Ni( double *vect_aux, double **pg, int tipo_elem, int dimension, int j, int k );
double N ( int tipo_elem, int i, int dimension, double p, double n, double c );
void trascendencia ( struct elem *elems, struct tipo_2 *tabla_ady, double *T, double alfa, double rho,
                     double delta_t, int pasos, int n_elem, int nodos, int tipo_elem, int nodos_elem,
                     char *r, double *det, int pg_elem, int TI, int **ti_surf_ind, double *ti_surf_val,
                     int dimension );
void escribe_T ( char *r, int pg_elem, int dimension, int tipo_elem, int nodos, double *T,
                 double t, int imprime );
void escribe_flujos ( char *r, double **q_pg, double **q_nodos, double **q_pg2, double **error_PG,
                      int pg_elem, int dimension, int nodos,  double t );
void incorpora_temp ( double *T, int TP, int *t_points_ind, double *t_points_val, int TL,
                      int *t_line_ind, double *t_line_val, int TS, int **t_surf_ind,
                      double *t_surf_val, int TI, int **ti_surf_ind, double *ti_surf_val,
                      int tipo_elem, int dimension );

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
// Variables globales
struct elem *elems;
struct tipo_2 *tabla_ady;
double **coord;
int n_elem;
int nodos;

int tipo_elem;
int dimension;
int nodos_elem;
int pg_elem;

int iter;
double tol;

// Para las condiciones de contorno
int nQ, TP, FP, TL, FL, TS, FS, TI;
int *t_points_ind, *f_points_ind;
double *t_points_val, *f_points_val;
int *t_line_ind, *f_line_ind;
double *t_line_val, **f_line_val;
int **t_surf_ind, **f_surf_ind, **ti_surf_ind;
double *t_surf_val, **f_surf_val, *ti_surf_val;

// Para el caso transitorio
int pasos, transient;
double alfa, rho, delta_t;

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

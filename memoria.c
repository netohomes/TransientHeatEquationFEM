/** Asigna memoria **/
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double **mem_matriz ( int n, int  m, char *nombre )
{
	/** Descripción **/
	// Pide memoria para una matriz de n x m con valores tipo doble

	/** Variables de entrada **/
	// n		Número de renglones
	// m		Número de columnas
	// nombre	Nombre de la matriz

	/** Variables de salida **/
	// k		Apuntador que contiene la dirección de la memoria reservada

	int i;
	double **k;


	k = ( double ** ) calloc ( n, sizeof ( double * ) );
	if ( k == NULL ){
		printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
		exit ( 1 );
	}
	for ( i = 0; i < n; i++ ){
		k[i] = ( double * ) calloc ( m, sizeof ( double ) );
		if ( k == NULL ){
			printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
			exit ( 1 );
		}
	}

	return k;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
//void redim_matriz ( double **A, int n, int  m, char *nombre )
//{
	/** Descripción **/
	// Resimensiona la memoria para una matriz de n x m con valores tipo doble

	/** Variables de entrada **/
	// n		Nuevo número de renglones
	// m		Nuevo número de columnas
	// nombre	Nombre de la matriz

	/** Variables de salida **/
	// k		Apuntador que contiene la dirección de la memoria reservada
/*
	int i;
	puts ( " 5 ");
	A = ( double ** ) realloc ( A, n * sizeof ( double * ) );
	if ( A == NULL ){
		printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
		exit ( 1 );
	}
	puts ( " 6 ");
	for ( i = 0; i < n; i++ ){
		A[i] = ( double * ) realloc ( A[i], m * sizeof ( double ) );
		if ( A[i] == NULL ){
			printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
			exit ( 1 );
		}
	}
}*/
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
int **mem_matriz_int ( int n, int  m, char *nombre )
{
	/** Descripción **/
	// Pide memoria para una matriz de n x m con valores tipo entero

	/** Variables de entrada **/
	// n		Número de renglones
	// m		Número de columnas
	// nombre	Nombre de la matriz

	/** Variables de salida **/
	// k		Apuntador que contiene la dirección de la memoria reservada

	int i;
	int **k;
	k = ( int ** ) calloc ( n, sizeof ( int * ) );
	if ( k == NULL ){
		printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
		exit ( 1 );
	}
	for ( i = 0; i < n; i++ ){
		k[i] = ( int * ) calloc ( m, sizeof ( int ) );
		if ( k == NULL ){
			printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
			exit ( 1 );
		}
	}
	return k;
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double *mem_vector ( int n, char  *nombre )
{
	/** Descripción **/
	// Pide memoria para un vector de tamaño n con valores tipo doble

	/** Variables de entrada **/
	// n			Tamaño del vector
	// nombre		Nombre del vector

	/** variables de salida **/
	// v			Apuntador que contiene la dirección de la memoria reservada

	double *v;
	v = ( double *) calloc ( n, sizeof ( double ) );
	if ( v == NULL ){
		printf ( "\nMemoria insuficiente para almacenar el vector %s", nombre );
		exit ( 1 );
	}
	return v;
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
int *mem_vector_int ( int n, char  *nombre )
{
	/** Descripción **/
	// Pide memoria para un vector de tamaño n con valores tipo entero

	/** Variables de entrada **/
	// n			Tamaño del vector
	// nombre		Nombre del vector

	/** variables de salida **/
	// v			Apuntador que contiene la dirección de la memoria reservada

	int *v;
	v = ( int *) calloc ( n, sizeof ( int ) );
	if ( v == NULL ){
		printf ( "\nMemoria insuficiente para almacenar el vector %s", nombre );
		exit ( 1 );
	}
	return v;
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
struct elem *mem_vector_elem ( int n, char  *nombre )
{
	/** Descripción **/
	// Pide memoria para un vector de tamaño n con valores tipo struct elem

	/** Variables de entrada **/
	// n			Tamaño del vector
	// nombre		Nombre del vector

	/** variables de salida **/
	// v			Apuntador que contiene la dirección de la memoria reservada

	struct elem *v;
	v = ( struct elem *) calloc ( n, sizeof ( struct elem ) );
	if ( v == NULL ){
		printf ( "\nMemoria insuficiente para almacenar el vector %s", nombre );
		exit ( 1 );
	}
	return v;
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void redim_vector_elem ( struct elem *a, int n, char  *nombre )
{
	/** Descripción **/
	// Redimensiona la memoria para un vector de tamaño n con valores tipo struct elem

	/** Variables de entrada **/
	// n			Nuevo tamaño del vector
	// nombre		Nombre del vector

	/** variables de salida **/
	// a			Apuntador que contiene la dirección de la memoria reservada

	a = ( struct elem *) realloc ( a, n * sizeof ( struct elem ) );
	if ( a == NULL ){
		printf ( "\nMemoria insuficiente para almacenar el vector %s", nombre );
		exit ( 1 );
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
struct tipo_2 *mem_vector_tipo_2 ( int n, char  *nombre )
{
	/** Descripción **/
	// Pide memoria para un vector de tamaño n con valores tipo struct tipo_2

	/** Variables de entrada **/
	// n			Tamaño del vector
	// nombre		Nombre del vector

	/** variables de salida **/
	// v			Apuntador que contiene la dirección de la memoria reservada

	struct tipo_2 *v;
	v = ( struct tipo_2 *) calloc ( n, sizeof ( struct tipo_2 ) );
	if ( v == NULL ){
		printf ( "\nMemoria insuficiente para almacenar el vector %s", nombre );
		exit ( 1 );
	}
	return v;
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void redimensiona_vector_tipo_2 ( struct tipo_2 *A, int m,int n, char *nombre )
{
    /** Descripción **/
    // Aumenta o disminuye el tamaño de memoria de un vector

    /** Variables de entrada **/
    // A        Arreglo de estructuras al que pertenece el vector m
    // n        Nuevo tamaño
    // m        Vector de la estructura que será redimensionado

    A[m].ind = (int *) realloc ( A[m].ind, n*sizeof ( int ) );

    if ( A[m].ind == NULL ){
		printf ( "\nMemoria insuficiente para redimensionar el vector %s", nombre );
		exit ( 1 );
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void redimensiona_vector ( int **A, int m, int n, char *nombre )
{
    /** Descripción **/
    // Aumenta o disminuye el tamaño de memoria de un vector

    /** Variables de entrada **/
    // A        Matriz a la que pertenece el vector a redimensionar
    // n        Nuevo tamaño
    // m        Vector de la matriz que será redimensionado

    A[m] = (int *) realloc ( A[m], n*sizeof ( int ) );
    puts ( "jeje");
    if ( A[m] == NULL ){
		printf ( "\nMemoria insuficiente para redimensionar el vector %s", nombre );
		exit ( 1 );
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inicializa_matriz ( double **A, int reng, int col )
{
	/** Descripción **/
	// Inicializa una matriz [A] de tamaño reng x col a 0.0

	int i, j;

	for ( i = 0; i < reng; i++){
		for ( j = 0; j < col; j++ ){
			A[i][j] = 0.0;
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inicializa_vector ( double *A, int reng )
{
	/** Descripción **/
	// Inicializa un vector {A} de tamaño reng a 0.0 tipo doble

	int i;

	for ( i = 0; i < reng; i++){
		A[i] = 0.0;
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inicializa_vector_int ( int *A, int reng )
{
	/** Descripción **/
	// Inicializa un vector {A} de tamaño reng a 0 tipo entero

	int i;

	for ( i = 0; i < reng; i++){
		A[i] = 0;
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inicializa_mat_identidad ( double **A, int n )
{
	/** Descripción **/
	// Guarda una matriz indentidad de tamaño n x n en [A]

	int i, j;

	for ( i = 0; i < n; i++){
		for ( j = 0; j < n; j++ ){
			if ( j == i ){
				A[i][j] = 1.0;
			}
			else{
				A[i][j] = 0.0;
			}
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_matriz ( double **a, int reng )
{

    /** Descripción **/
	// Libera memoria para una matriz de n x m con valores tipo doble

	/** Variables de entrada **/
	// reng		Número de renglones de la columna ( n )
	// f		Dirección al primer elemento de la matriz

	/** Variables locales **/
	// a[][]	Apuntador a la matriz

	int i;

	for ( i = 0; i < reng; i++){
		free ( a[i] );
	}

	free ( a );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_matriz_int ( int **a, int reng )
{

    /** Descripción **/
	// Libera memoria para una matriz de n x m con valores tipo entero

	/** Variables de entrada **/
	// reng		Número de renglones de la columna ( n )
	// f		Dirección al primer elemento de la matriz

	/** Variables locales **/
	// a[][]	Apuntador a la matriz

	int i;

	for ( i = 0; i < reng; i++){
		free ( a[i] );
	}

	free ( a );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_vector ( double *v )
{
    /** Descripción **/
	// Libera memoria de un vector con valores tipo doble

	/** Variables de entrada **/
	// col		Número de columnas
	// v		Dirección al primer elemento de la matriz

	/** Variables locales **/
	// v[]	Apuntador al vector

    free ( v );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_vector_int ( int *v )
{
    /** Descripción **/
	// Libera memoria de un vector con valores tipo entero

	/** Variables de entrada **/
	// col		Número de columnas
	// v		Dirección al primer elemento de la matriz

	/** Variables locales **/
	// v[]	Apuntador al vector

    free ( v );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_elem ( struct elem *a, int n_elem )
{
    /** Descripción **/
    // Libera memoria de un arreglo de elementos tipo elem

    int m;

    for ( m = 0; m < n_elem; m++ ){
        free ( a[ m ].conec );
        free ( a[ m ].prop );
    }

    free ( a );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_tipo_2 ( struct tipo_2 *a, int nodos )
{
    /** Descripción **/
    // Libera memoria de un arreglo de elementos tipo elem

    int i;

    for ( i = 0; i < nodos; i++ ){
        free ( a[ i ].val );
        free ( a[ i ].ind );
    }

    free ( a );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_mem1 ( double *T, double **q_pg, double **q_pg2, double **q_nodos, double **error_PG,
                   int dimension, double *det )
{
    /** Descripción **/
    // Libera memoria de los arreglos de main

    // Vector de temperaturas
     libera_vector( T );

     // Vectores de flujos en P.G.
    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        libera_matriz ( q_pg, pg_elem * n_elem );
        libera_matriz ( q_pg2, pg_elem * n_elem );
        libera_matriz ( error_PG, pg_elem * n_elem );
        libera_vector ( det );
    }

    // Vector de flujos en nodos
    libera_matriz ( q_nodos, nodos );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void libera_mem2 ( void )
{
    /** Descripción **/
    // Libera memoria de los arreglos globales

    // Matriz de coordenadas
    libera_matriz ( coord, nodos );

    // Estructura de elementos
    libera_elem ( elems, n_elem );

    // Estructura tipo_2
    libera_tipo_2 ( tabla_ady, nodos );

    // De las condiciones de contorno

    // Aplicable a 1D, 2D y 3D
    if ( TP > 0 ){
        libera_vector_int ( t_points_ind );
        libera_vector ( t_points_val );

    }

    // Aplicable a 1D
    if ( FP > 0 ){
        libera_vector_int ( f_points_ind );
        libera_vector ( f_points_val );
    }

    // Aplicable a 1D, 2D y 3D
    if ( TL > 0 ){
        libera_vector_int ( t_line_ind );
        libera_vector ( t_line_val );
    }

    // Aplicable a 2D
    if ( FL > 0 ){
        libera_vector_int ( f_line_ind );
        libera_matriz ( f_line_val, FL );
    }

    // Aplicable a 3D
    if ( TS > 0 ){
        libera_matriz_int ( t_surf_ind, TS );
         libera_vector ( t_surf_val );
    }

    // Aplicable a 3D
    if ( FS > 0 ){
        libera_matriz_int ( f_surf_ind, FS );
        libera_matriz ( f_surf_val, FS );
    }

    // Aplicable a 3D
    if ( TI > 0 ){
        libera_matriz_int ( ti_surf_ind, TI );
        libera_vector ( ti_surf_val );
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


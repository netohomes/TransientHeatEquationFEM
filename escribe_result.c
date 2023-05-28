//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void escribe_result( char *r, int pg_elem, int tipo_elem,  int nodos, double *T, double **q_pg,
                     double **q_pg2, double **q_nodos, double **error_PG )
{
    /** Descripción **/
    // Escribe el archivo de resultados

    /** Variables de entrada **/
    // r                Dirección al archivo con el que se está trabajando, sin extensión.


    /** Variables locales **/
    // ruta             Dirección al archivo con el que se está trabajando, con extensión.
    // nombre_elem      Tipo de elemento

    int i, j, aux, num_menu;
    char *ruta, *nombre_elem;
    FILE *ft;

    ruta = ( char *) calloc ( 256, sizeof ( char ) );

    // Tipo de elemento
    if ( dimension == 2 ){
        if ( tipo_elem == 1 ){
            nombre_elem = "Triangle";
        }
        else{
            nombre_elem = "Quadrilateral";
        }
    }
    else{
        if ( tipo_elem == 3 ){
            nombre_elem = "Tetrahedra";
        }
        else{
            nombre_elem = "Hexahedra";
        }
    }

    strcpy ( ruta, r );
    strcat ( ruta, ".post.res" );
    //printf ( "\n\n%s\n\n", ruta );

    ft = fopen ( ruta , "w" );
    fprintf( ft, "Gid Post Results File 1.0\n\n" );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){

        fprintf( ft, "GaussPoints BoardGaussInternal ElemType %s\n", nombre_elem );
        fprintf( ft, "Number Of Gauss Points: %d\n", pg_elem );
        fprintf( ft, "Natural Coordinates: internal\n" );
        fprintf( ft, "End Gausspoints\n\n" );

    }
    num_menu = 0;
    fprintf( ft, "Result Temperature Menu %d Scalar OnNodes\n", num_menu );
    fprintf( ft, "Values\n");
    for ( i = 0; i < nodos; i++){
        fprintf( ft, "%d\t", i + 1 );
        fprintf( ft, "%lf\n", T[ i ] );
    }
    fprintf( ft, "End Values\n\n" );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        num_menu++;
        fprintf( ft, "Result FluxOnGP_1 Menu %d Vector OnGaussPoints BoardGaussInternal\n", num_menu );
        fprintf( ft, "Values\n");
        for ( i = 0; i < n_elem; i++ ){
            aux = i * pg_elem;
            for ( j = 0; j < pg_elem; j++ ){
                if ( dimension == 2 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\n", q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ] );
                    }
                }
                if ( dimension == 3 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ], q_pg[ aux + j ][ 2 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\t%lf\n", q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ], q_pg[ aux + j ][ 2 ] );
                    }
                }
            }
        }
        fprintf( ft, "End Values\n\n" );
    }
    num_menu++;
    fprintf( ft, "Result FluxOnNodes Menu %d Vector OnNodes\n", num_menu );
    if ( dimension == 1 ){
        fprintf( ft, "Componentnames qx, qy\n" );
        fprintf( ft, "Values\n");
         for ( i = 0; i < nodos; i++ ){
            fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_nodos[ i ][ 0 ], 0.000 );
        }
    }
    if ( dimension == 2 ){
        fprintf( ft, "Componentnames qx, qy\n" );
        fprintf( ft, "Values\n");
        for ( i = 0; i < nodos; i++ ){
            fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_nodos[ i ][ 0 ], q_nodos[ i ][ 1 ] );
        }
    }
    if ( dimension == 3 ){
        fprintf( ft, "Componentnames qx, qy, qz\n" );
        fprintf( ft, "Values\n");
        for ( i = 0; i < nodos; i++ ){
            fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, q_nodos[ i ][ 0 ], q_nodos[ i ][ 1 ], q_nodos[ i ][ 2 ] );
        }
    }
    fprintf( ft, "End Values\n\n" );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        num_menu++;
        fprintf( ft, "Result FluxOnGP_2 Menu %d Vector OnGaussPoints BoardGaussInternal\n", num_menu );
        fprintf( ft, "Values\n");
        for ( i = 0; i < n_elem; i++ ){
            aux = i * pg_elem;
            for ( j = 0; j < pg_elem; j++ ){
                if ( dimension == 2 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\n", q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ] );
                    }
                }
                if ( dimension == 3 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ], q_pg2[ aux + j ][ 2 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\t%lf\n", q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ], q_pg2[ aux + j ][ 2 ] );
                    }
                }
            }
        }
        fprintf( ft, "End Values\n\n" );
    }

     if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        num_menu++;
        fprintf( ft, "Result ErrorOnGP Menu %d Vector OnGaussPoints BoardGaussInternal\n", num_menu );
        fprintf( ft, "Values\n");
        for ( i = 0; i < n_elem; i++ ){
            aux = i * pg_elem;
            for ( j = 0; j < pg_elem; j++ ){
                if ( dimension == 2 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\n", i + 1, error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\n", error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ] );
                    }
                }
                if ( dimension == 3 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ], error_PG[ aux + j ][ 2 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\t%lf\n", error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ], error_PG[ aux + j ][ 2 ] );
                    }
                }
            }
        }
        fprintf( ft, "End Values\n\n" );
    }

    fclose ( ft );

    // Libera nombre de archivo
    free ( ruta );


}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void escribe_T ( char *r, int pg_elem, int dimension, int tipo_elem, int nodos, double *T,
                 double t, int imprime )
{
    /** Descripción **/
    // Escribe en el archivo de resultados los valores de T en el tiempo

    /** Variables de entrada **/
    // r                Dirección al archivo con el que se está trabajando, sin extensión.


    /** Variables locales **/
    // ruta             Dirección al archivo con el que se está trabajando, con extensión.
    // nombre_elem      Tipo de elemento

    int i;
    char *ruta, *nombre_elem;
    FILE *ft;

    ruta = ( char *) calloc ( 256, sizeof ( char ) );

    // Tipo de elemento
    if ( dimension == 2 ){
        if ( tipo_elem == 1 ){
            nombre_elem = "Triangle";
        }
        else{
            nombre_elem = "Quadrilateral";
        }
    }
    else{
        if ( tipo_elem == 3 ){
            nombre_elem = "Tetrahedra";
        }
        else{
            nombre_elem = "Hexahedra";
        }
    }

    strcpy ( ruta, r );
    strcat ( ruta, ".post.res" );
    //printf ( "\n\n%s\n\n", ruta );

    ft = fopen ( ruta , "a+" );

    if ( imprime == 1 ){
        fprintf( ft, "Gid Post Results File 1.0\n\n" );
        if ( ( dimension == 2 ) || ( dimension == 3 ) ){

            fprintf( ft, "GaussPoints BoardGaussInternal ElemType %s\n", nombre_elem );
            fprintf( ft, "Number Of Gauss Points: %d\n", pg_elem );
            fprintf( ft, "Natural Coordinates: internal\n" );
            fprintf( ft, "End Gausspoints\n\n" );

        }
    }

    fprintf( ft, "Result Temperature Menu %lf Scalar OnNodes\n", t );
    fprintf( ft, "Values\n");
    for ( i = 0; i < nodos; i++){
        fprintf( ft, "%d\t", i + 1 );
        fprintf( ft, "%lf\n", T[ i ] );
    }
    fprintf( ft, "End Values\n\n" );

    fclose ( ft );

    // Libera nombre de archivo
    free ( ruta );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void escribe_flujos ( char *r, double **q_pg, double **q_nodos, double **q_pg2, double **error_PG,
                      int pg_elem, int dimension, int nodos,  double t )
{
    /** Descripción **/
    // Escribe en el archivo de resultados los flujos y el error para el caso transitorio

    int i, j, aux;
    char *ruta;
    FILE *ft;

    ruta = ( char *) calloc ( 256, sizeof ( char ) );

    strcpy ( ruta, r );
    strcat ( ruta, ".post.res" );

    ft = fopen ( ruta , "a+" );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        fprintf( ft, "Result FluxOnGP_1 Menu %lf Vector OnGaussPoints BoardGaussInternal\n", t );
        fprintf( ft, "Values\n");
        for ( i = 0; i < n_elem; i++ ){
            aux = i * pg_elem;
            for ( j = 0; j < pg_elem; j++ ){
                if ( dimension == 2 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\n", q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ] );
                    }
                }
                if ( dimension == 3 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ], q_pg[ aux + j ][ 2 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\t%lf\n", q_pg[ aux + j ][ 0 ], q_pg[ aux + j ][ 1 ], q_pg[ aux + j ][ 2 ] );
                    }
                }
            }
        }
        fprintf( ft, "End Values\n\n" );
    }

    fprintf( ft, "Result FluxOnNodes Menu %lf Vector OnNodes\n", t );
    if ( dimension == 1 ){
        fprintf( ft, "Componentnames qx, qy\n" );
        fprintf( ft, "Values\n");
         for ( i = 0; i < nodos; i++ ){
            fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_nodos[ i ][ 0 ], 0.000 );
        }
    }
    if ( dimension == 2 ){
        fprintf( ft, "Componentnames qx, qy\n" );
        fprintf( ft, "Values\n");
        for ( i = 0; i < nodos; i++ ){
            fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_nodos[ i ][ 0 ], q_nodos[ i ][ 1 ] );
        }
    }
    if ( dimension == 3 ){
        fprintf( ft, "Componentnames qx, qy, qz\n" );
        fprintf( ft, "Values\n");
        for ( i = 0; i < nodos; i++ ){
            fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, q_nodos[ i ][ 0 ], q_nodos[ i ][ 1 ], q_nodos[ i ][ 2 ] );
        }
    }
    fprintf( ft, "End Values\n\n" );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        fprintf( ft, "Result FluxOnGP_2 Menu %lf Vector OnGaussPoints BoardGaussInternal\n", t );
        fprintf( ft, "Values\n");
        for ( i = 0; i < n_elem; i++ ){
            aux = i * pg_elem;
            for ( j = 0; j < pg_elem; j++ ){
                if ( dimension == 2 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\n", i + 1, q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\n", q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ] );
                    }
                }
                if ( dimension == 3 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ], q_pg2[ aux + j ][ 2 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\t%lf\n", q_pg2[ aux + j ][ 0 ], q_pg2[ aux + j ][ 1 ], q_pg2[ aux + j ][ 2 ] );
                    }
                }
            }
        }
        fprintf( ft, "End Values\n\n" );
    }

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        fprintf( ft, "Result ErrorOnGP Menu %lf Vector OnGaussPoints BoardGaussInternal\n", t );
        fprintf( ft, "Values\n");
        for ( i = 0; i < n_elem; i++ ){
            aux = i * pg_elem;
            for ( j = 0; j < pg_elem; j++ ){
                if ( dimension == 2 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\n", i + 1, error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\n", error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ] );
                    }
                }
                if ( dimension == 3 ){
                    if ( j == 0 ){
                        fprintf( ft, "%d\t%lf\t%lf\t%lf\n", i + 1, error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ], error_PG[ aux + j ][ 2 ] );
                    }
                    else{
                        fprintf( ft, "\t%lf\t%lf\t%lf\n", error_PG[ aux + j ][ 0 ], error_PG[ aux + j ][ 1 ], error_PG[ aux + j ][ 2 ] );
                    }
                }
            }
        }
        fprintf( ft, "End Values\n\n" );
    }

    fclose ( ft );

    // Libera nombre de archivo
    free ( ruta );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

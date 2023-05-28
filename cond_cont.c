//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void dirichlet ( struct tipo_2 *a, int TP, int TL, int TS, int *t_points_ind, double *t_points_val,
                 int *t_line_ind, double *t_line_val, int **t_surf_ind, double *t_surf_val, int nodos );
void neumann ( struct tipo_2 *a, double **coord, int FL, int FS, int FP, int *f_line_ind,
               double **f_line_val, int **f_surf_ind, double **f_surf_val, int nodos,
               int *f_points_ind, double *f_points_val );
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void cond_cont( struct tipo_2 *a, int TP, int TL, int FL, int TS, int FS, int FP, int *t_points_ind,
                double *t_points_val, int *f_points_ind, double *f_points_val, int *t_line_ind,
                double *t_line_val, int *f_line_ind, double **f_line_val, int **t_surf_ind,
                double *t_surf_val, int **f_surf_ind, double **f_surf_val, int nodos, double **coord )
{
    /** Descripción **/
    // Impone las condiciones de contorno: Dirichlet y Neumann

    // Compara que no haya un flujo y una temperatura definidas sobre el mismo nodo
    //for ( i = 0; i < P)

    neumann ( a, coord, FL, FS, FP, f_line_ind, f_line_val, f_surf_ind, f_surf_val, nodos,
              f_points_ind, f_points_val );

    dirichlet ( a, TP, TL, TS, t_points_ind, t_points_val, t_line_ind, t_line_val, t_surf_ind,
                t_surf_val, nodos );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void dirichlet ( struct tipo_2 *a, int TP, int TL, int TS, int *t_points_ind, double *t_points_val,
                 int *t_line_ind, double *t_line_val, int **t_surf_ind, double *t_surf_val, int nodos )
{
     /** Descripción **/
    // Impone las condiciones de contorno: Dirichlet definidas en líneas o puntos sobre la
    // la geometría o la malla

    int i, j, k, c, aux1, aux3, continua ;
    double aux2;
    double *T_val;
    int *T_ind;

    // UNIFICA TODAS LAS TEMPERATURAS

    // Pide memoria para el máximo posible
    if ( ( dimension == 1 ) || ( dimension == 2 ) ){
        T_val = mem_vector ( TP + TL, "auxiliar para guardar temp. impuestas");
        T_ind = mem_vector_int ( TP + TL, "auxiliar para guardar temp. impuestas");
        inicializa_vector ( T_val, TP + TL );
        inicializa_vector_int ( T_ind, TP + TL );
    }
    else{
        if ( tipo_elem == 3 ){
            T_val = mem_vector ( TP + TL + TS * 3, "auxiliar para guardar temp. impuestas");
            T_ind = mem_vector_int ( TP + TL + TS * 3, "auxiliar para guardar temp. impuestas");
            inicializa_vector ( T_val, TP + TL + TS * 3 );
            inicializa_vector_int ( T_ind, TP + TL + TS * 3 );
        }
        if ( tipo_elem == 4 ){
            T_val = mem_vector ( TP + TL + TS * 4, "auxiliar para guardar temp. impuestas");
            T_ind = mem_vector_int ( TP + TL + TS * 4, "auxiliar para guardar temp. impuestas");
            inicializa_vector ( T_val, TP + TL + TS * 4 );
            inicializa_vector_int ( T_ind, TP + TL + TS * 4 );
        }
    }

    // Incorpora las temperaturas de puntos
    for ( i = 0; i < TP; i++ ){
        T_ind[ i ] = t_points_ind[ i ];
        T_val[ i ] = t_points_val[ i ];
    }

    // Suma las temperaturas si hay de línea y de punto sobre un mismo punto,
    // sino solo la incorpora.
    aux1 = TP;

    //printf( "\nTP = %d", TP );

    //printf( "\nTL = %d", TL );

    // Incorpora las temperaturas sobre lineas
    for ( i = 0; i < TL; i++ ){
        continua = 0;
        for ( j = 0; j < TP; j++ ){
            if ( t_line_ind[ i ] == t_points_ind[ j ] ){
                T_val[ j ] += t_line_val[ i ];
                //printf( "\nT_val[ %d ] = % lf", T_ind[ t_points_ind[j] ], T_val[ t_points_ind[j] ] );
                //printf( "\nt_line_val[ %d ] = % lf", t_line_ind[ i ], t_line_val[ i ] );
                continua = 1;
                break;
            }
        }
        if ( continua == 1 ) continue;
        T_ind[ aux1 ] = t_line_ind[ i ];
        T_val[ aux1 ] = t_line_val[ i ];
        //printf( "\nt_line_ind = %d", t_line_ind[ i ] );
        //printf( "\nt_line_val = %lf", t_line_val[ i ] );
        aux1 ++; // cuenta los nodos con temp. impuesta
    }

    // Integra para 3D las temperaturas de superficie
    // Suma los valores
    aux3 = aux1;
    if ( dimension == 3 ){
        for ( i = 0; i < TS; i++ ){
            if ( tipo_elem == 3 ){
                for ( j = 0; j < 3; j++ ){
                    continua = 0;
                    for ( k = 0; k < aux1; k++ ){
                        if ( T_ind[ k ] == t_surf_ind[ i ][ j ] ){
                            T_val[ k ] += t_surf_val[ i ];
                            continua = 1;
                            break;
                        }
                    }
                    if ( continua == 1 ) continue;
                    T_ind[ aux3 ] = t_surf_ind[ i ][ j ];
                    T_val[ aux3 ] += t_surf_val[ i ];
                    aux3++;
                }
            }
            if ( tipo_elem == 4 ){
                for ( j = 0; j < 4; j++ ){
                    continua = 0;
                    for ( k = 0; k < aux1; k++ ){
                        if ( T_ind[ k ] == t_surf_ind[ i ][ j ] ){
                            T_val[ k ] += t_surf_val[ i ];
                            continua = 1;
                            break;
                        }
                    }
                    if ( continua == 1 ) continue;
                    T_ind[ aux3 ] = t_surf_ind[ i ][ j ];
                    T_val[ aux3 ] += t_surf_val[ i ];
                    aux3++;
                }
            }
        }
    }

    /*printf ( "\nCondidiones de temperatura");
    for ( i = 0; i < aux3; i++ ){
        printf( "\n%d\t%lf", T_ind[ i ], T_val[ i ] );

    }*/

    // Ordena el vector de temperaturas
    //for ( i = 0; i < aux1, i++){
    //printf ( "\naux1 = %d", aux1 );

    // IMPONE LAS CONDICIONES SOBRE EL VECTOR DE FLUJOS Y MODIFICA LA MATRIZ DE RIGIDEZ
    // COMPRIMIDA DE ACUERDO A LAS CONDICONES DIRICHLET IMPUESTAS
    if ( dimension == 1 ){
        aux1 = aux1;
    }
    else if ( dimension == 2 ){
        aux1 = aux1;
    }
    else{
        aux1 = aux3;
    }

    for ( i = 0; i < aux1; i++ ){
        c = T_ind[ i ]; // columna
        for ( j = 0; j < nodos; j++ ){
            for ( k = 0; k < a[ j ].cols; k++ ){
                if ( a[ j ].ind[ k ] == c ){
                    aux2 = a[ j ].val[ k ] * T_val[ i ];
                    a[ j ].f -= aux2;
                    a[ j ].val [ k ] = 0.0;

                    //printf ( "\nT_val = %lf", T_val[ i ] );
                    //printf ( "\na[%d].val[%d] = %lf", j, k, a[ j ].val[ k ] );

                }
            }
        }
    }

    for ( i = 0; i < aux1; i++ ){
        c = T_ind[ i ]; // columna
        for ( j = 0; j < nodos; j++ ){
            if ( j == c ){
                a[ j ].f = T_val[ i ];
                for ( k = 0; k < a[ j ].cols; k++){
                    if ( k == 0 ){
                        a[ j ].val[ k ] = 1.0;
                    }
                    else{
                         a[ j ].val[ k ] = 0.0;
                    }

                }
            }
        }
    }

    /*printf ( "\n\nVECTOR DE FLUJOS CON C.C. DE TEMP." );
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%lf", tabla_ady[ i ].f );
    }*/

    /*printf ( "\n\nMATRIZ DE RIGIDEZ COMPRIMIDA CON C.C. DE TEMP. \n" );
    for ( i = 0; i < nodos; i++){
        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){
            printf ( "%lf\t", tabla_ady[ i ].val[ j ] );
            if ( j == tabla_ady[ i ].cols -1 ) printf( "\n");
        }
    }*/

     // Libera memoria
     libera_vector ( T_val );
     libera_vector_int ( T_ind );
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void neumann ( struct tipo_2 *a, double **coord, int FL, int FS, int FP, int *f_line_ind,
               double **f_line_val, int **f_surf_ind, double **f_surf_val, int nodos,
               int *f_points_ind, double *f_points_val )
{
    /** Descripción **/
    // Impone las condiones tipo Neumann en el vector de flujos

    /** Variables locales **/

    int b, c, d, e, i, j, k;
    double m1, m2, Fx, Fy, Fn;
    double aux1, aux2;
    double Ma, Mb, Mc, Mn, h1, h2, b1, b2, a_b, b_c, area;
    double l;
    double ax, ay, az, bx, by, bz, cx, cy, cz, fx, fy, fz;
    double nx, ny, nz;

    /*printf ( "\n\nFLUJOS DE LINEA" );
    for ( i = 0; i < FL; i++ ){
        printf ( "\n%d\t%lf\t%lf", f_line_ind[ i ], f_line_val[ i ][ 0 ], f_line_val[ i ][ 1 ] );
    }*/


    if ( dimension == 1 ){
        for ( i = 0; i < FP; i++ ){
            c = f_points_ind[ i ];
            a[ c ].f += f_points_val[ i ];
        }
    }
    else if ( dimension == 2 ){
        // Para flujos sobre líneas, caso 2D
        for ( i = 0; i < FL; i++ ){
            c = f_line_ind[ i ]; // Nodo
            for ( j = i + 1; j < FL; j++ ){
                b = f_line_ind[ j ];
                for ( k = 0; k < a[ c ].cols; k++ ){
                    if ( a[ c ].ind[ k ] == b ){
                        aux1 = coord[ b ][ 1 ] - coord[ c ][ 1 ];
                        aux2 = coord[ b ][ 0 ] - coord[ c ][ 0 ];
                        l = sqrt ( aux1 * aux1 + aux2 * aux2 );
                        Fx = f_line_val[ i ][ 0 ];
                        Fy = f_line_val[ i ][ 1 ];
                        if ( aux2 < 0.0001 ){
                            Fn = Fx;
                        }
                        else if ( aux1 < 0.0001 ){
                            Fn = Fy;
                        }
                        else{
                            m1 = fabs( aux1 / aux2 );
                            m2 = 1.0 / m1;
                            Fn = ( 1.0 * Fx + m2 * Fy ) / ( sqrt ( 1.0 + m2 * m2 ) );
                        }

                        //printf ( "\nm1 = %lf", m1 );
                        //printf ( "\nm2 = %lf", m2 );
                        //printf ( "\nl = %lf", l );
                        //printf ( "\nFn = %lf", Fn );

                        a[ c ].f += 0.5 * Fn * l;
                        a[ b ].f += 0.5 * Fn * l;
                    }
                }
            }
        }
    }
    else{
        // Para flujos sobre superficies, caso 3D
        for ( i = 0; i < FS; i++ ){

            b = f_surf_ind[ i ][ 0 ];
            c = f_surf_ind[ i ][ 1 ];
            d = f_surf_ind[ i ][ 2 ];

            if ( tipo_elem == 4 ) e = f_surf_ind[ i ][ 3 ];

            // Vector de flujos sobre la cara
            fx = f_surf_val[ i ][ 0 ];
            fy = f_surf_val[ i ][ 1 ];
            fz = f_surf_val[ i ][ 2 ];

            // Coordenadas de vectores a y b
            ax = coord[ c ][ 0 ] - coord[ b ][ 0 ];
            ay = coord[ c ][ 1 ] - coord[ b ][ 1 ];
            az = coord[ c ][ 2 ] - coord[ b ][ 2 ];
            bx = coord[ d ][ 0 ] - coord[ b ][ 0 ];
            by = coord[ d ][ 1 ] - coord[ b ][ 1 ];
            bz = coord[ d ][ 2 ] - coord[ b ][ 2 ];

            // Producto punto c = a x b
            nx = ay * bz - az * by;
            ny = az * bx - ax * bz;
            nz = ax * by - ay * bx;

            // Producto punto Fn = c . f
            // f es el vector de flujos actuando sobre la cara
            Mn = sqrt ( nx * nx + ny * ny + nz * nz );
            Fn = ( nx * fx + ny * fy + nz * fz ) / Mn;

            // Area del elemento
            //a.b
            a_b = ax * bx + ay * by + az * bz;

            // Magnitudes
            Ma = sqrt ( ax * ax + ay * ay + az * az );
            Mb = sqrt ( bx * bx + by * by + bz * bz );

            // Bases
            b1 = a_b / Mb;

            // Alturas
            h1 = sqrt ( Ma * Ma - b1 * b1 );

            // Areas
            area = 0.5 * Mb * h1;

            if ( tipo_elem == 4 ){

                // Componentes vector
                cx = coord[ e ][ 0 ] - coord[ b ][ 0 ];
                cy = coord[ e ][ 1 ] - coord[ b ][ 1 ];
                cz = coord[ e ][ 2 ] - coord[ b ][ 2 ];

                // Magnitud
                Mc = sqrt ( cx * cx + cy * cy + cz * cz );

                //b.c
                b_c = bx * cx + by * cy + bz * cz;

                // Bases
                b2 = b_c / Mb;

                // Alturas
                h2 = sqrt ( Mc * Mc - b2 * b2 );

                area += 0.5 * Mb * h2;

            }

            // Integra sobre la superficie en cuestión
            if ( tipo_elem == 3 ){
                for ( j = 0; j < 3; j++ ){
                    b = f_surf_ind[ i ][ j ];
                    a[ b ].f += 1.0/3.0 * Fn * area;
                }
            }
            else{
                 for ( j = 0; j < 4; j++ ){
                    b = f_surf_ind[ i ][ j ];
                    a[ b ].f += 1.0/4.0 * Fn * area;
                }
            }
        }
    }

    /*printf ( "\n\nVECTOR DE FLUJOS CON CC NEUMANN." );

    for ( i = 0; i < nodos; i++ ){
        printf ( "\n %d\t %lf", i, a[ i ].f  );
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void incorpora_temp ( double *T, int TP, int *t_points_ind, double *t_points_val, int TL,
                      int *t_line_ind, double *t_line_val, int TS, int **t_surf_ind,
                      double *t_surf_val, int TI, int **ti_surf_ind, double *ti_surf_val,
                      int tipo_elem, int dimension )
{
     /** Descripción **/
     // Establece el vector de temperaturas inicial

    int i, j, k, a, aux1, aux3, continua, *T_ind;
    double *T_val;


     // Pide memoria para el máximo posible
    if ( ( dimension == 1 ) || ( dimension == 2 ) ){
        T_val = mem_vector ( TP + TL, "auxiliar para guardar temp. impuestas");
        T_ind = mem_vector_int ( TP + TL, "auxiliar para guardar temp. impuestas");
        inicializa_vector ( T_val, TP + TL );
        inicializa_vector_int ( T_ind, TP + TL );
    }
    else{
        if ( tipo_elem == 3 ){
            T_val = mem_vector ( TP + TL + TS * 3, "auxiliar para guardar temp. impuestas");
            T_ind = mem_vector_int ( TP + TL + TS * 3, "auxiliar para guardar temp. impuestas");
            inicializa_vector ( T_val, TP + TL + TS * 3 );
            inicializa_vector_int ( T_ind, TP + TL + TS * 3 );
        }
        if ( tipo_elem == 4 ){
            T_val = mem_vector ( TP + TL + TS * 4, "auxiliar para guardar temp. impuestas");
            T_ind = mem_vector_int ( TP + TL + TS * 4, "auxiliar para guardar temp. impuestas");
            inicializa_vector ( T_val, TP + TL + TS * 4 );
            inicializa_vector_int ( T_ind, TP + TL + TS * 4 );
        }
    }

    // Incorpora las temperaturas de puntos
    for ( i = 0; i < TP; i++ ){
        T_ind[ i ] = t_points_ind[ i ];
        T_val[ i ] = t_points_val[ i ];
    }

    // Suma las temperaturas si hay de línea y de punto sobre un mismo punto,
    // sino solo la incorpora.
    aux1 = TP;

    //printf( "\nTP = %d", TP );

    //printf( "\nTL = %d", TL );

    // Incorpora las temperaturas sobre lineas
    for ( i = 0; i < TL; i++ ){
        continua = 0;
        for ( j = 0; j < TP; j++ ){
            if ( t_line_ind[ i ] == t_points_ind[ j ] ){
                T_val[ j ] += t_line_val[ i ];
                //printf( "\nT_val[ %d ] = % lf", T_ind[ t_points_ind[j] ], T_val[ t_points_ind[j] ] );
                //printf( "\nt_line_val[ %d ] = % lf", t_line_ind[ i ], t_line_val[ i ] );
                continua = 1;
                break;
            }
        }
        if ( continua == 1 ) continue;
        T_ind[ aux1 ] = t_line_ind[ i ];
        T_val[ aux1 ] = t_line_val[ i ];
        //printf( "\nt_line_ind = %d", t_line_ind[ i ] );
        //printf( "\nt_line_val = %lf", t_line_val[ i ] );
        aux1 ++; // cuenta los nodos con temp. impuesta
    }

    // Integra para 3D las temperaturas de superficie
    // Suma los valores
    aux3 = aux1;
    if ( dimension == 3 ){
        for ( i = 0; i < TS; i++ ){
            if ( tipo_elem == 3 ){
                for ( j = 0; j < 3; j++ ){
                    continua = 0;
                    for ( k = 0; k < aux1; k++ ){
                        if ( T_ind[ k ] == t_surf_ind[ i ][ j ] ){
                            T_val[ k ] += t_surf_val[ i ];
                            continua = 1;
                            break;
                        }
                    }
                    if ( continua == 1 ) continue;
                    T_ind[ aux3 ] = t_surf_ind[ i ][ j ];
                    T_val[ aux3 ] += t_surf_val[ i ];
                    aux3++;
                }
            }
            if ( tipo_elem == 4 ){
                for ( j = 0; j < 4; j++ ){
                    continua = 0;
                    for ( k = 0; k < aux1; k++ ){
                        if ( T_ind[ k ] == t_surf_ind[ i ][ j ] ){
                            T_val[ k ] += t_surf_val[ i ];
                            continua = 1;
                            break;
                        }
                    }
                    if ( continua == 1 ) continue;
                    T_ind[ aux3 ] = t_surf_ind[ i ][ j ];
                    T_val[ aux3 ] += t_surf_val[ i ];
                    aux3++;
                }
            }
        }
    }

    for ( i = 0; i < aux3; i++ ){
        a = T_ind[ i ];
        T[ a ] = T_val[ i ];
    }

    for ( i = 0; i < TI; i++ ){
        if ( tipo_elem == 3 ){
            for ( j = 0; j < 3; j++ ){
                a = ti_surf_ind[ i ][ j ];
                T[ a ] = ti_surf_val[ j ];
            }
        }
        if ( tipo_elem == 4 ){
            for ( j = 0; j < 4; j++ ){
                a = ti_surf_ind[ i ][ j ];
                T[ a ] = ti_surf_val[ j ];
            }
        }
    }

    /*printf ( "\nCondidiones de temperatura");
    for ( i = 0; i < aux3; i++ ){
        printf( "\n%d\t%lf", T_ind[ i ], T_val[ i ] );

    }*/

    /*printf ( "\nT inicial" );
    for ( i = 0; i < nodos; i++ ){
        printf( "\n%d\t%lf", i, T[ i ] );

    }*/

    // Ordena el vector de temperaturas
    //for ( i = 0; i < aux1, i++){
    //printf ( "\naux1 = %d", aux1 );

    libera_vector_int ( T_ind );
    libera_vector ( T_val );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

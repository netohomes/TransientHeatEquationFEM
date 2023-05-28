//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void M_K ( struct tipo_2 *tabla_ady, struct elem *elems, double *m_elem, double alfa, int n_elem,
         int nodos_elem, double **C, double **Caux, double *det );
void b ( struct tipo_2 *tabla_ady, double **Caux, double *b1, double *T, int nodos );
void masas ( double *m_elem, double **pg, double **wpg, int nodos_elem, int tipo_elem, int dimension,
             double rho, int pg_elem );
void ensambla_rigid ( struct tipo_2 *tabla_ady, struct elem *elems, double *m_elem, int nodos_elem,
                      int m, double alfa, double **C, double **Caux );
double **mem_C ( struct tipo_2 *tabla_ady, int nodos, char *nombre );
void mult_MKT ( struct tipo_2 *A, double **val, double *x, double *b, int n );
void inicia_C ( struct tipo_2 *tabla_ady, double **A, int n );
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void trascendencia ( struct elem *elems, struct tipo_2 *tabla_ady, double *T, double alfa, double rho,
                     double delta_t, int pasos, int n_elem, int nodos, int tipo_elem, int nodos_elem,
                     char *r, double *det, int pg_elem, int TI, int **ti_surf_ind, double *ti_surf_val,
                     int dimension )
{
    /** Descripción **/
    // Realiza el proceso de trascendencia comprobado para 3D
    // Calcula las temperaturas y los flujos

    double **pg, **wpg, **C, **Caux, *m_elem, *F;
    double t;
    int i, iter;
    //int j;
    //int k;

    pg = mem_matriz ( pg_elem, dimension,  "de puntos de Gauss" );
    wpg = mem_matriz ( pg_elem, dimension,  "de pesos de Gauss" );
    puntos_gauss ( pg, wpg, tipo_elem );

    inicializa_vector ( T, nodos ); // Falta inicializar a lo indicado desde GID

    incorpora_temp ( T, TP, t_points_ind, t_points_val, TL, t_line_ind, t_line_val, TS, t_surf_ind,
                     t_surf_val, TI, ti_surf_ind, ti_surf_val, tipo_elem, dimension );

    t = 0.0;

    escribe_T ( r, pg_elem, dimension, tipo_elem, nodos, T, t, 1 );

    m_elem = mem_vector ( nodos_elem, "de masas." );
    masas ( m_elem, pg, wpg, nodos_elem, tipo_elem, dimension, rho, pg_elem );

    /*printf( "\nMasas\n" );
    for ( i = 0; i < nodos_elem; i++ ){
        printf( "\n%lf", m_elem[ i ] );
    }*/

    for ( i = 0; i < nodos_elem; i++ ){
        m_elem[ i ] *= 1.0/delta_t;
    }

    C = mem_C ( tabla_ady, nodos, "C" );
    Caux = mem_C ( tabla_ady, nodos, "Caux" );

    F = mem_vector ( nodos, "b1." );
    for ( i = 0; i < nodos; i++ ){
        F[ i ] = tabla_ady[ i ].f;
    }

    M_K ( tabla_ady, elems, m_elem, alfa, n_elem, nodos_elem, C, Caux, det );

    inicia_C ( tabla_ady, C, nodos );

    for ( i = 0; i < pasos; i++ ){

        b ( tabla_ady, Caux, F, T, nodos );

       /* printf ( "\n\nVECTOR FLUJOS" );
        for ( j = 0; j < nodos; j++ ){
            printf ( "\n%d\t%lf", j + 1, tabla_ady[ j].f );
        }*/

        iter = nodos;
        grad_conj ( tabla_ady, T, nodos, iter, tol );

        /*printf ( "\n\nVECTOR SOLUCION" );
        for ( j = 0; j < nodos; j++ ){
            printf ( "\n%d\t%lf", j + 1, T[ j ] );
        }*/

        t += delta_t;

        escribe_T ( r, pg_elem, dimension, tipo_elem, nodos, T, t, 0 );

    }

    libera_matriz ( pg, pg_elem );
    libera_matriz ( wpg, pg_elem );
    libera_matriz ( C, nodos );
    libera_matriz ( Caux, nodos );
    libera_vector ( m_elem );
    libera_vector ( F );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void M_K ( struct tipo_2 *tabla_ady, struct elem *elems, double *m_elem, double alfa, int n_elem,
         int nodos_elem, double **C, double **Caux, double *det )
{
    /** Descripción **/
    // Obtine las matrices C:
    // [C] = 1 / delta_t * [ M ] + alfa * [ K ]
    // [Caux] = 1 / delta_t * [ M ] - ( 1 - alfa ) * [ K ]

    int m, i, j;
    double aux1, aux2;

    for ( i = 0; i < nodos; i++){

        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){

            aux1 = alfa * tabla_ady[ i ].val[ j ];

            aux2 = - ( 1.0 - alfa ) * tabla_ady[ i ].val[ j ];

            C[ i ][ j ] = aux1;

            Caux[ i ][ j ] = aux2;
        }
    }

    for ( m = 0; m < n_elem; m++ ){

        for ( i = 0; i < nodos_elem; i++ ) m_elem[ i ] *= det[ m ];

        //mult_mat_escalar ( m_elem, det[ m ], nodos_elem, nodos_elem );

        ensambla_rigid( tabla_ady, elems, m_elem, nodos_elem, m, alfa, C, Caux );

    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void b ( struct tipo_2 *tabla_ady, double **Caux, double *b1, double *T, int nodos )
{
    int i;
    double *b2;

    b2 = mem_vector ( nodos, "b2." );

    mult_MKT ( tabla_ady, Caux, T, b2, nodos );

    /*printf( "\nFlujod b1" );
    for ( i = 0; i < nodos; i++ ){

        printf( "\n%lf", b1[ i ] );
    }*/

    for ( i = 0; i < nodos; i++ ){
        tabla_ady[ i ].f = b1[ i ] + b2[ i ];
    }

    libera_vector ( b2 );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void masas ( double *m_elem, double **pg, double **wpg, int nodos_elem, int tipo_elem, int dimension,
             double rho, int pg_elem )
{
    /** Descripción **/
    // Calcula la matriz de masas

     /** Variables locales **/
    double p = 0.0, n = 0.0, c = 0.0, wp = 0.0, wn = 0.0, wc = 0.0;
    double Nreng, Ncol, aux;
    int k, i, j;

    inicializa_vector ( m_elem, nodos_elem );

    for ( k = 0; k < pg_elem; k++ ){

        if ( dimension == 1 ){
            printf( "\nEntra a 1D, para asignar p,n,c... en trascencencia" );
            exit( 1 );
	    }
        else if ( dimension == 2 ){
            p  = pg[ k ][ 0 ];
            n  = pg[ k ][ 1 ];
            wp = wpg[ k ][ 0 ];
            wn = wpg[ k ][ 1 ];
	    }
	    else {
            p  = pg[ k ][ 0 ];
            n  = pg[ k ][ 1 ];
            c  = pg[ k ][ 2 ];
            wp = wpg[ k ][ 0 ];
            wn = wpg[ k ][ 1 ];
            wc = wpg[ k ][ 2 ];
	    }

        //inicializa_vector ( M_aux, nodos_elem );

        for ( i = 0; i < nodos_elem; i++ ){

            aux = 0.0;

            if ( dimension == 2 ){
                Nreng = N ( tipo_elem, i, dimension, p, n, 0.0 );
            }
            if ( dimension == 3 ){
                Nreng = N ( tipo_elem, i, dimension, p, n, c );
            }

            for ( j = 0; j < nodos_elem; j++ ){

                if ( dimension == 2 ){
                    Ncol = N ( tipo_elem, j, dimension, p, n, 0.0 );
                }
                if ( dimension == 3 ){
                    Ncol = N ( tipo_elem, j, dimension, p, n, c );
                }

                if (  dimension == 2 ){

                    //m_elem[ i ][ j ] += Nreng * Ncol * wp * wn;
                    //m_elem[ i ][ j ] *= rho;
                    aux += Nreng * Ncol * wp * wn * rho;
                }
                if (  dimension == 3 ){
                    //m_elem[ i ][ j ] += Nreng * Ncol * wp * wn * wc;
                    //m_elem[ i ][ j ] *= rho;
                     aux += Nreng * Ncol * wp * wn * wc * rho;

                }

            }

            m_elem[ i ] += aux;

        }

    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void ensambla_rigid ( struct tipo_2 *tabla_ady, struct elem *elems, double *m_elem, int nodos_elem,
                      int m, double alfa, double **C, double **Caux )
{
    /** Descripción **/
   // Obtiene las matricez Cder y Cizq, siendo:
   // [C] = 1 / delta_t * [ M ] + alfa * [ K ]
   // [Caux] = 1 / delta_t * [ M ] - ( 1 - alfa ) * [ K ]

    int i, a;
    //int b, j, k;
    //int j;

    for ( i = 0; i < nodos_elem; i++){

        a = elems[ m ].conec[ i ]; // renglón

        C[ a ][ 0 ] += m_elem[ i ];
        Caux[ a ][ 0 ] += m_elem[ i ];

    }

     /*printf ( "\n\nMATRIZ C \n\n" );
    for ( i = 0; i < nodos; i++){
        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){
            printf ( "%lf\t", C[ i ][ j ] );
            if ( j == tabla_ady[ i ].cols -1 ) printf( "\n");
        }
    }*/



    /*for ( j = 0; j < nodos_elem; j++ ){

            a = elems[ m ].conec[ i ]; // renglón
            b = elems[ m ].conec[ j ]; // columna

            for ( k = 0; k < tabla_ady[ a ].cols; k++){

                if ( tabla_ady[ a ].ind[ k ] == b ){
                    C[ a ][ k ] += m_elem[ i ][ j ];
                    Caux[ a ][ k ] =  m_elem[ i ][ j ] + Caux[ a ][ k ];
                    //printf ( "\nK [ %d ][ %d ] = %lf", i, j, K_elem[ i ][ j ]);
                    //printf ( "\ntabla_ady[ %d ].val[ %d ] = %lf\n", a, b, tabla_ady[ a ].val[ b ] );
                }
            }
        }
    }*/
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double **mem_C ( struct tipo_2 *tabla_ady, int nodos, char *nombre )
{
	/** Descripción **/
	// Pide memoria para las matrices C

	int i;
	double **k;

	k = ( double ** ) calloc ( nodos, sizeof ( double * ) );
	if ( k == NULL ){
		printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
		exit ( 1 );
	}
	for ( i = 0; i < nodos; i++ ){
		k[i] = ( double * ) calloc ( tabla_ady[ i ].cols, sizeof ( double ) );
		if ( k == NULL ){
			printf ( "\nMemoria insuficiente para almacenar matriz %s", nombre );
			exit ( 1 );
		}
	}

	return k;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void mult_MKT ( struct tipo_2 *A, double **val, double *x, double *b, int n )
{
    /** Descripción **/
    // Realiza la multiplicación:
    // {b2} = [Cder] * {T(t-1)}

    int i, j;
    double suma;

    inicializa_vector ( b, n );

    for ( i = 0; i < n; i++ ){
        suma = 0.0;
        for ( j = 0; j < A[ i ].cols; j++ ){
            suma +=  val[ i ][ j ] * x[ A[ i ].ind[ j ] ];
        }
        b[ i ] = suma;
    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inicia_C ( struct tipo_2 *tabla_ady, double **A, int n )
{
    /** Descripción **/
    // Los valores de tabla_ady[].val[] son inicializados nuevamente al valor contenido en Cizq

    int i, j;

    for ( i = 0; i < n; i++ ){
        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){
            tabla_ady[ i ].val[ j ] = A[ i ][ j ];
        }
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


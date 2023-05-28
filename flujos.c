//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void matriz_B ( double **B, double **jacobiano, double *vect_aux, double *vec_aux2, double **pg,
                int dimension, int nodos_elem, int k );
void dif_matriz ( double **difusion, double *prop, int dimension );
void q ( double **q_pg, double *prop, double **pg, double **wpg, double **jacobiano, double **B,
         double **derin, double *vect_aux, double *vect_aux2, double *vect_aux3, double *T_elem,
         int dimension, int tipo_elem, int nodos_elem, int pos, double **difusion, double **coord,
         double *det_aux );
void masas_flujos_pg ( double **pg, double **wpg, double *M, double *f_pg, double f,double nodos_elem,
                       int pg_elem, int tipo_elem, double *det_aux );
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void flujos_1D ( struct elem *elems, double *T, double **q_nodos, int nodos, int n_elem )
{

    double K_elem[ 2 ][ 2 ], T_elem[ 2 ];
    int a, b, m;
    //int i;
    double Q, K, lx, ly, l;

    for ( m = 0; m < n_elem; m++ ){

        a = elems[ m ].conec[ 0 ];
        b = elems[ m ].conec[ 1 ];
        K = elems[ m ].prop[ 0 ];
        Q = elems[ m ].prop[ 1 ];

        lx = coord[ b ][ 0 ] - coord[ a ][ 0 ];
        ly = coord[ b ][ 1 ] - coord[ a ][ 1 ];
        l = sqrt( lx * lx + ly * ly );

        K_elem[ 0 ][ 0 ] =  K / l;
        K_elem[ 0 ][ 1 ] = -K / l;
        K_elem[ 1 ][ 0 ] = -K / l;
        K_elem[ 1 ][ 1 ] =  K / l;

        T_elem[ 0 ] = T[ a ];
        T_elem[ 1 ] = T[ b ];

        q_nodos[ a ][ 0 ] =  K_elem[ 0 ][ 0 ] * T_elem[ 0 ] + K_elem[ 0 ][ 1 ] * T_elem[ 1 ]
                             - Q * l / 2.0;
        if ( m == n_elem - 1 ){
            q_nodos[ b ][ 0 ] = -( K_elem[ 1 ][ 0 ] * T_elem[ 0 ] + K_elem[ 1 ][ 1 ] * T_elem[ 1 ]
                                - Q * l / 2.0 );

        }

    }

    /*printf ( "\n\nFlujos en nodos 1D" );
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%lf", q_nodos[ i ][ 0 ] );
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void suavizado_flujos ( struct elem *a, double *T, int dimension, int tipo_elem, int n_elem,
                        int pg_elem, int nodos_elem, double **coord, int nodos, double **q_pg,
                        double **q_nodos, double **q_pg2, double **error_PG )
{
    /** Descripción **/
    // - Calcula los flujos  qx y qy de los elementos en los P.G.
    // - Calcula los flujos  qx y qy de los elementos en los nodos

    /** Variables de entrada **/

    /** Variables de salida **/

    /** Variables locales **/
    int m, pos, i, j, k, aux;
    double **difusion, **pg, **wpg, **jacobiano, **B, **derin, **coord_aux;
    double *T_elem, *prop, *vect_aux, *vect_aux2, *vect_aux3, *M, *f_pg, *det_aux;
    double f, p, n, Naux, *MASAS, *sigma_x, *sigma_y, *sigma_z;
    double sum_x, sum_y, sum_z, mag_error, mag_flujo, fx, fy, fz;
    double error_x, error_y, error_z, error_xy, error_xyz, b, c, d, e, x, y;
    char pctge;

    pctge = '%';

    // Pide memoria
    difusion =  mem_matriz ( dimension, dimension,  "difusion" );
    prop = mem_vector ( dimension, "de prop de elementos");
    pg = mem_matriz ( pg_elem, dimension,  "de puntos de Gauss" );
	wpg = mem_matriz ( pg_elem, dimension,  "de pesos de Gauss" );
	jacobiano = mem_matriz ( dimension, dimension,  "del Jacobiano" );
	B =  mem_matriz ( dimension, nodos_elem,  "B" );
	derin = mem_matriz ( dimension, nodos_elem,  "de derivadas naturales" );
	coord_aux = mem_matriz ( nodos_elem, dimension, "de coord_aux de nodos locales para elementos");
	T_elem = mem_vector ( nodos_elem, "temperatura para calc. qx y qy" );
	vect_aux = mem_vector ( dimension, "vectaux" );
	vect_aux2 = mem_vector ( dimension, "vectaux2" );
	vect_aux3 = mem_vector ( dimension, "vectaux3" );
	M = mem_vector ( nodos_elem, "matriz de masas");
	f_pg = mem_vector ( nodos_elem, "flujos en P.G.");
	det_aux = mem_vector ( pg_elem, "determinantes");
	MASAS = mem_vector ( nodos, "matriz de MASAS global" );
	sigma_x = mem_vector ( nodos, "flujos en P.G. en x" );
	sigma_y = mem_vector ( nodos, "flujos en P.G. en y" );
	sigma_z = mem_vector ( nodos, "flujos en P.G. en z" );

    // Determina los P.G. y sus respectivos pesos
	puntos_gauss ( pg, wpg, tipo_elem );

	inicializa_matriz ( q_pg, pg_elem * n_elem, dimension );
	inicializa_matriz ( q_pg2, pg_elem * n_elem, dimension );
	inicializa_matriz ( q_nodos, nodos, dimension );
	inicializa_vector ( sigma_x, nodos );
	inicializa_vector ( sigma_y, nodos );
	inicializa_vector ( sigma_z, nodos );
	inicializa_vector ( MASAS, nodos );

    for ( m = 0; m < n_elem; m++ ){

        //Obtienen las propiedades del elemento
        if ( dimension == 1 ){
            printf( "\nEntra a 1D, para asignar propiedades." );
            exit( 1 );
        }
        else if ( dimension == 2 ){
            prop[ 0 ] = a[ m ].prop[ 0 ];
            prop[ 1 ] = a[ m ].prop[ 1 ];
        }
        else{
            prop[ 0 ] = a[ m ].prop[ 0 ];
            prop[ 1 ] = a[ m ].prop[ 1 ];
            prop[ 2 ] = a[ m ].prop[ 2 ];
        }

        //printf( "\nprop");
        //printf( "\nprop[ 0 ] = %lf", prop[ 0 ] );
        //printf( "\nprop[ 1 ] = %lf", prop[ 1 ] );

        // Temperaturas correspondientes al elemento
        //printf( "\n");
        for ( i = 0; i < nodos_elem; i++ ){
            T_elem[ i ] = T[ a[m].conec[i] ];
            //printf( "\nTemperaturas");
            //printf( "\n%lf", T_elem[ i ] );
        }


        // Obtiene las coodenadas de  los nodos del elemento
        for ( i = 0; i < nodos_elem; i++){
            for ( j = 0; j < dimension; j++){
                aux = a[ m ].conec[ i ];
                coord_aux[ i ][ j ] = coord[ aux ][ j ];
            }
        }

        // Posición de los P.G. de acuerdo al elemento analizado
        pos = m * pg_elem;

        q ( q_pg, prop, pg, wpg, jacobiano, B, derin, vect_aux, vect_aux2, vect_aux3, T_elem,
            dimension, tipo_elem, nodos_elem, pos, difusion, coord_aux, det_aux );

        /*printf ( "\n\n" );
        printf ( "\nDet. para cada P.G." );
        for ( i = 0; i < pg_elem; i++ ){
            printf ( "\n%d\t%lf", i, det_aux[ i ] );
        }*/

    }

     for ( m = 0; m < n_elem; m++ ){
         pos = m * pg_elem;
        // Se genera la matriz de masas y los respectivos vectores y se ensamblan.
        for ( k = 0; k < pg_elem; k++ ){
            for ( j = 0; j < dimension; j++ ){
                f = q_pg[ pos + k ][ j ];
                masas_flujos_pg ( pg, wpg, M, f_pg, f, nodos_elem, pg_elem, tipo_elem, det_aux );

                /*printf ( "\nf_pg ");
                for ( i = 0; i < nodos_elem; i++ ){
                    printf( "\n%lf", f_pg[ i ] );
                }
                */
                /*printf ( "\nMASAS ");
                for ( i = 0; i < nodos_elem; i++ ){
                    printf( "\n%lf", M[ i ] );
                }*/

                // Ensambla flujos

                if ( dimension == 1 ){
                    printf( "\nEntra a 1D para ensamblar flujos de suavizado.");
                    exit( 1 );
                }
                else if ( dimension == 2 ){
                    for ( i = 0; i < nodos_elem; i++ ){
                        if ( j == 0 ){
                            sigma_x[ a[m].conec[i] ] += f_pg[ i ];
                        }
                        if ( j == 1 ){
                            sigma_y[ a[m].conec[i] ] += f_pg[ i ];
                        }
                    }
                }
                else{
                    for ( i = 0; i < nodos_elem; i++ ){
                        if ( j == 0 ){
                            sigma_x[ a[m].conec[i] ] += f_pg[ i ];
                        }
                        if ( j == 1 ){
                            sigma_y[ a[m].conec[i] ] += f_pg[ i ];
                        }
                        if ( j == 2 ){
                            sigma_z[ a[m].conec[i] ] += f_pg[ i ];
                        }
                    }
                }

                /*for ( i = 0; i < nodos_elem; i++ ){
                    q_nodos[ a[m].conec[i] ][ j ] += f_pg[ i ] ;
                }*/

            }
            // Ensambla matriz de masas
            for ( i = 0; i < nodos_elem; i++ ){
                MASAS[ a[m].conec[i] ] += M[ i ];
            }
        }
    }

    /*printf ( "\n\n" );
    printf ( "\nMatriz de masas" );
    k = 0;
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%d\t%lf", i, MASAS[ i ] );
    }*/

    // Calcula los flujos en nodos
    if ( dimension == 1 ){
        printf( "\nEntra a 1D para calcular flujos en nodos.");
        exit( 1 );
    }
    else if ( dimension == 2 ){
        for ( i = 0; i < nodos; i++ ){
            q_nodos[ i ][ 0 ] = sigma_x[ i ] / MASAS[ i ];
            q_nodos[ i ][ 1 ] = sigma_y[ i ] / MASAS[ i ];
        }
    }
    else{
        for ( i = 0; i < nodos; i++ ){
            q_nodos[ i ][ 0 ] = sigma_x[ i ] / MASAS[ i ];
            q_nodos[ i ][ 1 ] = sigma_y[ i ] / MASAS[ i ];
            q_nodos[ i ][ 2 ] = sigma_z[ i ] / MASAS[ i ];
        }
    }

    /*printf ( "\n\n" );
    printf ( "\nFlujos de P.G. en nodos" );
    k = 0;
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%d\t%lf\t%lf", i, sigma_x[ i ], sigma_y[ i ] );
    }*/

    /*
    printf ( "\n\n" );
    printf ( "\nFlujos pasados a nodos" );
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%d\t%lf\t%lf", i, q_nodos[ i ][ 0 ], q_nodos[ i ][ 1 ] );
    }

    printf ( "\n\n" );
    printf ( "\nFlujos en los puntos de Gauss antes" );
    k = 0;
    for ( i = 0; i < n_elem; i++ ){
        pos = i * pg_elem;
        printf ( "\nElemento %d", i );
        for ( j = 0; j < pg_elem; j++){
            printf ( "\n%d\t%lf\t%lf", j, q_pg[ pos + j ][ 0 ], q_pg[ pos + j ][ 1 ] );
        }
    }*/

    // Vuelve a calcular los flujos en los P.G.
    for ( m = 0; m < n_elem; m++ ){

        pos = m * pg_elem;

        if ( dimension == 1 ){
            printf( "\nEntra a 1D para calc. flujos en P.G. a partir de nodos");
            exit( 1 );
        }
        else if ( dimension == 2 ){
            switch ( tipo_elem ){
                case 1:
                {
                    p = pg[ 0 ][ 0 ];
                    n = pg[ 0 ][ 1 ];
                    for ( j = 0; j < nodos_elem; j++ ){
                        Naux = N ( tipo_elem, j, dimension, p, n, 0.0 );
                        q_pg2[ pos ][ 0 ] += Naux * q_nodos[ a[m].conec[j] ][ 0 ];
                        q_pg2[ pos ][ 1 ] += Naux * q_nodos[ a[m].conec[j] ][ 1 ];
                    }
                    break;
                }
                case 2:
                {
                    for ( i = 0; i < pg_elem; i++){
                        p = pg[ i ][ 0 ];
                        n = pg[ i ][ 1 ];
                        for ( j = 0; j < nodos_elem; j++ ){
                            Naux = N ( tipo_elem, j, dimension, p, n, 0.0 );
                            q_pg2[ pos + i ][ 0 ] += Naux * q_nodos[ a[m].conec[j] ][ 0 ];
                            q_pg2[ pos + i ][ 1 ] += Naux * q_nodos[ a[m].conec[j] ][ 1 ];
                        }
                    }
                    break;
                }
                default:
                {
                    printf ( "Error: no se pudo calcular flujos en P.G. a partir de flujos en nodos (2D).");
                }
            }
        }
        else{
            switch ( tipo_elem ){
                case 3:
                {
                   for ( i = 0; i < pg_elem; i++){
                        p = pg[ i ][ 0 ];
                        n = pg[ i ][ 1 ];
                        c = pg[ i ][ 2 ];
                        for ( j = 0; j < nodos_elem; j++ ){
                            Naux = N ( tipo_elem, j, dimension, p, n, c );
                            q_pg2[ pos + i ][ 0 ] += Naux * q_nodos[ a[m].conec[j] ][ 0 ];
                            q_pg2[ pos + i ][ 1 ] += Naux * q_nodos[ a[m].conec[j] ][ 1 ];
                            q_pg2[ pos + i ][ 2 ] += Naux * q_nodos[ a[m].conec[j] ][ 2 ];
                        }
                   }
                    break;
                }
                case 4:
                {
                    for ( i = 0; i < pg_elem; i++){
                        p = pg[ i ][ 0 ];
                        n = pg[ i ][ 1 ];
                        for ( j = 0; j < nodos_elem; j++ ){
                            Naux = N ( tipo_elem, j, dimension, p, n, 0.0 );
                            q_pg2[ pos + i ][ 0 ] += Naux * q_nodos[ a[m].conec[j] ][ 0 ];
                            q_pg2[ pos + i ][ 1 ] += Naux * q_nodos[ a[m].conec[j] ][ 1 ];
                            q_pg2[ pos + i ][ 2 ] += Naux * q_nodos[ a[m].conec[j] ][ 2 ];
                        }
                    }
                    break;
                }
                default:
                {
                    printf ( "Error: no se pudo calcular flujos en P.G. a partir de flujos en nodos (3D).");
                }
            }
        }
    }

    /*printf ( "\n\n" );
    printf ( "\nFlujos en los puntos de Gauss nuevos" );
    k = 0;
    for ( i = 0; i < n_elem; i++ ){
        pos = i * pg_elem;
        printf ( "\nElemento %d", i );
        for ( j = 0; j < pg_elem; j++){
            printf ( "\n%d\t%lf\t%lf", j, q_pg2[ pos + j ][ 0 ], q_pg2[ pos + j ][ 1 ] );
        }
    }*/

    // Obtiene el error
    sum_x = 0.0;
    sum_y = 0.0;
    sum_z = 0.0;

    fx = 0.0;
    fy = 0.0;
    fz = 0.0;

    for ( i = 0; i < pg_elem * n_elem; i++ ){

        if ( dimension == 2 ){

            b = q_pg[ i ][ 0 ] - q_pg2[ i ][ 0 ];
            c = q_pg[ i ][ 1 ] - q_pg2[ i ][ 1 ];

            error_PG[ i ][ 0 ] = b;
            error_PG[ i ][ 1 ] = c;

            d = q_pg2[ i ][ 0 ];
            e = q_pg2[ i ][ 1 ];

            sum_x += fabs( b );
            sum_y += fabs( c );

            fx += fabs( d );
            fy += fabs( e );

        }

        if ( dimension == 3 ){

            b = q_pg[ i ][ 0 ] - q_pg2[ i ][ 0 ];
            c = q_pg[ i ][ 1 ] - q_pg2[ i ][ 1 ];
            x = q_pg[ i ][ 2 ] - q_pg2[ i ][ 2 ];

            error_PG[ i ][ 0 ] = b;
            error_PG[ i ][ 1 ] = c;
            error_PG[ i ][ 2 ] = x;

            d = q_pg2[ i ][ 0 ];
            e = q_pg2[ i ][ 1 ];
            y = q_pg2[ i ][ 2 ];

            sum_x += fabs ( b );
            sum_y += fabs ( c );
            sum_z += fabs ( x );

            fx += fabs( d );
            fy += fabs( e );
            fz += fabs( y );



        }

    }

    mag_error = sqrt ( sum_x * sum_x + sum_y * sum_y + sum_z * sum_z );
    mag_flujo = sqrt ( fx * fx + fy * fy + fz * fz );

    if ( dimension == 2 ){
        error_x = ( sum_x / mag_error ) * ( sum_x / mag_error ) / mag_error * 100.0 / mag_flujo;
        error_y = ( sum_y / mag_error ) * ( sum_y / mag_error ) / mag_error * 100.0 / mag_flujo;
        error_xy = mag_error / mag_flujo * 100.0;
        printf( "\nError_x  ( %c ) = %lf", pctge, error_x );
        printf( "\nError_y  ( %c ) = %lf", pctge, error_y );
        printf( "\nError_xy ( %c ) = %lf\n", pctge, error_xy );
    }

    if ( dimension == 3 ){
        error_x = ( sum_x / mag_error ) * ( sum_x / mag_error ) * mag_error * 100.0 / mag_flujo;
        error_y = ( sum_y / mag_error ) * ( sum_y / mag_error ) * mag_error * 100.0 / mag_flujo;
        error_z = ( sum_z / mag_error ) * ( sum_z / mag_error ) * mag_error * 100.0 / mag_flujo;
        error_xyz = mag_error / mag_flujo * 100.0;
        //error_xyz = suma5 / suma6 * 100.0;
        printf( "\nError_x  ( %c ) = %lf", pctge, error_x );
        printf( "\nError_y  ( %c ) = %lf", pctge, error_y );
        printf( "\nError_z  ( %c ) = %lf", pctge, error_z );
        printf( "\nError_xyz  ( %c ) = %lf\n", pctge, error_xyz );
        //printf( "\nError_xyz ( %c ) = %lf\n\n", pctge, error_xyz );
    }

    libera_matriz ( difusion, dimension );
    libera_vector ( prop );
    libera_matriz ( pg, pg_elem );
    libera_matriz ( wpg, pg_elem );
    libera_matriz ( jacobiano, dimension );
    libera_matriz ( B, dimension );
    libera_matriz ( derin, dimension );
    libera_matriz ( coord_aux, nodos_elem );
    libera_vector ( T_elem );
    libera_vector ( vect_aux );
    libera_vector ( vect_aux2 );
    libera_vector ( vect_aux3 );
    libera_vector ( M );
    libera_vector ( f_pg );
    libera_vector ( det_aux );
    libera_vector ( MASAS );
    libera_vector ( sigma_x );
    libera_vector ( sigma_y );
    libera_vector ( sigma_z );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void q ( double **q_pg, double *prop, double **pg, double **wpg, double **jacobiano, double **B,
         double **derin, double *vect_aux, double *vect_aux2, double *vect_aux3, double *T_elem,
         int dimension, int tipo_elem, int nodos_elem, int pos, double **difusion, double **coord,
         double *det_aux )
{
    /** Descripción **/
    // Calcula qx y qy en los puntos especificados

    /** Variables de entrada **/

    /** Variables de salida **/

    /** Variables locales **/

    int i, k;
    double det_jacobiano, **inv_jacobiano;
    double p = 0.0, n = 0.0, c = 0.0, wp = 0.0, wn = 0.0, wc = 0.0;
    //int a, b;

    inv_jacobiano = mem_matriz ( dimension, dimension,  "del inv_jacobiano" );

    dif_matriz ( difusion, prop, dimension );//nueva

    for ( k = 0; k < pg_elem; k++ ){

        // Coordenadas de P.G. y sus respectivos pesos
	    if ( dimension == 1 ){
            printf( "\nEntra a 1D, para asignar p,n,c..." );
            exit( 1 );
	    }
        else if ( dimension == 2 ){
            p  = pg[ k ][ 0 ];
            n  = pg[ k ][ 1 ];
            wp = wpg[ k ][ 0 ];
            wn = wpg[ k ][ 1 ];
	    }
	    else if ( dimension == 3 ){
            p  = pg[ k ][ 0 ];
            n  = pg[ k ][ 1 ];
            c  = pg[ k ][ 2 ];
            wp = wpg[ k ][ 0 ];
            wn = wpg[ k ][ 1 ];
            wc = wpg[ k ][ 2 ];
	    }

		inicializa_matriz ( jacobiano, dimension, dimension );
		inicializa_matriz ( inv_jacobiano, dimension, dimension );
		inicializa_matriz ( B, dimension, nodos_elem );

		if ( dimension == 1){
            printf( "\nEntra a 1D para determinar kelem.");
            exit( 1 );
        }
        else if ( dimension == 2 ){
            fjacobiano ( jacobiano, coord, p, n, 0.0, nodos_elem, tipo_elem, dimension, derin );
        }
        else{
            fjacobiano ( jacobiano, coord, p, n, c, nodos_elem, tipo_elem, dimension, derin );
        }

		det_jacobiano = determinante ( jacobiano, dimension );

		det_aux[ k ] = det_jacobiano; // Para después emplearlos en la integración de la matriz de masas

		//printf ( "\ndet_jacobiano = %lf", det_jacobiano );

		inversa ( jacobiano, inv_jacobiano, dimension, det_jacobiano );

		/*printf ( "\ninversa_jacobiano\n");
		for ( a = 0; a < dimension; a++ ){
			for ( b = 0; b < dimension; b++ ){
				printf ( "%lf\t", inversa_jacobiano[a][b] );
				if ( b == dimension - 1 ) printf ( "\n" );
			}
		}
		printf ( "\n\n" );*/

        matriz_B ( B, inv_jacobiano, vect_aux, vect_aux2, pg, dimension, nodos_elem, k ); //nueva

		/*printf ( "\nMatriz B\n " );
		for ( a = 0; a <  dimension; a++ ){
			for ( b = 0; b < nodos_elem; b++ ){
				printf ( "%lf\t", B[a][b] );
			if ( b == nodos_elem - 1 ) printf ( "\n" );
			}
		}*/
        inicializa_vector ( vect_aux3, dimension );
		mult_mat_vector ( vect_aux3, B, T_elem, dimension, nodos_elem, nodos_elem );

		/*printf ( "\nB * T " );
		printf( "\n%lf", vect_aux3[ 0 ] );
		printf( "\n%lf", vect_aux3[ 1 ] );
		printf ( "\n" );*/

        inicializa_vector ( vect_aux, dimension );
		mult_mat_vector ( vect_aux, difusion, vect_aux3, dimension, dimension, dimension );

        /*printf ( "\nMatris difusion\n " );
		for ( a = 0; a <  dimension; a++ ){
			for ( b = 0; b < dimension; b++ ){
				printf ( "%lf\t", difusion[a][b] );
			if ( b == dimension - 1 ) printf ( "\n" );
			}
		}*/

		/*printf ( "\nK *B * T " );
		printf( "\n%lf", vect_aux[ 0 ] );
		printf( "\n%lf", vect_aux[ 1 ] );
		printf ( "\n" );*/

        for ( i = 0; i < dimension; i++ ){
            q_pg[ pos + k ][ i ] = - vect_aux[ i ];
        }

    }

    libera_matriz ( inv_jacobiano, dimension );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void dif_matriz ( double **difusion, double *prop, int dimension )
{
    /** Descripción **/
    // Determina la matriz de propiedades de difusión
    // Es decir la que contiene a Kx y Ky

    int i;

    inicializa_matriz ( difusion, dimension, dimension );

    for ( i = 0; i < dimension; i ++ ){
		difusion[ i ][ i ] = prop[ i ];
	}

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void matriz_B ( double **B, double **inv_jacobiano, double *vect_aux, double *vect_aux2, double **pg,
                int dimension, int nodos_elem, int k )
{
    /** Descripción **/
    // Calcula la matriz B para cada punto de Gauss

    int i, j;
    //int a, b;

    for ( j = 0; j < nodos_elem; j++ ){

        inicializa_vector ( vect_aux, dimension );
        inicializa_vector ( vect_aux2, dimension );

        deriv_Ni ( vect_aux, pg, tipo_elem, dimension, j, k );

        /*printf ( "\n\nvect_aux" );
        for ( a = 0; a <  dimension; a++ ){
            printf ( "\n%lf", vect_aux[a] );
        }*/

        mult_mat_vector ( vect_aux2, inv_jacobiano, vect_aux, dimension, dimension, dimension );

        /*printf ( "\n\nvect_aux2" );
        for ( a = 0; a <  dimension; a++ ){
            printf ( "\n%lf", vect_aux2[a] );
        }*/

        for ( i = 0; i < dimension; i++ ){
            B[i][j] = vect_aux2[i];
        }

    }

    /*printf ( "\nMatriz B\n " );
    for ( a = 0; a <  dimension; a++ ){
        for ( b = 0; b < nodos_elem; b++ ){
            printf ( "%lf\t", B[a][b] );
        if ( b == nodos_elem - 1 ) printf ( "\n" );
        }
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void masas_flujos_pg ( double **pg, double **wpg, double *M, double *f_pg, double f,double nodos_elem,
                       int pg_elem, int tipo_elem, double *det_aux )
{
    /** Descripción **/
    // Calcula la matriz de masas

    /** Variables de entrada **/

    /** Variables de salida **/

    /**Variables locales **/
    double p = 0.0, n = 0.0, c = 0.0, wp = 0.0, wn = 0.0, wc = 0.0;
    double Nreng, Ncol, M_aux;
    int k, i, j;

    inicializa_vector ( M, nodos_elem );
    inicializa_vector ( f_pg, nodos_elem );

    for ( k = 0; k < pg_elem; k++ ){

        if ( dimension == 1 ){
            printf( "\nEntra a 1D, para asignar p,n,c..." );
            exit( 1 );
	    }
        else if ( dimension == 2 ){
            p  = pg[ k ][ 0 ];
            n  = pg[ k ][ 1 ];
            wp = wpg[ k ][ 0 ];
            wn = wpg[ k ][ 1 ];
	    }
	    else if ( dimension == 3 ){
            p  = pg[ k ][ 0 ];
            n  = pg[ k ][ 1 ];
            c  = pg[ k ][ 2 ];
            wp = wpg[ k ][ 0 ];
            wn = wpg[ k ][ 1 ];
            wc = wpg[ k ][ 2 ];
	    }

        //inicializa_vector ( M_aux, nodos_elem );

        for ( i = 0; i < nodos_elem; i++ ){
            if ( dimension == 1 ){
                printf( "\nEntra a 1D, para calc. Nren." );
                exit( 1 );
            }
            else if ( dimension == 2 ){
                Nreng = N ( tipo_elem, i, dimension, p, n, 0.0 );
            }
            else{
                Nreng = N ( tipo_elem, i, dimension, p, n, c );
            }

            M_aux = 0.0;

            for ( j = 0; j < nodos_elem; j++ ){

                if ( dimension == 1 ){
                    printf( "\nEntra a 1D, para calc. Ncol." );
                    exit( 1 );
                }
                else if ( dimension == 2 ){
                    Ncol = N ( tipo_elem, j, dimension, p, n, 0.0 );
                }
                else{
                    Ncol = N ( tipo_elem, j, dimension, p, n, c );
                }

                M_aux += Nreng * Ncol;

            }

            if ( dimension == 2 ){
                M[ i ] += M_aux * wp * wn * det_aux[ k ];
                f_pg[ i ] += f * Nreng * wp * wn * det_aux[ k ];
            }
            if ( dimension == 3 ){
                M[ i ] += M_aux * wp * wn * wc * det_aux[ k ];
                f_pg[ i ] += f * Nreng * wp * wn * wc * det_aux[ k ];
            }
        }

    }

    /*printf ( "\n\n" );
    printf ( "\nMatriz de masas" );
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%d\t%lf", i, M[ i ] );
    }*/

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double N ( int tipo_elem, int i, int dimension, double p, double n, double c )
{
    /** Descripción **/
    // Evalúa las funciones de forma en un P.G.

    double N;

    if ( dimension == 1 ){
            printf( "\nEntra a 1D, para asignar evaluar funcion N." );
            exit( 1 );
    }
    else if ( dimension == 2 ){
        switch ( tipo_elem ){
            case 1:
            {
                if ( i == 0 ){
                    N = 1 - p - n;
                }
                else if ( i == 1 ){
                    N = p;
                }
                else{
                    N = n;
                }
                break;
            }
            case 2:
            {
                if ( i == 0){
                    N = 1.0/4.0 * ( 1 - p ) * ( 1 - n );
                }
                else if ( i == 1 ){
                    N = 1.0/4.0 * ( 1 + p ) * ( 1 - n );
                }
                else if ( i == 2 ){
                    N = 1.0/4.0 * ( 1 + p ) * ( 1 + n );
                }
                else{
                    N = 1.0/4.0 * ( 1 - p ) * ( 1 + n );
                }
                break;
            }
            default:
            {
                puts ( "Error: calculo de matriz de masas. Caso no incluido" );
                exit ( 1 );
            }
        }
    }
    else{
        switch ( tipo_elem ){
            case 3:
            {
                if ( i == 0 ){
                    N = 1 - p - n - c;
                }
                else if ( i == 1 ){
                    N = p;
                }
                else if ( i == 2 ){
                    N = n;
                }
                else{
                    N = c;
                }
                break;
            }
            case 4:
            {
                if ( i == 0){
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 - n ) * ( 1 - c );
                }
                else if ( i == 1 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 - n ) * ( 1 - c );
                }
                else if ( i == 2 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 + n ) * ( 1 - c );
                }
                else if ( i == 3 ){
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 + n ) * ( 1 - c );
                }
                else if ( i == 4 ){
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 - n ) * ( 1 + c );
                }
                else if ( i == 5 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 - n ) * ( 1 + c );
                }
                else if ( i == 6 ){
                    N = 1.0/8.0 * ( 1 + p ) * ( 1 + n ) * ( 1 + c );
                }
                else{
                    N = 1.0/8.0 * ( 1 - p ) * ( 1 + n ) * ( 1 + c );
                }
                break;
            }
            default:
            {
                puts ( "Error: calculo de matriz de masas. Caso no incluido" );
                exit ( 1 );
            }
        }
    }

    return N;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

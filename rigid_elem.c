//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void calor_interno ( double *f_elem, double p, double n, double c, double wp, double wn, double wc,
                     double Q, double det_jacobiano, int tipo_elem, int nodos_elem, double l );
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void fk_elem ( double **K_elem, double *f_elem, double **coord, double *prop, int tipo_elem,
               int dimension, int nodos_elem, int pg_elem, double **pg, double **wpg,
               double **jacobiano, double *vect_aux, double *vect_aux2, double **B,
               double **BT, double **Baux, double **K_elemaux, double **derin,
               double **difusion, double *det, int m )
{
    /** Descripción **/
    // Calcula la matriz de rigidez de un elemento y su vector de flujos interno (Q).

    /** Variables de entrada **/

    /** Variables de salida **/

	int i, j, k;
    double Q, p = 0.0, n = 0.0, c = 0.0, wp = 0.0, wn = 0.0, wc = 0.0, det_jacobiano, **inv_jacobiano;
    double lx ,ly, l, K;
    //int a, b;
	//FILE *fp;


	if ( dimension == 1 ){
        lx = coord[ 1 ][ 0 ] - coord[ 0 ][ 0 ];
        ly = coord[ 1 ][ 1 ] - coord[ 0 ][ 1 ];
        l = sqrt( lx * lx + ly * ly );
        Q = prop[ 1 ];
        calor_interno ( f_elem, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Q, 0.0, tipo_elem, nodos_elem, l );
        K = prop[ 0 ];
        K_elem[ 0 ][ 0 ] =  1.0 * K / l;
        K_elem[ 0 ][ 1 ] = -1.0 * K / l;
        K_elem[ 1 ][ 0 ] = -1.0 * K / l;
        K_elem[ 1 ][ 1 ] =  1.0 * K / l;
	}
    else{

        inv_jacobiano = mem_matriz ( dimension, dimension,  "del inv_jacobiano" );

        inicializa_matriz ( difusion, dimension, dimension );

        for ( i = 0; i < dimension; i ++ ){
            difusion[i][i] = prop[i];
        }

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

            /** Cálculo de la matriz de rigidez elemental **/

            inicializa_matriz ( jacobiano, dimension, dimension );
            inicializa_matriz ( inv_jacobiano, dimension, dimension );
            inicializa_matriz ( B, dimension, nodos_elem );
            inicializa_matriz ( BT, nodos_elem, dimension );
            inicializa_matriz ( Baux, dimension, nodos_elem );
            inicializa_matriz ( K_elemaux, nodos_elem, nodos_elem );
            inicializa_matriz ( derin, dimension, nodos_elem );

            if ( dimension == 2 ){
                fjacobiano ( jacobiano, coord, p, n, 0.0, nodos_elem, tipo_elem, dimension, derin );
            }
            if ( dimension == 3){
                fjacobiano ( jacobiano, coord, p, n, c, nodos_elem, tipo_elem, dimension, derin );
            }

            det_jacobiano = determinante ( jacobiano, dimension );

            det[ m ] = det_jacobiano;

            //printf ( "\ndet_jacobiano = %lf\n\n", det_jacobiano );

            inversa ( jacobiano, inv_jacobiano, dimension, det_jacobiano );

            /*printf ( "\njacobiano\n");
            for ( a = 0; a < dimension; a++ ){
                for ( b = 0; b < dimension; b++ ){
                    printf ( "%lf\t", jacobiano[a][b] );
                    if ( b == dimension - 1 ) printf ( "\n" );
                }
            }
            printf ( "\n\n" );*/

            /*printf ( "\ninversa_jacobiano\n");
            for ( a = 0; a < dimension; a++ ){
                for ( b = 0; b < dimension; b++ ){
                    printf ( "%lf\t", inv_jacobiano[a][b] );
                    if ( b == dimension - 1 ) printf ( "\n" );
                }
            }
            printf ( "\n\n" );*/


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

            transponer_mat ( B, BT, dimension, nodos_elem );

            /*printf ( "\nBT\n " );
            for ( a = 0; a <  nodos_elem; a++ ){
                for ( b = 0; b < dimension; b++ ){
                    printf ( "%lf\t", BT[a][b] );
                if ( b == dimension - 1 ) printf ( "\n" );
                }
            }*/

            multiplica ( Baux, difusion, B, dimension, dimension, nodos_elem );

            /*printf ( "\nBaux\n " );
            for ( a = 0; a <  dimension; a++ ){
                for ( b = 0; b < nodos_elem; b++ ){
                    printf ( "%lf\t",Baux[a][b] );
                if ( b == nodos_elem - 1 ) printf ( "\n" );
                }
            }*/

            multiplica ( K_elemaux, BT, Baux, nodos_elem, dimension, nodos_elem );

            mult_mat_escalar ( K_elemaux, det_jacobiano, nodos_elem, nodos_elem );

            mult_mat_escalar ( K_elemaux, wp, nodos_elem, nodos_elem );
            mult_mat_escalar ( K_elemaux, wn, nodos_elem, nodos_elem );

            if ( dimension == 3 ){
                mult_mat_escalar ( K_elemaux, wc, nodos_elem, nodos_elem );
            }

            /*printf ( "\nMatriz k_elemaux\n " );
            for ( a = 0; a <  nodos_elem; a++ ){
                for ( b = 0; b < nodos_elem; b++ ){
                    printf ( "%lf\t", K_elemaux[a][b] );
                if ( b == nodos_elem - 1 ) printf ( "\n" );
                }
            }*/

            // Suma las matrices elementales corresponcientes a cada P.G.
            suma_mat ( K_elem, K_elemaux, nodos_elem, nodos_elem );

            if ( dimension == 1 ){
                printf( "\nEntra a 1D, para asignar calcular Q" );
                exit( 1 );
            }
            else if ( dimension == 2 ){
                Q = prop[ 2 ];
                calor_interno ( f_elem, p, n, 0.0, wp, wn, 0.0, Q, det_jacobiano, tipo_elem, nodos_elem, 0.0 );
            }
            else{
                 Q = prop[ 3 ];
                calor_interno ( f_elem, p, n, c, wp, wn, wc, Q, det_jacobiano, tipo_elem, nodos_elem, 0.0 );
            }

            // Imprime en pantalla f_elem para cada P.G.
            /*printf ( "\n\nf_elem " );
            for ( a = 0; a < nodos_elem; a++ ){
                printf( "\n%lf", f_elem[ a ] );
            }*/
        }

        libera_matriz ( inv_jacobiano, dimension );
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void fjacobiano ( double **jacobiano, double **coord, double p, double n, double c, int nodos_elem,
                  int tipo_elem, int dimension, double **derin )
{
	/** Descripción **/
	// Evalúa el jacobiano con las coordenadas globales de los nodos
	// y un punto de Gauss que se empleará para la integración.

	/** Variables de entrada **/
	// dimension	Dimensión del problema: 1D, 2D o 3D
	// nodos_elem	Número de nodos del elemento
	// coord[][]	Coordenadas globales de los nodos del elemento
	// p			Coordenada psi del P.G.
	// n			Coordenada eta del P.G.
	// tipo_elem	Tipo de elemento empleado

	/** Variables de salida **/
	// jacobiano	Jacobiano evaluado en P.G. ( p , n )

	/** Variables locales **/
	// derin		Derivadas naturales evaluadas en un punto de Gauss
	//			en P.G. ( p , n )

	//int a, b;
	//double **derin;


	//derin = mem_matriz ( dimension, nodos_elem,  "de derivadas naturales" );

	fderin ( derin, p, n, c, tipo_elem );
    //printf( "Despues de fderin \n");
	/*printf ( "\nderin\n");
	for ( a = 0; a < dimension; a++ ){
		for ( b = 0; b < nodos_elem; b++ ){
			printf ( "%lf\t", derin[a][b] );
			if ( b == nodos_elem - 1 ) printf ( "\n" );
		}
	}
	printf ( "\n\n" );*/

	multiplica ( jacobiano, derin, coord, dimension, nodos_elem, dimension );

	/*printf ( "\ncoord\n");
	for ( a = 0; a < nodos_elem; a++ ){
		for ( b = 0; b < dimension; b++ ){
			printf ( "%lf\t", coord[a][b] );
			if ( b == dimension - 1 ) printf ( "\n" );
		}
	}
	printf ( "\n\n" );*/


}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void fderin ( double **derin, double p, double n, double c, int tipo_elem )
{
	/** Descripción **/
	// Evalúa las derivadas naturales en los puntos de Gauss de acuerdo al tipo de
	// elemento que se esté empleando

	/** Variables de entrada */
	// p		Coordenada del P.G.
	// n		Coordenada del P.G.

	/**Variables de salida **/
	// derin	Derivadas naturales evaluadas en el P.G. ( p, n )

	/*
		n
		^
		|
		|
		|
		|--------------> p
	*/

    if ( dimension == 1 ){
        printf( "\nEntra a 1D, para derivadas naturales." );
        exit( 1 );
    }
    else if ( dimension == 2 ){
        switch ( tipo_elem ){
            case 1:
            {
                derin[0][0] = -1.0;
                derin[0][1] =  1.0;
                derin[0][2] =  0.0;

                derin[1][0] = -1.0;
                derin[1][1] =  0.0;
                derin[1][2] =  1.0;

                break;
            }
            case 2:
            {
                derin[0][0] = - 1.0/4.0 * ( 1.0 - n );
                derin[0][1] =   1.0/4.0 * ( 1.0 - n );
                derin[0][2] =   1.0/4.0 * ( 1.0 + n );
                derin[0][3] = - 1.0/4.0 * ( 1.0 + n );

                derin[1][0] = - 1.0/4.0 * ( 1.0 - p );
                derin[1][1] = - 1.0/4.0 * ( 1.0 + p );
                derin[1][2] =   1.0/4.0 * ( 1.0 + p );
                derin[1][3] =   1.0/4.0 * ( 1.0 - p );

                break;
            }
            default:
            {
                puts ( "Error: no se pudieron calcular las derivadas naturales. Caso no incluido (2D)." );
                exit ( 1 );
            }
        }
    }
    else{
        switch ( tipo_elem ){
            case 3:
            {
                derin[0][0] = -1.0;
                derin[0][1] =  1.0;
                derin[0][2] =  0.0;
                derin[0][3] =  0.0;

                derin[1][0] = -1.0;
                derin[1][1] =  0.0;
                derin[1][2] =  1.0;
                derin[1][3] =  0.0;

                derin[2][0] = -1.0;
                derin[2][1] =  0.0;
                derin[2][2] =  0.0;
                derin[2][3] =  1.0;

                break;
            }
            case 4:
            {
                derin[0][0] = -1.0/8.0 * ( 1 - n ) * ( 1 - c );
                derin[0][1] =  1.0/8.0 * ( 1 - n ) * ( 1 - c );
                derin[0][2] =  1.0/8.0 * ( 1 + n ) * ( 1 - c );
                derin[0][3] = -1.0/8.0 * ( 1 + n ) * ( 1 - c );
                derin[0][4] = -1.0/8.0 * ( 1 - n ) * ( 1 + c );
                derin[0][5] =  1.0/8.0 * ( 1 - n ) * ( 1 + c );
                derin[0][6] =  1.0/8.0 * ( 1 + n ) * ( 1 + c );
                derin[0][7] = -1.0/8.0 * ( 1 + n ) * ( 1 + c );

                derin[1][0] = -1.0/8.0 * ( 1 - p ) * ( 1 - c );
                derin[1][1] = -1.0/8.0 * ( 1 + p ) * ( 1 - c );
                derin[1][2] =  1.0/8.0 * ( 1 + p ) * ( 1 - c );
                derin[1][3] =  1.0/8.0 * ( 1 - p ) * ( 1 - c );
                derin[1][4] = -1.0/8.0 * ( 1 - p ) * ( 1 + c );
                derin[1][5] = -1.0/8.0 * ( 1 + p ) * ( 1 + c );
                derin[1][6] =  1.0/8.0 * ( 1 + p ) * ( 1 + c );
                derin[1][7] =  1.0/8.0 * ( 1 - p ) * ( 1 + c );

                derin[2][0] = -1.0/8.0 * ( 1 - p ) * ( 1 - n );
                derin[2][1] = -1.0/8.0 * ( 1 + p ) * ( 1 - n );
                derin[2][2] = -1.0/8.0 * ( 1 + p ) * ( 1 + n );
                derin[2][3] = -1.0/8.0 * ( 1 - p ) * ( 1 + n );
                derin[2][4] =  1.0/8.0 * ( 1 - p ) * ( 1 - n );
                derin[2][5] =  1.0/8.0 * ( 1 + p ) * ( 1 - n );
                derin[2][6] =  1.0/8.0 * ( 1 + p ) * ( 1 + n );
                derin[2][7] =  1.0/8.0 * ( 1 - p ) * ( 1 + n );

                break;
            }
            default:
            {
                puts ( "Error: no se pudieron calcular las derivadas naturales. Caso no incluido (3D)." );
                exit ( 1 );
            }
        }
    }
}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
double determinante ( double **A, int dimension )
{
	/** Descripción **/
	// Calcula el determinante de una matriz  [A] de tamaño 2 x 2

	double det;

	switch ( dimension ){
		case 2:
		{
			det = A[0][0]*A[1][1] - A[1][0]*A[0][1];
			break;
		}
		case 3:
		{
		    det =   A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*( A[1][0]*A[2][2]
                  - A[1][2]*A[2][0] ) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
            break;
		}
		default:
		{
			puts ( "Error: no se pudo calcular el determinante. Caso no incluido." );
			exit ( 1 );
		}
	}

    if ( det == 0 ){
        puts ( "Determinante del jacobiano = 0" );
        exit ( 1 );
    }

	return det;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void inversa ( double **A, double **inversa_A, int n, double det )
{
	/** Descripción **/
	// Calcula la inversa de matrices de 2 x 2

	/** Variables de entrada **/
	// A			Matriz a invertir
	// n			Tamaño de la matriz [A]

	/** Variables de salida **/
	// inversa_A	Inversa de la mat.[A]


	double cte;

	switch ( n ){
		case 2:
		{
			cte = A[0][0];
			inversa_A[0][0] =  1.0/det * A[1][1];
			inversa_A[0][1] = -1.0/det * A[0][1];
			inversa_A[1][0] = -1.0/det * A[1][0];
			inversa_A[1][1] =  1.0/det * cte;
			break;
		}
		case 3:
		{
            inversa_A[0][0] =  1.0/det * ( A[1][1]*A[2][2] - A[1][2]*A[2][1] );
            inversa_A[0][1] = -1.0/det * ( A[0][1]*A[2][2] - A[0][2]*A[2][1] );
            inversa_A[0][2] =  1.0/det * ( A[0][1]*A[1][2] - A[0][2]*A[1][1] );

            inversa_A[1][0] = -1.0/det * ( A[1][0]*A[2][2] - A[1][2]*A[2][0] );
            inversa_A[1][1] =  1.0/det * ( A[0][0]*A[2][2] - A[0][2]*A[2][0] );
            inversa_A[1][2] = -1.0/det * ( A[0][0]*A[1][2] - A[0][2]*A[1][0] );

            inversa_A[2][0] =  1.0/det * ( A[1][0]*A[2][1] - A[1][1]*A[2][0] );
            inversa_A[2][1] = -1.0/det * ( A[0][0]*A[2][1] - A[0][1]*A[2][0] );
            inversa_A[2][2] =  1.0/det * ( A[0][0]*A[1][1] - A[0][1]*A[1][0] );

		    break;
		}
		default:
		{
			puts ( "Error: no se pudo calcular la inversa. Caso no incluido." );
			exit ( 1 );
		}
	}
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void deriv_Ni( double *vect_aux, double **pg, int tipo_elem, int dimension, int j, int k )
{
    /** Descripción **/
    // Calcula las derivadas de la función de forma Ni respecto a las coordenadas del espacio
    // normalizado (p,n,c)


    if ( dimension == 1){
            printf( "\nEntra a 1D para determinar deriv_Ni.");
            exit( 1 );
    }
    else if ( dimension == 2){
        switch ( tipo_elem ){
            case 1:
            {
                if ( j == 0 ){
                    vect_aux[0] = -1.0;
                    vect_aux[1] = -1.0;
                }
                else if ( j == 1 ){
                    vect_aux[0] = 1.0;
                    vect_aux[1] = 0.0;
                }
                else{
                    vect_aux[0] = 0.0;
                    vect_aux[1] = 1.0;
                }
                break;
            }
            case 2:
            {
                if ( j == 0){
                    vect_aux[0] = - 1.0/4.0 * ( 1.0 - pg[k][1] );
                    vect_aux[1] = - 1.0/4.0 * ( 1.0 - pg[k][0] );
                }
                else if ( j == 1 ){
                    vect_aux[0] =   1.0/4.0 * ( 1.0 - pg[k][1] );
                    vect_aux[1] = - 1.0/4.0 * ( 1.0 + pg[k][0] );
                }
                else if ( j == 2 ){
                    vect_aux[0] =   1.0/4.0 * ( 1.0 + pg[k][1] );
                    vect_aux[1] =   1.0/4.0 * ( 1.0 + pg[k][0] );
                }
                else{
                    vect_aux[0] = - 1.0/4.0 * ( 1.0 + pg[k][1] );
                    vect_aux[1] =   1.0/4.0 * ( 1.0 - pg[k][0] );
                }
                break;
            }
            default:
            {
                puts ( "Error: no se pudo calcular el vec_aux para B. Caso no incluido" );
                exit ( 1 );
            }
        }
    }
    else{
         switch ( tipo_elem ){
            case 3:
            {
                if ( j == 0 ){
                    vect_aux[0] = -1.0;
                    vect_aux[1] = -1.0;
                    vect_aux[2] = -1.0;
                }
                else if ( j == 1 ){
                    vect_aux[0] = 1.0;
                    vect_aux[1] = 0.0;
                    vect_aux[2] = 0.0;
                }
                else if ( j == 2){
                    vect_aux[0] = 0.0;
                    vect_aux[1] = 1.0;
                    vect_aux[2] = 0.0;
                }
                else{
                    vect_aux[0] = 0.0;
                    vect_aux[1] = 0.0;
                    vect_aux[2] = 1.0;
                }
                break;
            }
            case 4:
            {
                if ( j == 0){
                    vect_aux[0] = -1.0/8.0 * ( 1.0 - pg[k][1] ) * ( 1.0 - pg[k][2] );
                    vect_aux[1] = -1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 - pg[k][2] );
                    vect_aux[2] = -1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 - pg[k][1] );
                }
                else if ( j == 1 ){
                    vect_aux[0] =  1.0/8.0 * ( 1.0 - pg[k][1] ) * ( 1.0 - pg[k][2] );
                    vect_aux[1] = -1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 - pg[k][2] );
                    vect_aux[2] = -1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 - pg[k][1] );
                }
                else if ( j == 2 ){
                    vect_aux[0] =  1.0/8.0 * ( 1.0 + pg[k][1] ) * ( 1.0 - pg[k][2] );
                    vect_aux[1] =  1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 - pg[k][2] );
                    vect_aux[2] = -1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 + pg[k][1] );
                }
                else if ( j == 3 ){
                    vect_aux[0] = -1.0/8.0 * ( 1.0 + pg[k][1] ) * ( 1.0 - pg[k][2] );
                    vect_aux[1] =  1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 - pg[k][2] );
                    vect_aux[2] = -1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 + pg[k][1] );
                }
                else if ( j == 4 ){
                    vect_aux[0] = -1.0/8.0 * ( 1.0 - pg[k][1] ) * ( 1.0 + pg[k][2] );
                    vect_aux[1] = -1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 + pg[k][2] );
                    vect_aux[2] =  1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 - pg[k][1] );
                }
                else if ( j == 5 ){
                    vect_aux[0] =  1.0/8.0 * ( 1.0 - pg[k][1] ) * ( 1.0 + pg[k][2] );
                    vect_aux[1] = -1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 + pg[k][2] );
                    vect_aux[2] =  1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 - pg[k][1] );
                }
                else if ( j == 6 ){
                    vect_aux[0] =  1.0/8.0 * ( 1.0 + pg[k][1] ) * ( 1.0 + pg[k][2] );
                    vect_aux[1] =  1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 + pg[k][2] );
                    vect_aux[2] =  1.0/8.0 * ( 1.0 + pg[k][0] ) * ( 1.0 + pg[k][1] );
                }
                else if ( j == 7 ){
                    vect_aux[0] = -1.0/8.0 * ( 1.0 + pg[k][1] ) * ( 1.0 + pg[k][2] );
                    vect_aux[1] =  1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 + pg[k][2] );
                    vect_aux[2] =  1.0/8.0 * ( 1.0 - pg[k][0] ) * ( 1.0 + pg[k][1] );
                }

                break;
            }
            default:
            {
                puts ( "Error: no se pudo calcular el vec_aux para B. Caso no incluido" );
                exit ( 1 );
            }
        }
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void calor_interno ( double *f_elem, double p, double n, double c, double wp, double wn, double wc,
                     double Q, double det_jacobiano, int tipo_elem, int nodos_elem, double l )
{
    /** Descripción **/
    // Calcula del vector de flujos debido a la generación interna de calor Q

    /** Variables de entrada **/
    // tipo_elem
    // nodos_elem
    // f_elem[]
    // pg[][]
    // wpg
    // det_jacobiano

    int j;

    //printf ( "\ndet_jacobiano = %lf", det_jacobiano );
    if ( dimension == 1 ){
        for ( j = 0; j < nodos_elem; j++ ){
            f_elem[ j ] = Q * l / 2.0;
        }
    }
    else if ( dimension == 2 ){
        switch ( tipo_elem ){
            case 1:
            {
                for ( j = 0; j < nodos_elem; j++){
                    if ( j == 0){
                        f_elem [ j ] += ( 1 - p - n ) * Q * det_jacobiano
                                        * wp * wn;
                    }
                    else if ( j == 1 ){
                        f_elem [ j ] += p * Q * det_jacobiano * wp * wn;
                    }
                    else{
                        f_elem [ j ] += n * Q * det_jacobiano * wp * wn;
                    }
                }
                break;
            }
            case 2:
            {
                for ( j = 0; j < nodos_elem; j++){
                    if ( j == 0){
                        f_elem[ j ] += 1.0 / 4.0 * ( 1 - p ) * ( 1 - n )
                                    * Q * det_jacobiano * wp * wn;
                    }
                    else if ( j == 1 ){
                        f_elem[ j ] += 1.0 / 4.0 * ( 1 + p ) * ( 1 - n )
                                    * Q * det_jacobiano * wp * wn;
                    }
                    else if ( j == 2 ){
                        f_elem[ j ] += 1.0 / 4.0 * ( 1 + p ) * ( 1 + n )
                                    * Q * det_jacobiano * wp * wn;
                    }
                    else{
                        f_elem[ j ] += 1.0 / 4.0 * ( 1 - p ) * ( 1 + n )
                                    * Q * det_jacobiano * wp * wn;
                    }
                }

                break;
            }
            default:
            {
                printf( "\nError: no se pudo calcular el vector de flujos debido a Q.");
            }
        }
    }
    else{

        switch ( tipo_elem ){
            case 3:
            {
                for ( j = 0; j < nodos_elem; j++){
                    if ( j == 0){
                        f_elem [ j ] += ( 1 - p - n - c ) * Q * det_jacobiano
                                        * wp * wn * wc;
                    }
                    else if ( j == 1 ){
                        f_elem [ j ] += p * Q * det_jacobiano * wp * wn * wc;
                    }
                    else if ( j == 2 ){
                        f_elem [ j ] += n * Q * det_jacobiano * wp * wn * wc;
                    }
                    else{
                        f_elem [ j ] += c * Q * det_jacobiano * wp * wn * wc;
                    }
                }
                break;
            }
            case 4:
            {
                for ( j = 0; j < nodos_elem; j++){
                    if ( j == 0){
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 - p ) * ( 1 - n ) * ( 1 - c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else if ( j == 1 ){
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 + p ) * ( 1 - n ) * ( 1 - c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else if ( j == 2 ){
                        f_elem[ j ] +=1.0 / 8.0 * ( 1 + p ) * ( 1 + n ) * ( 1 - c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else if ( j == 3 ){
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 - p ) * ( 1 + n ) * ( 1 - c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else if ( j == 4 ){
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 - p ) * ( 1 - n ) * ( 1 + c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else if ( j == 5 ){
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 + p ) * ( 1 - n ) * ( 1 + c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else if ( j == 6 ){
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 + p ) * ( 1 + n ) * ( 1 + c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                    else{
                        f_elem[ j ] += 1.0 / 8.0 * ( 1 - p ) * ( 1 + n ) * ( 1 + c )
                                        * Q * det_jacobiano * wp * wn * wc ;
                    }
                }
                break;
            }
            default:
            {
                printf( "\nError: no se pudo calcular el vector de flujos debido a Q.");
            }
        }
    }
}

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

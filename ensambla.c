//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void ady ( struct elem *elems, int n_elem, int nodos, int nodos_elem )
{
    /** Descripción **/
    // Obtiene la tabla de adyacencias a partir de la conectividad proporcionada
    // por GID. La conectividad es guardada en el vector conec de la estructura
    // elems.
    // IMPORTANTE: EL ARREGLO DE ESTRUCTURA elems[] CONTIENE LAS CONECTIVIDADES
    // TAL CUAL DADAS POR EL GID. NO SE LE HA RESTADO 1 A CADA NODO. Por alguna
    // razón a la hora de restar 1 a los subí­ndices variaban un poco las propie-
    // dades Kx y Ky del materia. Para implementar la resta solo hay que cambiar
    // los subÃ­ndices a y b de este archivo y activar el ciclo que resta.
    //
    // Lo anterior se resolvió restando 1 desde la lectura de los datos y ya to-
    // do esta iniccializado a cero. Esto es elems[].num, elems[].conec.
    //

    /** Variables de entrada **/
    // n_elem               Num. de elementos de la malla
    // nodos                Num. de nodos de la malla
    // nodos_elem           Num. de nodos del elemento que se empleó para mallar
    // elems[].conec        Contiene los nodos de los elementos
    // coord[][]            Guarda las coordenadas de todos los nodo

    /** Variables de salida **/
    // tabla_ady[].cols     Guarda el total de columnas de la matriz rala comprimida
    // tabla_ady[].ind[]    Contiene los í­ndices de los valores no cero de la matriz
    //                      de rigidez sin comprimir

    int i, j, k, h, a;
    int pos, no;
    //struct tipo_2 *tabla_ady;

    tabla_ady = mem_vector_tipo_2 ( nodos, "de tabla de adyacencias");
    for ( i = 0; i < nodos; i++ ){
        tabla_ady[ i ].ind = mem_vector_int ( 1, " vector de i­ndices de tabla_ady");
        tabla_ady[ i ].cols = 1;
    }

    for ( i = 0; i < n_elem; i++){
        for ( j = 0; j < nodos_elem; j++ ){
            for ( k = 0; k < nodos_elem; k++){
                if ( k == j ) continue;


                // Compara con todos los elementos existentes en la
                // tabla de adyacencias el valor que se quiere incor-
                // porar para no repetirlo
                pos = 1;
                no = 0;
                //a = elems[ i ].conec[ j ] - 1;
                a = elems[ i ].conec[ j ];
                for ( h = 1; h < tabla_ady[ a ].cols; h++ ){
                    if ( elems[ i ].conec[ k ] != tabla_ady[ a ].ind[ h ] ){
                        pos++;
                    }
                    else{
                        no = 1;
                    }
                }


                if ( no == 1 ) continue;

                pos++;

                // Si el valor no existe, pide memoria para el nuevo elemento
                redimensiona_vector_tipo_2 ( tabla_ady, a, pos, "adyacencias" );

                // Incorpora el nuevo elemento
                tabla_ady[ a ].ind[ pos - 1 ] = elems[i].conec[k];

                tabla_ady[ a ].cols = pos;

            }
        }
    }

    // Inicializa la primer columna de la tabla_ady que correponde al núm. de nodo
    for ( i = 0; i < nodos; i++){
        tabla_ady[ i ].ind[ 0 ] = i;
    }

    /*printf ( "\n\nTABLA ADYACENCIAS \n\n" );
    for ( i = 0; i < nodos; i++){
        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){
            printf ( "%d\t", tabla_ady[ i ].ind[ j ] );
            if ( j == tabla_ady[ i ].cols -1 ) printf( "\n");
        }
    }*/
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void rigid ( struct tipo_2 *tabla_ady, int tipo_elem, int nodos_elem, int pg_elem, int dimension,
             double *det )
{
    /** Descripción **/
    // Obtiene las matrices elementales y las ensambla en el formato comprimido

    /** Variables de entrada **/
    // tipo_elem            Tipo de elemento a emplear definido desde main
    // nodos_elem           Núm. de nodos del elemento que se empleó para mallar
    // pg_elem              Núm. de puntos de Gauss para integrar de acuerdo al
    //                      tipo de elemento
    // dimension            Dimensión
    // tabla_ady[].cols     Guarda el total de columnas de la matriz rala comprimida
    // tabla_ady[].ind[]    Contiene los índices de los valores no cero de la matriz
    //                      de rigidez sin comprimir

    /** Variables de salida **/
    // tabla_ady[].val[]    Contiene los valores de la matriz de rigidez en el for-
    //                      mato comprimido

    /** Variables locales **/
    // K_elem[][]           Matriz de rigidez elemental
    // coord_aux            Matriz con las coordenadas globales de los nodos de un
    //                      elemento
    // prop_aux[]           Vector con las propiedades de un elemento
    // f_elem[]             Vector de flujos del elemento causados por Q

    int i, j, m, k;
    int a, b;
    double **K_elem;
    double **coord_aux;
    double *prop_aux;
    double *f_elem;
    double **pg, **wpg;
	double **jacobiano;
	double *vect_aux, *vect_aux2;
	double **B, **BT, **Baux;
	double **K_elemaux, **derin;
	double **difusion;


    // Memoria para las variables locales
    if ( dimension == 1 ){
        K_elem =  mem_matriz ( nodos_elem, nodos_elem,  "K_elem" );
        coord_aux = mem_matriz ( nodos_elem, 2, "de coord_aux de nodos locales para elementos");
        prop_aux = mem_vector ( dimension + 1, "de prop_aux de elementos");
        f_elem = mem_vector ( nodos_elem, "vector de flujos de Q." );
    }
    else{
        K_elem =  mem_matriz ( nodos_elem, nodos_elem,  "K_elem" );
        coord_aux = mem_matriz ( nodos_elem, dimension, "de coord_aux de nodos locales para elementos");
        prop_aux = mem_vector ( dimension + 1, "de prop_aux de elementos");
        f_elem = mem_vector ( nodos_elem, "vector de flujos de Q." );
    }

	if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        pg = mem_matriz ( pg_elem, dimension,  "de puntos de Gauss" );
        wpg = mem_matriz ( pg_elem, dimension,  "de pesos de Gauss" );
        jacobiano = mem_matriz ( dimension, dimension,  "del Jacobiano" );
        vect_aux = mem_vector ( dimension, "vectaux" );
        vect_aux2 = mem_vector ( dimension, "vectaux2" );
        B =  mem_matriz ( dimension, nodos_elem,  "B" );
        BT =  mem_matriz ( nodos_elem, dimension,  "BT" );
        K_elemaux =  mem_matriz ( nodos_elem, nodos_elem,  "K_elemaux" );
        difusion =  mem_matriz ( dimension, dimension,  "difusion" );
        derin = mem_matriz ( dimension, nodos_elem,  "de derivadas naturales" );
        Baux =  mem_matriz ( dimension, nodos_elem,  "Baux" );
	}

    // Memoria para los valores de la matriz de rigidez
    // Inicializa matriz de rigidez y vector de flujos globales
    for ( i = 0; i < nodos; i++ ){
        tabla_ady[ i ].val = mem_vector ( tabla_ady[ i ].cols, "valores de rigideces" );
        inicializa_vector ( tabla_ady[ i ].val, tabla_ady[ i ].cols );
        tabla_ady[ i ].f = 0.0;
    }

    // Inicializa P.G. y pesos
    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        puntos_gauss ( pg, wpg, tipo_elem );
    }

    for ( m = 0; m < n_elem; m++){

        // Inicializa la matriz de rigidez
        inicializa_matriz ( K_elem, nodos_elem, nodos_elem );

        // Inicializa vector de flujos
        inicializa_vector ( f_elem, nodos_elem );

        if ( dimension == 1 ){
            for ( i = 0; i < nodos_elem; i++){
                for ( j = 0; j < 2; j++){
                    // Se usó 2 porque el dibujo de la varilla puede no ser
                    // una recta
                    //a = elems[ m ].conec[ i ] - 1;
                    a = elems[ m ].conec[ i ];
                    coord_aux[ i ][ j ] = coord[ a ][ j ];
                }
            }
        }
        else{
            // Obtiene las coordenadas de  los nodos del elemento
            for ( i = 0; i < nodos_elem; i++){
                for ( j = 0; j < dimension; j++){
                    //a = elems[ m ].conec[ i ] - 1;
                    a = elems[ m ].conec[ i ];
                    coord_aux[ i ][ j ] = coord[ a ][ j ];
                }
            }
        }

        /*printf( "\n\n ELEMENTO %d\n", m );
        printf( "\n Coordenadas\n");
        for ( i = 0; i < nodos_elem; i++){
            for ( j = 0; j < dimension; j++){
                printf ( "%lf\t", coord_aux[i][j]);
                if ( j == dimension - 1) printf( "\n");
            }
        }*/

        //Obtienen las propiedades del elemento
        if ( dimension == 1 ){
            prop_aux[ 0 ] = elems[ m ].prop[ 0 ];
            prop_aux[ 1 ] = elems[ m ].prop[ 1 ];
        }
        else if ( dimension == 2 ){
            prop_aux[ 0 ] = elems[ m ].prop[ 0 ];
            prop_aux[ 1 ] = elems[ m ].prop[ 1 ];
            prop_aux[ 2 ] = elems[ m ].prop[ 2 ];
        }
        else{
            prop_aux[ 0 ] = elems[ m ].prop[ 0 ];
            prop_aux[ 1 ] = elems[ m ].prop[ 1 ];
            prop_aux[ 2 ] = elems[ m ].prop[ 2 ];
            prop_aux[ 3 ] = elems[ m ].prop[ 3 ];
        }


        // Evalua la matriz de rigidez elemental para el elemento m
        fk_elem ( K_elem, f_elem, coord_aux, prop_aux, tipo_elem, dimension, nodos_elem, pg_elem,
                  pg, wpg, jacobiano, vect_aux, vect_aux2, B, BT, Baux, K_elemaux, derin, difusion,
                  det, m );

        /*printf ( "\nMatriz k_elem\n " );
        for ( a = 0; a <  nodos_elem; a++ ){
            for ( b = 0; b < nodos_elem; b++ ){
                printf ( "%lf\t", K_elem[a][b] );
            if ( b == nodos_elem - 1 ) printf ( "\n" );
            }
        }

        printf ( "\nf_elem " );
        for ( a = 0; a <  nodos_elem; a++ ){
            printf ( "\n%lf", f_elem[ a ] );

        }*/

        // Ensambla las matrices elementales
        for ( i = 0; i < nodos_elem; i++){
            for ( j = 0; j < nodos_elem; j++ ){
                //a = elems[ m ].conec[ i ] - 1;
                //b = elems[ m ].conec[ j ] - 1;
                a = elems[ m ].conec[ i ]; // renglón
                b = elems[ m ].conec[ j ]; // columna
                for ( k = 0; k < tabla_ady[ a ].cols; k++){
                    if ( tabla_ady[ a ].ind[ k ] == b ){
                        tabla_ady[ a ].val[ k ] += K_elem[i][j];
                        //printf ( "\nK [ %d ][ %d ] = %lf", i, j, K_elem[ i ][ j ]);
                        //printf ( "\ntabla_ady[ %d ].val[ %d ] = %lf\n", a, b, tabla_ady[ a ].val[ b ] );
                    }
                }
            }
            tabla_ady[ a ].f += f_elem[ i ];
        }
    }

    /*printf ( "\n\nMATRIZ DE RIGIDEZ COMPRIMIDA \n\n" );
    for ( i = 0; i < nodos; i++){
        for ( j = 0; j < tabla_ady[ i ].cols; j++ ){
            printf ( "%lf\t", tabla_ady[ i ].val[ j ] );
            if ( j == tabla_ady[ i ].cols -1 ) printf( "\n");
        }
    }*/

    /*printf ( "\n\nVECTOR DE FLUJOS \n" );
    for ( i = 0; i < nodos; i++ ){
        printf ( "\n%lf", tabla_ady[ i ].f );
    }*/

    libera_matriz ( K_elem, nodos_elem );
    libera_matriz ( coord_aux, nodos_elem );
    libera_vector ( f_elem );
    libera_vector ( prop_aux );

    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        libera_matriz ( pg, pg_elem );
        libera_matriz ( wpg, pg_elem );
        libera_matriz ( B, dimension );
        libera_matriz ( BT, nodos_elem );
        libera_matriz ( Baux, dimension );
        libera_matriz ( K_elemaux, nodos_elem );
        libera_matriz ( derin, dimension );
        libera_matriz ( difusion, dimension );
        libera_vector ( vect_aux );
        libera_vector ( vect_aux2 );
        libera_matriz ( jacobiano, dimension );
    }

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void puntos_gauss ( double **pg, double **wpg, int tipo_elem )
{
	/** Descripción **/
	// Devuelve los puntos de Gauss y los pesos que deberán ser evaluados según
	// el tipo de elemento que se emplee

	/** Variables de entrada **/
	// pg[][]		Apuntador a la matriz pg
	// wpg[][]		Apuntador a la matriz wpg
	// tipo_elem	Tipo de elemento que depende del número de nodos y el número de
	//			puntos de Gauss que se emplearán

	/** Variables de salida **/
	// pg[][]		Contine los puntos de Gauss
	// wpg[][]		Contiene los pesos

    if ( dimension == 1 ){
        printf( "\nEntra a 1D para determinar P-G.");
        exit( 1 );
    }
    else if ( dimension == 2 ){
        switch ( tipo_elem ){
            case 1:
            {
                // falta comprobar
                pg[0][0] = 1.0 / 3.0;
                pg[0][1] = 1.0 / 3.0;
                wpg[0][0] = 1.0 / 2.0;
                wpg[0][1] = 1.0;

                break;
            }
            case 2:
            {
                pg[0][0] = - sqrt ( 3.0 ) / 3.0;
                pg[0][1] = - sqrt ( 3.0 ) / 3.0;
                wpg[0][0] = 1.0;
                wpg[0][1] = 1.0;

                pg[1][0] =   sqrt ( 3.0 ) / 3.0;
                pg[1][1] = - sqrt ( 3.0 ) / 3.0;
                wpg[1][0] = 1.0;
                wpg[1][1] = 1.0;

                pg[2][0] =   sqrt ( 3.0 ) / 3.0;
                pg[2][1] =   sqrt ( 3.0 ) / 3.0;
                wpg[2][0] = 1.0;
                wpg[2][1] = 1.0;

                pg[3][0] = - sqrt ( 3.0 ) / 3.0;
                pg[3][1] =   sqrt ( 3.0 ) / 3.0;
                wpg[3][0] = 1.0;
                wpg[3][1] = 1.0;

                break;
            }
            default:
            {
                puts ( "Error: no se pudieron determinar los puntos de Gauss. Caso no incluido. " );
                exit ( 1 );
            }
        }
    }
    else{
        switch ( tipo_elem ){
            case 3:
            {

                /*pg[0][0] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                pg[0][1] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                pg[0][2] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                wpg[0][0] = 1.0/24.0;
                wpg[0][1] = 1.0;
                wpg[0][2] = 1.0;

                pg[1][0] = ( 5.0 + 3 * sqrt( 5.0 ) ) / 20.0;
                pg[1][1] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                pg[1][2] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                wpg[1][0] = 1.0/24.0;
                wpg[1][1] = 1.0;
                wpg[1][2] = 1.0;

                pg[2][0] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                pg[2][1] = ( 5.0 + 3.0 * sqrt( 5.0 ) ) / 20.0;
                pg[2][2] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                wpg[2][0] = 1.0/24.0;
                wpg[2][1] = 1.0;
                wpg[2][2] = 1.0;

                pg[3][0] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                pg[3][1] = ( 5.0 - sqrt( 5.0 ) ) / 20.0;
                pg[3][2] = ( 5.0 + 3.0 * sqrt( 5.0 ) ) / 20.0;
                wpg[3][0] = 1.0/24.0;
                wpg[3][1] = 1.0;
                wpg[3][2] = 1.0;*/

                pg[0][0] = 0.25;
                pg[0][1] = 0.25;
                pg[0][2] = 0.25;
                wpg[0][0] = 1.0 / 6.0;
                wpg[0][1] = 1.0;
                wpg[0][2] = 1.0;

                break;
            }
            case 4:
            {

                pg[0][0] = - sqrt ( 3.0 ) / 3.0;
                pg[0][1] = - sqrt ( 3.0 ) / 3.0;
                pg[0][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[0][0] = 1.0;
                wpg[0][1] = 1.0;
                wpg[0][2] = 1.0;

                pg[1][0] =   sqrt ( 3.0 ) / 3.0;
                pg[1][1] = - sqrt ( 3.0 ) / 3.0;
                pg[1][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[1][0] = 1.0;
                wpg[1][1] = 1.0;
                wpg[1][2] = 1.0;

                pg[2][0] =   sqrt ( 3.0 ) / 3.0;
                pg[2][1] =   sqrt ( 3.0 ) / 3.0;
                pg[2][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[2][0] = 1.0;
                wpg[2][1] = 1.0;
                wpg[2][2] = 1.0;

                pg[3][0] = - sqrt ( 3.0 ) / 3.0;
                pg[3][1] =   sqrt ( 3.0 ) / 3.0;
                pg[3][2] = - sqrt ( 3.0 ) / 3.0;
                wpg[3][0] = 1.0;
                wpg[3][1] = 1.0;
                wpg[3][2] = 1.0;

                pg[4][0] = - sqrt ( 3.0 ) / 3.0;
                pg[4][1] = - sqrt ( 3.0 ) / 3.0;
                pg[4][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[4][0] = 1.0;
                wpg[4][1] = 1.0;
                wpg[4][2] = 1.0;

                pg[5][0] =   sqrt ( 3.0 ) / 3.0;
                pg[5][1] = - sqrt ( 3.0 ) / 3.0;
                pg[5][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[5][0] = 1.0;
                wpg[5][1] = 1.0;
                wpg[5][2] = 1.0;

                pg[6][0] =   sqrt ( 3.0 ) / 3.0;
                pg[6][1] =   sqrt ( 3.0 ) / 3.0;
                pg[6][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[6][0] = 1.0;
                wpg[6][1] = 1.0;
                wpg[6][2] = 1.0;

                pg[7][0] = - sqrt ( 3.0 ) / 3.0;
                pg[7][1] =   sqrt ( 3.0 ) / 3.0;
                pg[7][2] =   sqrt ( 3.0 ) / 3.0;
                wpg[7][0] = 1.0;
                wpg[7][1] = 1.0;
                wpg[7][2] = 1.0;

                break;
            }
            default:
            {
                puts ( "Error: no se pudieron determinar los puntos de Gauss. Caso no incluido. " );
                exit ( 1 );
            }
        }
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


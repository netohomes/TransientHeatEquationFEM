//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void carac_elem ( int *n_prop );
void ordena_nodos ( struct elem *elems, int n_elem, int nodos_elem );
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void lectura ( char *r )
{
    /** Descripción **/
    // Lee los datos proporcionados por GID

    /** Variables de entrada **/
    // tipo_elem            Tipo de elemento a emplear definido desde main
    // nodos_elem           Núm. de nodos del elemento que se empleó para mallar
    // dimension            Dimensión

    /** Variables de salida **/
    // nodos            Núm. de nodos de la malla
    // nodos_elem       Núm. de nodos del elemento que se empleó para mallar
    // elems[].conec    Contiene los nodos de los elementos
    // coord[][]        Guarda las coordenadas de todos los nodos
    // n_prop           Núm. de propiedades a guardar dependiendo de la dimension en la que
    //                  se trabaje



    int i, j, aux1, aux2, n_prop = 0;
    char *ruta;
    FILE *fp;

    // Inicializa el núm. de las condicones de contorno
    nQ = 0;
    TP = 0;
    FP = 0;
    TL = 0;
    FL = 0;
    TS = 0;
    FS = 0;

    ruta = ( char *) calloc ( 256, sizeof ( char ) );

    strcpy ( ruta, r );
    strcat ( ruta, ".dat" );
    //printf ( "\n%s", ruta );

    fp = fopen ( ruta , "r" );
	if (fp == NULL){
		printf ( "Archivo de datos no encontrado\n");
		exit (1 );
	}

    // Dimensión
    fscanf ( fp, "NDIME %d", &dimension );

    // SOLO SE CUMPLE PARA CUANDO LOS ELEMENTOS QUE SE EMPLEAN SON LOS MISMOS
    // 2 -> Triángulo; 3-> Cuadrilátero
    fscanf ( fp, "\nELEMTYPE %d", &tipo_elem );
    tipo_elem -= 1;

     // Ajusta la dimensión para 1D
    if ( tipo_elem == 0 ) dimension = 1;

    // Determina los nodos del elemento y los P.G. de acuerdo al tipo de
    // elemento. También el número de propiedades a guardar.
    carac_elem( &n_prop );

    // Número de elementos
    fscanf ( fp, "\nELEMENTS %d", &n_elem );

	elems = mem_vector_elem ( n_elem, "de conectividad y prop. de elementos");

	for ( i = 0; i < n_elem; i++ )
	{
        elems[i].conec = mem_vector_int ( nodos_elem, "de conectividad." );
        elems[i].prop = mem_vector ( n_prop, "de propiedades." );

	}

    // Conectividades y propiedades de los elementos
    if ( dimension == 1 ){
        fscanf ( fp, "\n\nELEMENT  CONEC  K  Q\n" );
        // Nodos de los elementos (conectividad)
        for ( i = 0; i < n_elem; i++ ){
            fscanf ( fp, "%d", & elems[i].num );
            elems[i].num -= 1;
            for ( j = 0; j < nodos_elem; j++ ){
                fscanf ( fp, "%d", & elems[i].conec[j] );
                elems[i].conec[j] -= 1;
            }

            // Propiedades de los elementos. Se incluye Q por si no se quiere
            // determinarlo para cada elemento. De todas formas se puede cambiar
            // los elementos que sean distintos o todos.
            fscanf ( fp, "%lf", & elems[i].prop[ 0 ] );
            fscanf ( fp, "%lf", & elems[i].prop[ 1 ] );
        }
    }
    else if ( dimension == 2 ){
        fscanf ( fp, "\n\nELEMENT  CONEC  Kx  Ky  Q\n" );
        // Nodos de los elementos (conectividad)
        for ( i = 0; i < n_elem; i++ ){
            fscanf ( fp, "%d", & elems[i].num );
            elems[i].num -= 1;
            for ( j = 0; j < nodos_elem; j++ ){
                fscanf ( fp, "%d", & elems[i].conec[j] );
                elems[i].conec[j] -= 1;
            }

            // Propiedades de los elementos. Se incluye Q por si no se quiere
            // determinarlo para cada elemento. De todas formas se puede cambiar
            // los elementos que sean distintos o todos.
            fscanf ( fp, "%lf", & elems[i].prop[ 0 ] );
            fscanf ( fp, "%lf", & elems[i].prop[ 1 ] );
            fscanf ( fp, "%lf", & elems[i].prop[ 2 ] );
        }
    }
    else{
        fscanf ( fp, "\n\nELEMENT  CONEC  Kx  Ky  Kz  Q\n" );
        // Nodos de los elementos (conectividad)
        for ( i = 0; i < n_elem; i++ ){
            fscanf ( fp, "%d", & elems[i].num );
            elems[i].num -= 1;
            for ( j = 0; j < nodos_elem; j++ ){
                fscanf ( fp, "%d", & elems[i].conec[j] );
                elems[i].conec[j] -= 1;
            }

            // Propiedades de los elementos. Se incluye Q por si no se quiere
            // determinarlo para cada elemento. De todas formas se puede cambiar
            // los elementos que sean distintos o todos.
            fscanf ( fp, "%lf", & elems[i].prop[ 0 ] );
            fscanf ( fp, "%lf", & elems[i].prop[ 1 ] );
            fscanf ( fp, "%lf", & elems[i].prop[ 2 ] );
            fscanf ( fp, "%lf", & elems[i].prop[ 3 ] );
        }
    }

    // Número de nodos
	fscanf ( fp, "\n\nPOINTS %d", &nodos);

    if ( dimension == 1 ){
        coord = mem_matriz ( nodos, 2, "de coordenadas");
    }
    else{
         coord = mem_matriz ( nodos, dimension, "de coordenadas");
    }

	fscanf ( fp, "\nNODE COORD\n" );

	if ( dimension == 1 ){
        for ( i = 0; i < nodos; i++ ){
            for ( j = 0; j < 2 + 1; j++ ){
                if ( j == 0 ){
                    fscanf ( fp, "%d", &aux1 );
                }
                else{
                    // Matriz de coordenadas
                    fscanf ( fp, "%lf", &coord[ i ][ j -1 ] );
                }
            }
        }
	}
    else{
        for ( i = 0; i < nodos; i++ ){
            for ( j = 0; j < dimension + 1; j++ ){
                if ( j == 0 ){
                    fscanf ( fp, "%d", &aux1 );
                }
                else{
                    // Matriz de coordenadas
                    fscanf ( fp, "%lf", &coord[ i ][ j -1 ] );
                }
            }
        }
    }

    // Nuevo Q de los elementos
    fscanf ( fp, "\nSOURCE %d", &nQ );
	fscanf ( fp, "\n\nELEMENT Q\n" );
    for ( i = 0; i < nQ; i++ ){
        fscanf( fp, "%d", &aux1 );
        fscanf( fp, " %d", &aux2 );
        elems[ aux1 - 1 ].prop[ 2 ] = aux2;
    }

    // Temperaturas sobre puntos
    fscanf ( fp, "\nTEMPERATURE_POINTS %d", &TP );
    fscanf( fp, "\nNODE T\n" );
    t_points_ind = mem_vector_int ( TP, "condicion t_points_ind" );
    t_points_val = mem_vector ( TP, "condicion t_points_val" );
    for ( i = 0; i < TP; i++ ){
        fscanf( fp, "%d",  &t_points_ind[ i ] );
        fscanf( fp, "%lf", &t_points_val[ i ] );
        t_points_ind [ i ] -= 1;
    }

    if ( dimension == 1 ){
        //Flujo sobre puntos
        fscanf ( fp, "\nFLUX_POINTS %d", &FP );
        fscanf( fp, "\nNODE F\n" );
        f_points_ind = mem_vector_int ( FP, "condicion t_points_ind" );
        f_points_val = mem_vector ( FP, "condicion t_points_val" );
        for ( i = 0; i < FP; i++ ){
            fscanf( fp, "%d",  &f_points_ind[ i ] );
            fscanf( fp, "%lf", &f_points_val[ i ] );
            f_points_ind [ i ] -= 1;
        }
    }

    // Temperaturas sobre líneas
    fscanf ( fp, "\nTEMPERATURE_LINE %d", &TL );
    fscanf( fp, "\nNODE T\n" );
    t_line_ind = mem_vector_int ( TL, "condicion t_points_ind" );
    t_line_val = mem_vector ( TL, "condicion t_points_val" );
    for ( i = 0; i < TL; i++ ){
        fscanf( fp, "%d",  &t_line_ind[ i ] );
        fscanf( fp, "%lf", &t_line_val[ i ] );
        t_line_ind [ i ] -= 1;
    }

    if ( dimension == 2 ){
        // Flujos sobre líneas
        fscanf ( fp, "\nFLUX_LINE %d", &FL );
        fscanf ( fp, "\nNODE F(X) F(Y)\n" );
        f_line_ind = mem_vector_int ( FL, "condicion t_points_ind" );
        f_line_val = mem_matriz ( FL, dimension, "condicion t_points_val" );
        for ( i = 0; i < FL; i++ ){
            fscanf( fp, "%d",  &f_line_ind[ i ] );
            fscanf( fp, "%lf", &f_line_val[ i ][ 0 ] );
            fscanf( fp, "%lf", &f_line_val[ i ][ 1 ] );
            f_line_ind [ i ] -= 1;
        }
    }

    if ( dimension == 3 ){
        // Temperaturas sobre superficies
        fscanf ( fp, "\nTEMPERATURE_SURFACE %d", &TS );
        fscanf( fp, "\nNODE T\n" );
        if ( tipo_elem == 3 ){
            t_surf_ind = mem_matriz_int ( TS, 3, "condicion t_surf_ind" );
        }
        else{
            t_surf_ind = mem_matriz_int ( TS, 4, "condicion t_surf_ind" );
        }
        t_surf_val = mem_vector ( TS, "condicion t_surf_val" );
        for ( i = 0; i < TS; i++ ){
            if ( tipo_elem == 3 ){
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 0 ] );
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 1 ] );
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 2 ] );
                fscanf( fp, "%lf", &t_surf_val[ i ] );
                t_surf_ind[ i ][ 0 ] -= 1;
                t_surf_ind[ i ][ 1 ] -= 1;
                t_surf_ind[ i ][ 2 ] -= 1;
            }
            else{
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 0 ] );
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 1 ] );
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 2 ] );
                fscanf( fp, "%d",  &t_surf_ind[ i ][ 3 ] );
                fscanf( fp, "%lf", &t_surf_val[ i ] );
                t_surf_ind[ i ][ 0 ] -= 1;
                t_surf_ind[ i ][ 1 ] -= 1;
                t_surf_ind[ i ][ 2 ] -= 1;
                t_surf_ind[ i ][ 3 ] -= 1;
            }
        }

        // Flujos sobre superficies
        fscanf ( fp, "\nFLUX_SURFACE %d", &FS );
        fscanf ( fp, "\nNODE F(X) F(Y) F(Z)\n" );
        if ( tipo_elem == 3 ){
            f_surf_ind = mem_matriz_int ( FS, 3, "condicion f_surf_ind" );
        }
        else{
            f_surf_ind = mem_matriz_int ( FS, 4, "condicion f_surf_ind" );
        }
        f_surf_val = mem_matriz ( FS, dimension, "condicion f_surf_val" );
        for ( i = 0; i < FS; i++ ){
            if ( tipo_elem == 3 ){
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 0 ] );
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 1 ] );
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 2 ] );
                f_surf_ind[ i ][ 0 ] -= 1;
                f_surf_ind[ i ][ 1 ] -= 1;
                f_surf_ind[ i ][ 2 ] -= 1;
            }
            else{
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 0 ] );
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 1 ] );
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 2 ] );
                fscanf( fp, "%d",  &f_surf_ind[ i ][ 3 ] );
                f_surf_ind[ i ][ 0 ] -= 1;
                f_surf_ind[ i ][ 1 ] -= 1;
                f_surf_ind[ i ][ 2 ] -= 1;
                f_surf_ind[ i ][ 3 ] -= 1;
            }
            fscanf( fp, "%lf", &f_surf_val[ i ][ 0 ] );
            fscanf( fp, "%lf", &f_surf_val[ i ][ 1 ] );
            fscanf( fp, "%lf", &f_surf_val[ i ][ 2 ] );
        }

        // Temperaturas iniciales sobre superficies
        fscanf ( fp, "\nT_INITIAL %d", &TI );
        fscanf ( fp, "\nNODE T\n" );
        if ( tipo_elem == 3 ){
            ti_surf_ind = mem_matriz_int ( TI, 3, "condicion f_surf_ind" );
        }
        else{
            ti_surf_ind = mem_matriz_int ( TI, 4, "condicion f_surf_ind" );
        }
        ti_surf_val = mem_vector ( TI, "condicion t_surf_val" );
        for ( i = 0; i < TI; i++ ){
            if ( tipo_elem == 3 ){
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 0 ] );
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 1 ] );
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 2 ] );
                ti_surf_ind[ i ][ 0 ] -= 1;
                ti_surf_ind[ i ][ 1 ] -= 1;
                ti_surf_ind[ i ][ 2 ] -= 1;
            }
            else{
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 0 ] );
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 1 ] );
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 2 ] );
                fscanf( fp, "%d",  &ti_surf_ind[ i ][ 3 ] );
                ti_surf_ind[ i ][ 0 ] -= 1;
                ti_surf_ind[ i ][ 1 ] -= 1;
                ti_surf_ind[ i ][ 2 ] -= 1;
                ti_surf_ind[ i ][ 3 ] -= 1;
            }
            fscanf( fp, "%lf", &ti_surf_val[ i ] );
        }

    }

    // Suma las temperaturas impuestas

    //T0[  ]

    // Iteraciones y tolerancia
    fscanf ( fp, "\nSIMULATION" );
    fscanf ( fp, "\nITERATIONS %d", &iter );
    fscanf ( fp, "\nERROR %lf", &tol );
    fscanf ( fp, "\nTRANSIENT %d", &transient );
    fscanf ( fp, "\nRHO %lf", &rho );
    fscanf ( fp, "\nALFA %lf", &alfa );
    fscanf ( fp, "\nSTEPS %d", &pasos );
    fscanf ( fp, "\nTIME_STEP %lf", &delta_t );

    //printf ( "\niter = %d", iter );
    //printf ( "\ntol = %lf", tol );

    // Acaba lectura de datos
	fclose ( fp );

    /*printf ( "\nLECTURA DE DATOS\n");
	printf ( "\nn_elem = %d", n_elem);
	printf ( "\nnodos = %d", nodos);
	printf ( "\n" );
    printf ( "\nConectividades y propiedades");
	for ( i = 0; i < n_elem; i++ ){
	    printf ( "\n%d", elems[i].num);
	    for ( j = 0; j < nodos_elem; j++){
            printf ( "\t%d", elems[i].conec[j] );
	    }
	    if ( dimension == 1 ){
            printf( "\nEntra a 1D, para imprimir en pantalla propiedades" );
            exit( 1 );
        }
	    if ( dimension == 2 ){
            printf ( "\t%lf", elems[i].prop[ 0 ] );
            printf ( "\t%lf", elems[i].prop[ 1 ] );
            printf ( "\t%lf", elems[i].prop[ 2 ] );
	    }
	    else{
            printf ( "\t%lf", elems[i].prop[ 0 ] );
            printf ( "\t%lf", elems[i].prop[ 1 ] );
            printf ( "\t%lf", elems[i].prop[ 2 ] );
            printf ( "\t%lf", elems[i].prop[ 3 ] );
	    }
	}
	printf ( "\n" );
	printf ( "\nCoordenadas de nodos\n");
	for ( i = 0; i < nodos; i++ ){
	    for ( j = 0; j < dimension; j++ )
	    {
	        printf ( "%lf\t", coord[i][j]);
	        if ( j == dimension -1 ) printf ( "\n" );
	    }
	}*/

    // Ordena nodos de los elementos
    //if ( dimension == 3 ) ordena_nodos( elems, n_elem, nodos_elem );

    // Libera nombre de archivo
    free ( ruta );

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void carac_elem ( int *n_prop )
{
    /** Descripción **/
    // Contiene los tipos de elemento que se pueden emplear


    if ( dimension == 1 ){
        *n_prop = 2;
        nodos_elem = 2;
        pg_elem = 0;
    }
    else if ( dimension == 2 ){
        *n_prop = 3;
        switch ( tipo_elem ){
            case 1:
            {
                nodos_elem = 3;
                pg_elem = 1;
                break;
            }
            case 2:
            {
                nodos_elem = 4;
                pg_elem = 4;
                break;
            }
            default:
            {
                printf( "\nError: no se pudo determinar nodos_elem y pg_elem. Caso no incluido (2D)." );
                exit ( 1 );
            }
        }
    }
    else{
         *n_prop = 4;
         switch ( tipo_elem ){
            case 3:
            {
                nodos_elem = 4;
                pg_elem = 1;
                break;
            }
            case 4:
            {
                nodos_elem = 8;
                pg_elem = 8;
                break;
            }
            default:
            {
                printf( "\nError: no se pudo determinar nodos_elem y pg_elem. Caso no incluido (3D)." );
                exit ( 1 );
            }
        }
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
void ordena_nodos ( struct elem *elems, int n_elem, int nodos_elem )
{
    /** Descripción **/
    //Ordena los nodos de los elementos de mayor a menor

    int m, i, j;
    double a, b, c;

    for ( m = 0; m < n_elem; m++ ){
        for ( i = 0; i < nodos_elem; i++ ){
            a = elems[ m ].conec[ i ];
            for ( j = i; j < nodos_elem; j++ ){
                b = elems[ m ].conec[ j ];
                if ( b < a ){
                    c = a;
                    a = b;
                    b = c;
                    elems[ m ].conec[ i ] = a;
                    elems[ m ].conec[ j ] = b;
                }
            }
        }
    }
}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

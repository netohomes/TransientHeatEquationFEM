//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

            /*** SOLUCI�N A LA ECUACI�N DE CALOR TRASCENDENTE POR ELEMENTO FINITO ***/

                            // Para 3D
    // Para elementos: tetraedros y hexaedros.
    // Dim. No.      Tipo de elemento    P.G.    C.C v�lidas
    //                                      Temperatura Flujos
    // 1    0        Lineales                    P,L         P
    // 2    1        Tri�ngulares        1       P,L         L
    // 2    2        Cuadril�teros       4       P,L         L
    // 3    3        Tetraedro           4       P,L,S       S
    // 3    4        Hexaedro            8       P,L,S       S
    // Nota: P -> Sobre puntos, L -> Sobre l�neas, S -> Sobre superficies
    // Elementos con el m�nimo de P.G. necesarios para la integraci�n.

//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
#include "funciones.h"
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
int main( int argc, const char *argv[] )
{
    //----------------------------------------------------------------------------------------------

    /** Descripci�n **/
    // Funci�n principal del programa

    /** Variables de entrada **/
    // argc     N�mero de argumentos para el programa
    // argv[]   Direcci�n al archivo con el que se est� trabajando


    /** Variables locales **/
    // T0       Vector temperaturas en t-1
    // T1       Vector temperaturas en t
    // q_pg     Flujos en P.G.
    // q_pg2    Flujos en P.G. calculado a partir de flujos en nodos
    // q_nodos  Flujos en nodos (suavizado)

    //----------------------------------------------------------------------------------------------

	double *T, **q_pg, **q_pg2, **q_nodos, **error_PG, *det;

	//int i;
	//int j;

    //----------------------------------------------------------------------------------------------
    // Lectura de datos
	lectura ( (char *)argv[ 1 ] );

	//----------------------------------------------------------------------------------------------

    // Determina la tabla de adyacencias
    ady ( elems, n_elem, nodos, nodos_elem );

    //----------------------------------------------------------------------------------------------

    // Calcula las matriz de rigidez global en el formato comprimido.
    // Tambi�n calcula el vector de flujos debido a Q.
    det = mem_vector ( n_elem, "det. de todos los elementos." );
	rigid ( tabla_ady,tipo_elem, nodos_elem, pg_elem, dimension, det );

	//----------------------------------------------------------------------------------------------

    // Impone las condiciones de contorno
	cond_cont ( tabla_ady, TP, TL, FL, TS, FS, FP, t_points_ind, t_points_val, f_points_ind,
                f_points_val, t_line_ind, t_line_val, f_line_ind, f_line_val, t_surf_ind,
                t_surf_val, f_surf_ind, f_surf_val, nodos, coord );

    //----------------------------------------------------------------------------------------------

    // Memoria para vectores soluci�n del problema
    T = mem_vector ( nodos, "VECTOR DE SOLUCION DE TEMPERATURAS");
    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        q_pg = mem_matriz ( pg_elem * n_elem, dimension, "flujos en P.G." );
        q_pg2 = mem_matriz ( pg_elem * n_elem, dimension, "flujos en P.G. desde nodos" );
        error_PG = mem_matriz ( pg_elem * n_elem, dimension,  "error en PG." );
    }
    q_nodos = mem_matriz ( nodos, dimension, "flujos en P.G." );

    // Inicializa las temperaturas a 0.0 ( falta inicializar a la cc que el usuario indique )
    inicializa_vector ( T, nodos );
    if ( ( dimension == 2 ) || ( dimension == 3 ) ){
        inicializa_matriz ( q_pg, pg_elem * n_elem, dimension );
        inicializa_matriz ( q_pg2, pg_elem * n_elem, dimension );
        inicializa_matriz ( error_PG, pg_elem * n_elem, dimension );
    }
    inicializa_matriz ( q_nodos, nodos, dimension );

    //----------------------------------------------------------------------------------------------

    //Evalua si el caso es transitorio y resuleve para temperaturas

    if ( transient == 1 ){

        if ( ( dimension == 1 ) || ( dimension == 2 ) ){
            printf( " \n Caso transitorio aun no listo para esta dimension, solo 3D " );
            exit(1);
        }

        trascendencia ( elems, tabla_ady, T, alfa, rho, delta_t, pasos, n_elem, nodos,
                        tipo_elem, nodos_elem, (char *)argv[ 1 ], det, pg_elem, TI,
                        ti_surf_ind, ti_surf_val, dimension );

    }
    else {

        // Resuelve y devuelte el vetor de temperaturas T
        grad_conj ( tabla_ady, T, nodos, nodos, tol );

    }

    //----------------------------------------------------------------------------------------------

    // Suaviza flujos

    if ( dimension == 1 ){
        flujos_1D ( elems, T, q_nodos, nodos, n_elem );
    }
    else{
        // C�lcula y suaviza los flujos
        suavizado_flujos ( elems, T, dimension, tipo_elem, n_elem, pg_elem, nodos_elem, coord, nodos,
                           q_pg, q_nodos, q_pg2, error_PG );
    }

    //----------------------------------------------------------------------------------------------

    // Escribe resultados. Crea el archivo para el postroceso en GID

    if ( transient == 1 ){

        escribe_flujos ( (char *)argv[ 1 ], q_pg, q_nodos, q_pg2, error_PG, pg_elem, dimension,
                          nodos, pasos * delta_t );

    }
    else{

        escribe_result ( (char *)argv[ 1 ], pg_elem, tipo_elem, nodos, T, q_pg, q_pg2, q_nodos,
                         error_PG );
    }

    //----------------------------------------------------------------------------------------------

    // Libera memoria para variables globales

    // Arreglos de main
    libera_mem1 ( T, q_pg, q_pg2, q_nodos, error_PG, dimension, det );

    // Variables globales
    libera_mem2 ( );

    //----------------------------------------------------------------------------------------------

	return 0;

}
//{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

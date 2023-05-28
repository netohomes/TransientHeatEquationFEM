# VARIABLES DE COMPILACIÓN
# Compilador
CC = gcc
# Opciones de compilación
FLAGS = -Wall -g
# Nombre del programa
EXE = output
# Source files
SOURCE_FILES = Heat_equation.c cond_cont.c ensambla.c flujos.c lectura.c memoria.c oper_matrices.c rigid_elem.c solver.c escribe_result.c trascendencia.c

all:
	@echo Compilando...
	$(CC) $(FLAGS) $(SOURCE_FILES) -o $(EXE)
clean:
	@echo Borrando ejecutable...
	rm $(EXE)

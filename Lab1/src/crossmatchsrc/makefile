CC   = gcc		# the C compiler
COPT = -O3 -Wall 	# options for the C compiler
LIBS = -lm		# libraries

OBJS = tri_gen.o triangles.o bright_end.o xy_lin.o indexx.o refine.o main.o

FUNCS = bright_end.o main.o triangles.o tri_gen.o 

DEFS = main.o triangles.o

.c.o :
	${CC} ${COPT} -c $<


crossmatch: ${OBJS}
	${CC} ${OBJS} -o $@ ${LIBS} ${COPT}

${OBJS} : makefile

${FUNCS}: funcs.h

${DEFS} : defs.h

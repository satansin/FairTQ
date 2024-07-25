CC = g++
CFLAGS = -w -std=c++11
OBJS = cellTree.o collection.o Array.o QuadProg++.o filemem.o global.o main.o param.o point.o queryGen.o rnode.o rtree.o skyline.o tgs.o utk.o virtualRNode.o rentry.o liblpsolve55.a hypercube.cpp libglpk.a
INCLUDES = -I "./headers/" 
LIBS = 
FOLDER = ../bin/

main: ${OBJS} ${FOLDER}
	${CC} ${CFLAGS} ${INCLUDES} -o ../bin/ftq ${OBJS}  ${LIBS} -lpthread -llpsolve55 -ldl -g ;\
	

clean: 
	rm -f *.o
	rm -rf ${FOLDER}

.cpp.o:
	${CC} ${CFLAGS} ${INCLUDES} -c -g $<

.cc.o:
	${CC} ${CFLAGS} ${INCLUDES} -c -g $<


${FOLDER}:
	mkdir ${FOLDER}

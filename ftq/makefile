all: make_qhull make_ftq

make_qhull:
	make -C qhull/

make_ftq:
	make -C FairTQ-Exact/

clean:
	make cleanall -C qhull/
	make clean -C FairTQ-Exact/

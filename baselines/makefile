all: make_qhull make_greedy make_prop make_sr_adapt

make_qhull:
	make -C qhull/

make_greedy:
	make -C Greedy/src/

make_prop:
	make -C Prop/src/

make_sr_adapt:
	make -C SR-Adapt/src/

clean:
	make cleanall -C qhull/
	make clean -C Greedy/src/
	make clean -C Prop/src/
	make clean -C SR-Adapt/src/

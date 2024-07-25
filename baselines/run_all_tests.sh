#!/bin/bash

for f in test1 test2; do
	cd UTK/bin/
	cp ../../$f/* .
	./ftq
	./ftq appr1
	cp stat_* ../../$f/
	cp time_* ../../$f/
	find . ! -name 'ftq' -type f -exec rm -f {} +
	cd ../../
done

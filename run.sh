#!/bin/bash

NUMA_PATH=../numactl-2.0.11/build

gcc -O2 -std=c99 -lpthread -lnuma -L$NUMA_PATH/lib -I$NUMA_PATH/include benchmark.c \
  -DNUM_THREADS=4 -DMEM_OFF=0 -o benchmark

g++ -O2 -std=c++11 -pthread nonuma.cc -DNUM_THREADS=4 -o nonuma

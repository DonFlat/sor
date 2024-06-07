### Compile C
`mpicc -O3 -o sor_c sor_sendrecv.c -lm`
`mpirun -n 2 sor_c 10`
`mpicc -O3 src/sor_rma_c.c -o sor_rma_c -lm`

### Run program
`mpirun -n 2 ./target/debug/sor rma 1`
`mpirun -n 2 sor_c 10`
`prun -np 2 -1 OMPI_OPTS="--mca btl tcp,self --mca btl_tcp_if_include ib0" -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor 10`
`prun -np 2 -1 OMPI_OPTS="--mca btl tcp,self --mca btl_tcp_if_include ib0" -script $PRUN_ETC/prun-openmpi `pwd`/./sor_rma_c 16384 1`
`prun -np 2 -1 -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor  10`
`prun -np 2 -1 -script $PRUN_ETC/prun-openmpi `pwd`/./sor_rma_c 16384 2`
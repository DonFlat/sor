### Compile C
`mpicc -O3 -o sor_c sor_sendrecv.c -lm`
`mpirun -n 2 sor_c 10`

### Run program
`mpirun -n 2 ./target/debug/pingpong 10`
`prun -np 2 -1 OMPI_OPTS="--mca btl tcp,self --mca btl_tcp_if_include ib0" -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor 10`


```
// TODO: revisit RMA paper, why small message is quicker?
// TODO: try ping-pong test
// TODO: SOR with different numbers of iterations, problem size. don't need to get correct result
// Some found for experiment with: 8 nodes, 10 iters,
// matrix size: 50.000, RMA: 8460.146122999999 ms, sendrev too long to get result
// Even C sor -O3 on 50K took too long to get result, but 30k took 3.166942 s, RMA is 3562.165272 ms
// matrix size: 40.000, RMA: 5720.041005 ms, sendrev also too long to get result, C sor is same
//
// TODO: seq C vs seq Rust

// TODO: library, refactor test, revisit RMA paper, know synchronization
// TODO: see CPU cycle
// start from simple case regardless of displacement unit or so.
```
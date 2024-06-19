### Compile C
`mpicc -O3 -o sor_c sor_sendrecv.c -lm`
`mpirun -n 2 sor_c 10 2`
`mpicc -O3 src/sor_rma_c.c -o sor_rma_c -lm`
`mpicc -O3 src/sor_rma_lock.c -o c_lock -lm`

### Run program
`mpirun -n 3 ./target/debug/sor rma 6 3`
`mpirun -n 2 sor_c 10 2`
`mpirun -n 2 c_lock 10 2`
`prun -np 2 -1 OMPI_OPTS="--mca btl tcp,self --mca btl_tcp_if_include ib0" -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor 6 3`
`prun -np 2 -1 OMPI_OPTS="--mca btl tcp,self --mca btl_tcp_if_include ib0" -script $PRUN_ETC/prun-openmpi `pwd`/./sor_rma_c 16384 1`
`prun -np 2 -1 -script $PRUN_ETC/prun-openmpi `pwd`/./target/release/sor 6 3`
`prun -np 2 -1 -script $PRUN_ETC/prun-openmpi `pwd`/./sor_rma_c 16384 2`
`prun -np 2 -1 -script $PRUN_ETC/prun-openmpi `pwd`/./c_lock 16384 2`

### The prompt describing data tables
I have five files and their label used for legends:
sendrecv_data_node.csv  -> sendrecv
huge_data_node.csv -> big window
raw_data_node.csv -> unsafe
sor_c_rma_node.csv -> c_rma
split_data_node.csv -> split
lock_data_node.csv -> lock

Each file has 14 columns, first column is the matrix size, second column is node number.
The rest columns are program running times.

Help me write Python program to draw line charts. x-axis is growing node number, y-axis is the mean running time, it should be calculated
by taking out the max and min value from 12 repetitions. Use the first row as the base, which has 2 nodes. Draw speed up graph by increasing nodes.

Different lines should in different colour, marker, style to make them distinguishable.
Don't show, save graph
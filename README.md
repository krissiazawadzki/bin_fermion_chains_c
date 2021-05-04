# bin_spin_chains_c
Package provides modules to exactly diagonalize small finite fermionic chains in the reduced Hilbert space (Q,S,Sz=S). 

- examples
reduced_hilbert_QS: generates the rotation matrix between basis (Q,Sz) -> (Q,S) and saves all the possible configurations of the chain with size nsites with charge = Q = Nf (filling) and spin z component = Sz
    usage:
    ./reduced_hilbert_QS nsites Q Sz path_to_folder_to_save_results

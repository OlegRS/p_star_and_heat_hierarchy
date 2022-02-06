# p-star model
In this project we analyse an exponential random graph model (ERGM) in which the sufficient statistics are the total numbers of p-stars in the network. (Note that the term p-star should not be confused with p*, which sometimes stands for ERGMs in general.) We also introduce its mean field analogue ("MF model") defined by the Hamiltonian that depends solely on the powers of the connectance. It turns out that the MF model can be solved using the PDE-based techniques (apart from more conventional methods) and these analytical results are in excellent agreement with Monte Carlo simulations implemented in this repo.

For better understanding see my PhD thesis (phd_thesis_senkevich.pdf) and https://github.com/OlegRS/ERGMs.

## Code
Images from the paper stored in ./figures directory were produced using the Python scripts which can be found in ./scripts directory with the same names. Most of these scripts rely on the data stored in ./data directory in the form of .csv files with descriptive names. The data needed for all the scripts apart from 1000_5000_nodes_low_temp_degree_distr.py and 1000_5000_nodes_low_temp_lcc_distr.py (for which the files were too large for github) are already in the right location within ./data directory, so running the scripts with Python 3 should reproduce images from the paper apart from ./figures/low_temp_LCC_and_degree_distr.png (for which the data is missing). 

C++ implementations of Monte Carlo simulations that produced all the data can be found in ./models/main_sim directory. To compile them cd to ./main_sim directory and run "./compile.sh". (To accelerate compilation comment out unneeded lines from compile.sh.) After compilation cd to ./models/bin (which will be created by the compile.sh), and execute the binary form there. The data in the form of .csv files with descriptive names will be saved in the right location in ./data directory (see the code to find out where exactly).

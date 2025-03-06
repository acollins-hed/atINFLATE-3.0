# atINFLATE-3.0
February-March 2025

Welcome to atINFLATE 3.0. AtINFLATE is the (a)aRS-(t)RNA (I)nteraction
(N)etwork (F)itness (LA)ndscape (T)opographer (E)xpress, a dynamic model of 
its predecessor, atINFLAT, published in Collins-Hed and Ardell 2019. 
Installation requires g++ because there is the g++ sepcific instruction 
__builtin_popcount() in Genotype.cpp and Evolver.cpp which finds the weight 
of a number in the coding theory sense (example: 5 in binary is 101 and 
would have a weight of 2 since there are 2 ones in its binary 
representation), and the Eigen 3.4.0 library is also required.

To install simply type

make EIGEN=/PATH/TO/EIGEN

for example, it might look like: make EIGEN=/home/user/Eigen/eigen-3.4.0

and it will compile all the code. Type

./atinflate_3 -h

or

./atinflate_3 --help

and something like a "man" file will echo to the screen. An example of
a simulation you might want to run is

./atinflate_3 --num-traj=10 --halting-fixation=1000 --halting-fitness=0.9 --uniform-amino-acids --phi=0.8 --Mu=0.1 --proofreading --seed=137 -0 -f=5e-7 --Site-Types=50,50,50,50,50 -n=12 -T=5 -A=5 -o=test

The above will run 3 trajectories (--num-traj=3) for 1000 fixations each
(--halting-fixation=1000) and it will only halt if it reaches 1000 fixations
or the fitness is > 0.9 (--halting-fitness=0.9). --uniform-amino-acids sets
the amino acid physicochemical values to be distributed uniformly along [0,1],
--phi=0.8 sets the missense tolerance to 0.8, --Mu=0.1 sets the per codon
mutation rate to 0.1, --proofreading enables kinetic proofreading, --seed=137
sets the random seed to 137, -0 sets the p parameter for the bernoulli
distribution for each bit in the interaction interface to 0, -f=5e-7 sets the
rate selection constant to 5x10^(-7), --Site-Types=50,50,50,50,50 sets all
5 site-type frequencies to 50, -n=12 sets the interaction interface bit
length to 12, -T=5 and -A=5 set the number of tRNAs and aaRSs to 5 each (it
also sets the number of site-types = 5 by default, unless a different number
of site-types is specified). The final parameter -o=test will cause all the
data and log file names to be prepended with "test", for example the
trajectory datafile will be called test_traj.dat and the log file will be
named test.log.

Adding this line to test commit verification.
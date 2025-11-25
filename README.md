# atINFLATE-3.0
February-March 2025

Welcome to atINFLATE 3.0. AtINFLATE is the (a)aRS-(t)RNA (I)nteraction
(N)etwork (F)itness (LA)ndscape (T)opographer (E)xpress, a dynamic model of 
its predecessor, atINFLAT, published in [Collins-Hed and Ardell 2019](https://doi.org/10.1016/j.tpb.2019.03.007). 

## Getting Started
Clone the atINFLATE-3.0 git repository and change directory into your local
clone to compile atinflate_3

## Pre-Installation of Dependencies
Compiling atinflate_3 requires a local copy of Eigen 3.4.0, the
[Eigen C++ template library for linear algebra](https://eigen.tuxfamily.org/dox/GettingStarted.html). `atinflate_3` uses the KroneckerProduct module from the
`unsupported` collection of modules of Eigen.

Compilation requires `g++` because there is the `g++` sepcific instruction 
`__builtin_popcount()` in Genotype.cpp and Evolver.cpp which finds the weight 
of a number in the coding theory sense. Example: 5 in binary is 101 and 
would have a weight of 2 since there are 2 ones in its binary 
representation.

## Compiling
Once you have confirmed you have all the dependencies, open a terminal and
navigate to the atINFLATE-3.0 directory that you cloned from GitHub and
type

`make EIGEN=/PATH/TO/EIGEN`

for example, it might look like: `make EIGEN=/home/user/Eigen/eigen-3.4.0`

and it will compile all the code. If you want to run the code through a debugger, use

`make EIGEN=PATH/TO/EIGEN DEBUG=TRUE`

and it will compile using `-g3` and `-O0` instead of with `-O2` and no `-g3`.

## Help with Running
Once compiled, type

`bin/atinflate_3 -h`

or

`bin/atinflate_3 --help`

and something like a "man" page will echo to the screen.

## Example of Use
An example of a simulation you might want to run is

`bin/atinflate_3 --num-traj=10 --halting-fixation=1000 --halting-fitness=0.9 --uniform-site-types \`
`--phi=0.8 --Mu=0.1 --proofreading --seed=137 -0 -f=5e-7 --Site-Type-freq=50,50,50,50,50 -n=12 -T=5 \`
`-A=5 --Nthread=8 -o=test`

The above will run 3 trajectories (`--num-traj=3`) for 1000 fixations each
(`--halting-fixation=1000`) and it will only halt if it reaches 1000 fixations
or the fitness is > 0.9 (`--halting-fitness=0.9`). `--uniform-site-types` sets
the site-type physicochemical values to be distributed uniformly along [0,1],
`--phi=0.8` sets the missense tolerance to 0.8, `--Mu=0.1` sets the per codon
mutation rate to 0.1, `--proofreading` enables kinetic proofreading,
`--seed=137` sets the random seed to 137, `-0` sets the p parameter for the
bernoulli distribution for each bit in the interaction interface to 0, `-f=5e-7`
sets the rate selection constant to $5\times 10^{-7}$,
`--Site-Type-freq=50,50,50,50,50` sets
all 5 site-type frequencies to 50, `-n=12` sets the interaction interface bit
length to 12, `-T=5` and `-A=5` set the number of tRNAs and aaRSs to 5 each (it
also sets the number of site-types = 5 by default, unless a different number
of site-types is specified). `--Nthread=8` sets the number of threads for
multithreading to 8 threads. The final parameter `-o=test` will cause all the
data and log file names to be prepended with "test", for example the
trajectory datafile will be called test_traj.dat and the log file will be
named test.log.

## Arguments
What follows is a description of all the arguments for atINFLATE 3.0

### `-0`
Initial genotype is all 0 genotype. `-0` sets the p parameter for the Bernoulli distribution used for each bit in the interaction interfaces equal to 0. Equivalent to `bin/atinflate_3 --bp=0`.

By default each bit follows a Bernoulli distribution.

### `-1`
Initial genotype is all 1 genotype. `-1` sets the p parameter for the Bernoulli distribution used for eaach bin in the interaction interfaces equal to 1. Equivalent to `bin/atinflate_3 --bp=1`.

By default each bit follows a Bernoulli distribution.

### `-A=x`, `--aaRS=x`
Sets the number of aaRSs to `x`. Obviously, `x` must be a positive integer that is at most the number of site-types. If `x` is less than the number of site-types, then only the first `x` site-type physicochemical values will have corresponding aaRSs. To avoid a bias in aaRS physicochemical values, like in the event the `--uniform-site-types` argument is also used, the site-type physicochemical values are randomly permuted first and then the first `x` are chosen to determine which sites have a corresponding amino acid and aaRS.

By default `x` is equal to the number of site-types. The default number of site-types is 4.

### `--bp=x`
Sets the p parameter for the Bernoulli distribution used for each bit in the interaction interfaces equal to `x`. That is, `x` is the probability of a 1 and 1-`x` is the probability of 0. Because `x` is a probability, it must be in [0,1].

By default `x` = 0.5.

### `--codon-space-x`
where `x` is in {2,4} representing the number of bases. This also switches from the default space. `--codon-space-2` can only be used with `-T=2`, `-T=4`, or `-T=8`. `--codon-space-4` can only be used with `-T=4`, `-T=16`, or `-T=64`.

By default the codon space is a ring space where each codon can only mutate to one of two neighbors on the ring.

### `-f=x`
Sets the rate selection parameter to `x`. `x` can be any positive real number where a number near 0 represents high rate selection and a larger number represents low rate selection.

By default `x` $= \frac{1}{44}$ for no proofreading and `x` $= \frac{1}{9680}$ for proofreading.

### `--halting-fitness=x`
Here, `x` is the fitness threshold at which the trajectory will stop and the simulation will either move onto the next trajectory or, if there are no trajectories left, cease the simulation. So, if `x` is 0.8, then once the fitness reaches 0.8 or higher, the trajectory is terminated and the simulation moves onto the next trajectory if there is one or ceases. `x` must be in (0,1], thus if the user wants each trajectory to halt at a certain fixation regardless of fitness, set `x` equal to 1, since 1 is unobtainable for any system.

By default `x` is 0.001.

### `--halting-fixation=x`
Here, `x` is the fixation threshold at which the trajectory will stop and the simulation will either move onto the next trajectory or, if there are no trajectories left, cease the simulation. So, if `x` is 1000, then once the trajectory reaches the ${\tt 1000}^{\text{th}}$ fixation, the trajectory is terminated and the simulation moves on the next trajectory if there is one or ceases. `x` must be a positive integer.

By default `x` is 1.

### `-i=file_name`, `--ifile=file_name`
Here, `file_name` is the name of the input data you want to continue running, make sure you use the name associated with the file_name_checkpoint.log file you want. Note that the parameter `-o=file_name` (also `--ofile=file_name`) will create several files named file_name_x.dat as well as file_name.log and file_name_checkpoint.log. To run `-i=file_name` (or `--ifile=file_name`) you only need to enter file_name NOT file_name_checkpoint.log. The latter will trigger an error. So, for example `-i=file_name` will work, but `-i=file_name_checkpoint.log` will not.

There is no default.

### `--kmax=x`
`x` is a positive number representing the maximum dissociation rate between tRNAs and aaRSs measured in ${\tt seconds}^{-1}$. `x` must be strictly larger than the value for `kmin`.

By default `x` is ${\tt 10}^4$.

### `--kmin=x`
`x` is a positive number representing the minimum dissociation rate between tRNAs and aaRSs measured in ${\tt seconds}^{-1}$. `x` must be strictly smaller than the value for `kmax`.

By default `x` is 220.

### `--mu=x`
`x` is the mutation rate of identifying features and epistatic modifiers. It must be in (0,1) and at least an order of magnitude smaller than the codon mutation rate, `--Mu`.

By default `x` is ${\tt 10}^{-6}$.

### `--Mu=x`
`x` is the codon mutation rate. This refers to the mutation rate of the codon quasispecies. `x` must be in (0,1) and be at least an order of magnitude larger than the mutation rate for identifying features and epistatic modifiers.

By default `x` is ${\tt 10}^{-4}$.

### `-n=x`
`x` is the interaction interface word length, that is, the number of bits that make up the interaction interfaces. `x` must be a positive integer.

By default `x` is 4.

### `-N=x`, `--popsize=x`
`x` is the population size which must be a positive integer. The larger `x` is the greater the chance of fixing on a genotype of higher fitness and the lower the chance of fixing on genotype of lower fitness.

By default `x` is 100.

### `--Nthread=x`
`x` is the number of threads to use for multithreading. Note that the optimal number will vary from system to system and the user may need to run time trials to determine the best value of `x`.

By default `x` is 1.

### `--no-double-mutants`
Double mutants are extremely unlikely under certain conditions, thus the user might desire to exclude the possibility of fixing on a mutant two mutations away from the current genotype. This can markedly speed up the simulations.

By default the simulation will check every mutant within two mutations of the current wildtype when determining the next fixed mutant.

### `--no-mask`
This argument will turn off masking. Masking is used for studying the impact and evolution of epistatic modifiers, but the user may wish to turn this off either because it is not relevant to their research or the user wishes to see the behavior of the network without epistatic modifiers.

By default masking is ON.

### `--no-rate`
This argument removes any selection on translation rate. Thus, the fitness will only depend on the accuracy of translation.

By default rate selection is ON.

### `--num-traj=x`
`x` is the number of trajectories the user wishes to run in a single simulation. `x` must be a positive integer.

By default `x` is 1.

### `-o=file_name`, `--ofile=file_name`
This names the output files. So if, for example, `file_name` is "run", then the output files will be named run.log, run_checkpoint.log, and run_x.dat where the x represents all the different datafiles.

By default `file_name` is "run".

### `--phi=x`
Here, `x` is the missense tolerance, in other words, how well a site can tolerate being translated as a non-target amino acid. `x` must be in (0,1). `x` near 1 has a high tolerance and thus a weaker selection on translation accuracy, `x` near 0 has a low tolerance and thus a stronger selection on translation accuracy.

By default `x` is 0.99.

### `--proofreading`
This argument applies kinetic proofreading. Kinetic proofreading is an error correction mechanism found in many biochemical reactions, including interactions between tRNAs and aaRSs.

By default kinetic proofreading is OFF.

### `--seed=x`
`x` is the random seed for the simulation.

By default the random seed is generated with a random device from std::random_device.

### `--trans-bias=x`
`x` is the transition bias, that is the bias in favor of transtion mutations over transversions. A transition mutation refers to mutations between bases of the same kind, that is from one purine to another or one pyrimidine to another, and transversions are mutations from purine to pyrimidine or visa versa. `x` must be at least 1, where `x` equal to 1 means there is no bias and `x`>1 implies a transition bias. This parameter can only be used with `--codon-space-4` as it has no meaning in either `--codon-space-2` or codon ring space.

By default `x` is 1.

### `-S=x`
`x` is the number of site-types. `x` must be a positive integer that is at least the number of aaRSs.

By default `x` is 4.

### `--Site-Type-freq=x,y,z,...`
Here, `x,y,z,...` refer to the
absolute frequency of each site-type. All of `x`, `y`, `z`, and `...`
must be positive integers.

By default all site-types have a frequency of 1.

### `--Site-Type-pvs=x,y,z,...`
Here, `x,y,z,...` refer to the physicochemical values (pvs) of each site-type. All of `x`, `y`, `z`, and `...` must be unique, nonnegative, and at most 1. When using this argument with `--Site-Type-freq=` the site-type frequencies are assigned in the same order as `--Site-Type-pvs=`. For example,

`bin/atinflate_3 --Site-Type-freq=1,2,3,4 --Site-Type-pvs=0.1,0.49,0.51,0.999`

would result in a simulation where there is a site-type with pv 0.1 and frequency 1, a site-type with pv 0.49 and frequency 2, a site-type with pv 0.51 and frequency 3, and a site-type with pv 0.999 and frequency 4. If the user wants fewer aaRSs than site-types, then the user must also specify the number of site-types. For example,

`bin/atinflate_3 --Site-Type-freq=1,2,3,4 --Site-Type-pvs=0.1,0.49,0.51,0.999 -S=4 -A=3`

In this case the system would have site-types of pv 0.1, 0.49, 0.51, and 0.999 with frequencies 1, 2, 3, and 4 respectively, but only aaRSs with amino acids for the first three site-types (0.1, 0.49, and 0.51) would be present in the system. There would be a site-type with no corresponding amino acid, i.e. 0.999. 

By default all site-types are randomly drawn from the interval [0,1].

### `-T=x`, `--tRNAs=x`
`x` is the number of tRNAs in the system. `x` must be a positive integer. `x` will also determine the number of codons.

By default `x` is 4.

### `--uniform-site-types`
This argument distributes the physicochemical values of the site-types uniformly across [0,1]. So if `-S=5`, then there will be a site-type with physicochemcial value $0$, one with $\frac{1}{4}$, one with $\frac{1}{2}$, one with $\frac{3}{4}$, and one with $1$ as a physicochemical value.

By default site-type and amino acid physicochemical values are randomly drawn from the interval [0,1].

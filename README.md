#Reweighing-TASS-1.2

# Brief Description

Temperature Accelarated Sliced Sampling (TASS) method combines the temperature accelerated molecular dynamics with umbrella sampling and
metadynamics to sample the collective variable space in an efficient manner. \
[Ref :\
Exploring high dimensional free energy landscapes: Temperature accelerated sliced sampling J. Chem. Phys. 146, 094108 (2017).\
[![DOI] https://doi.org/10.1063/1.4977704 ]
Awasthi, S, Nair, NN. Exploring high‐dimensional free energy landscapes of chemical reactions. WIREs Comput Mol Sci. 2019.\
[![DOI]  https://doi.org/10.1002/wcms.1398 ]
Pal, A., Pal, S., Verma, S., Shiga, M., Nair, N. N., Mean force based temperature accelerated sliced sampling: Efficient reconstruction of high dimensional free energy landscapes .\
[![https://doi.org/10.1002/jcc.26727] ]

This Modular Fortran program computes the unbias probability of TASS output from CPMD/PLUMED run, which can be used to compute 1D/2D free energy via WHAM reweighting. This program also can compute 1D and 2D free enrgy using Mean Force method (MF). \

Basis Spline interpolation can be performed to find intermediate points in free energy .\

```diff
+ UPDATE    :: "1D and 2D free energy will be computed with Mean Force"
+ IMPORTANT :: "Input variables in the file are case sensitive (i.e. +-> can be turn on/off with upper/lower case.)"
[Ref : https://github.com/jacobwilliams/bspline-fortran]
```

# Modular Code Written by :- Rahul Verma
#---------------------------------------------- INPUT DESCRIPTION ---------------------------------------------------------------------------------------
```

NUMBER OF CV		 	# Total CV's Chosen in the Simulation	
5
CODE NAME			# Name of MD the Package (CPMD/PLUMED)
CPMD
NUMBER OF UMBRELLA		# Number of replica (total umbrella window during simulation)
22
UCV COLUMN			# Umbrella CV Index
1
MTD ENABLED			# IF Metadynamics performed durung simulation (y/n)
n
MTD CV COLUMN			# IF MTD ENABLED then Metadynamics CV Index [else=0]
2
MTD BIAS FACTOR			# Metadynamics Bias Factor
1200
SYSTEM TEMP			# Physical system Temperature
300
CV TEMP				# CV Temperature
1000
TMIN				# Minimum MD steps to Compute Probability
1000
TMAX				# Maximum MD steps to Compute Probability
14000
REWEIGHTING TOOL		# pmf/prob [pmf-->compute potential of mean force ; prob --> Unbias Probability]
pmf
PROBABILITY DIMENSION		# Dimension of Probability/Free Energy(MeanForce) to be generated [1/2]
2
PROBABILITY CV INDEX		# Probability/Free Energy(MeanForce) generated along CV indicis [1 --> 1D ; 1 2 --> 2D along CV1 and CV2 etc..]
1 3
STATISTICAL ERRORS BLOCK SIZE   # Compute Statistical Error in Block Data [optimal Block Size --> 4/5]
4
CV PRINT FREQUENCY		# Frequency of update in cvmdck_mtd file during run
1
MTD PRINT FREQUENCY		# Frequency of Hill Update during Metadynamics
10
GRIDS				# gridmin gridmax griddif for every CV [NCV]
1.8 6.0  0.01
0.5 5.0  0.05
1.0 8.0  0.01
1.0 10.0 0.01
1.0 8.0  0.01
```

#----------------
# HOW TO INSTALL
#----------------

```bash
./configure
[choose the options]
```

```Makefile
INSTALL :
Type...
make install   : create executabls
make bspline   : compile B-spline modules
make clean     : remove object and mod files
make distclean : clean the directory
```

```bash
How to Run -->
"-------------"
tass.x < $INPUT

A bash script is given along with the program (run.sh) 
Create execute permission by following command :
chmod 755 run.sh
./run.sh
```

# AUTHOUR

RAHUL VERMA \
DEPARTMENT OF CHEMISTRY \
IIT KANPUR, INDIA \
Email : vrahul@iitk.ac.in

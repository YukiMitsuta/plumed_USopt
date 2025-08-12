These codes are add-on programs for plumed (https://www.plumed.org/) to calculate the parameter optimization of the umbrella sampling method.

For more details, please read the article below.

Mitsuta, Y., & Asada, T. (2024). Parameter Optimization Method in Multidimensional Umbrella Sampling. Journal of Chemical Theory and Computation, 20(15), 6531-6548.

# How to use the umbrella sampling parameter optimization method. 

## 1. Installation 
This code works by adding it to plumed. First, download Plumed from the website. The current code has been confirmed to work with plumed-2.9.0. It will not work with versions before 2.8. Once downloaded, copy this program to the src/bias directory.
```
$ cp plumed_USopt/src/bias/* /path/to/plumed-2.9.0/src/bias 
```

This program also uses Eigen, a C++ matrix calculation library, so please install Eigen beforehand. If you are using Ubuntu, you can download it using apt. 
```
$ sudo apt install libeigen3-dev
```
 In addition, when installing plumed, you will need to specify the Eigen library in CPPFLAGS. 
```
$ ./configure --prefix=/installing/path/to/plumed-2.9.0 \
CPPFLAGS="-I/path/to/eigen3/include/eigen3/" \
(Other Settings) 
```
*If installed with Ubuntu's apt, the path is /usr/include/eigen3/. 

## 2. Tutorial (Alanine Dipeptide) 
You can use this code with any program that supports umbrella sampling in Plumed. Here I've attached a Gromacs calculation for alanine dipeptide as a tutorial file. The action uses OPTIMIZERESTRAINTHOTELLING. Its options are as follows: 
| Option  | Description |
| ------------- | ------------- |
| OFFDIAGONAL  | Whether to include off-diagonal terms in the bias potential.  |
| ARG  | CV to apply bias. |
| AT  | CV to apply bias. |
| TARGET  | Target point for optimization (optimized so that the mean point is at this position).  |
| KAPPA  | Initial value of the bias potential. Specify KAPPA0, KAPPA1... according to the matrix dimensions.  |
| OPTSTRIDE  | Optimization stride  |
| CVSTEPSIZE  | Optimization step size (learning rate) for CV  |
| KAPPASTEPSIZE  | Optimization step size (learning rate) for bias potential   |
| BETA1  | Stochastic optimization parameter (beta_1)  |
| BETA2  | Stochastic optimization parameter (beta_2)  |
| EPSILON  | Stochastic optimization parameter (epsilon)  |
| EXP_DECAYING_AVER  | Exponential decaying average term |
| HOTELLINGMAX  | Upper limit on Hotelling distance (specifying this parameter has been found not to improve optimization; turn it off by entering 99999)  |
| TARGETSIGMA  | Upper limit on optimization variance  |
| KAPIGNORECVPASTEPSIZE  | CV to ignore in optimization (enter 1.0 for the ignored dimension). |


Note that the dimension of TARGETSIGMA does not correspond to the dimension of CV, but rather to the order of the largest eigenvalues when the variance-covariance matrix is expanded. 
IGNORECV can cause optimization to fail to converge when there are upper or lower limits on the CV, so set it to ignore optimization near those limits (for example, the coordination number). 

The output options are as follows. Note that you must specify _label_._arg_\_cntr or _label_._arg1_\__arg2_\_kappa. 

| Option  | Description |
| ------------- | ------------- |
| _cntr  | Optimized reference point   |
| _kappa  | Optimized bias potential parameter   |
| _mean  |  Mean sampling point    |
| _Kgrad  | Optimized bias potential gradient    |
| _CVgrad  | Optimized reaction coordinate gradient    |


### Step 1: Calculate the optimization calculation. 
We have provided calljob.sh as the executable file for MD calculations. Modify it as needed. The script for performing the optimization calculation is as follows. 
```
$ gmx grompp -f ./npt.mdp -o ./npt.tpr -c ./min.gro -p ./topol.top
$ gmx mdrun -deffnm npt -plumed ./plumed_npt.dat
```
### Step 2: Analyze the optimization results and make plumed.dat. 
By printing the _cntr and _kappa values from the calculation results and analyzing their final time values, you can specify the bias potential of the optimized window. Here, we have provided a script called mkplumed.py, so please use it. Also, if you change the calculation system, please modify this script appropriately. 

### Step 3: Calculate window sampling. 
Once plumed.dat is created, all you need to do is run the main calculation to perform the umbrella sampling calculation. 
```
$ gmx grompp -f ./run.mdp -o ./run.tpr -c ./npt.gro -p ./topol.top 
$ gmx mdrun -deffnm run -plumed ./plumed.dat
```

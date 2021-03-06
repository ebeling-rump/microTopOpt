


# microTopOpt
Democratizing Microstructures - A Dolphin-Adjoint Code

# Install & Run

Step 1: Install Anaconda
https://www.anaconda.com/products/individual

Step 2: Update Anaconda, created a „fenicsproject“ virtual environment, install FEniCS, activate the virtual environment and install IPOPT

        conda update conda
        conda update anaconda
        conda create -n fenicsproject -c conda-forge fenics
        conda activate fenicsproject
        conda install -c conda-forge ipopt
       
Step 3: Navigate in the console to the python file and start the FEniCS program via
	
	      python microTopOpt.py
  
# Results

The first video shows the result of the program as is. The goal is to create a microstructure, which leads to a negative Poisson ratio of the homogenized macrostructure. In other words: If the structure is pulled apart horizontally, it will stretch vertically. 
The second video follows the same incentive. It is run with a larger Volume Constraint of 60% instead of 40%, the number of iteration steps is increased from 50 to 500 and the initial state has less holes.

![](git_vid.gif)

![](allcontrols_IPOPT_w1111_1.0_w1122_30_w2222_1.0_AT1111_0.2_AT1122_-0.1_AT2222_0.2_vc_0.6_GL_gamma_1e-05_GL_eps_1_niter_500_InEq_True_ndof_50_nu_T_-0.5_fst.gif)


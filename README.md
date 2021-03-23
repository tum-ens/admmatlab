# admMATLAB

*admMATLAB* is a MATLAB interface, which allows for decomposed optimization of any linear optimization problem. The decomposition is based on the consensus-based [alternating direction method of multipliers (ADMM)](https://en.wikipedia.org/wiki/Augmented_Lagrangian_method#Alternating_direction_method_of_multipliers).

## eXtremOS project
The *admMATLAB* interface has been developed within the scope of the BMWi (German Federal Ministry of Economics and Energy) funded research project eXtremOS. It is integrated to the energy system modelling framework [ISAaR](https://www.ffe.de/die-methoden/modelle/625-isaar-integriertes-simulationsmodell) of the [Forschungsstelle für Energiewirtschaft e.V. (FfE)](https://www.ffe.de/en), in order to regionally decompose the energy system models built with it in order to investigate the runtime and working memory benefits of decomposition.

## Contributors
The *admMATLAB* interface is developed by Soner Candas.

## Functionality
**Problem Input:**
 As input, *admMATLAB* interface takes an arbitrary linear optimization problem in its matrix/vector form, and applies the decomposition according to an accompanying *annotation vectors* of its variables. The *annotation vectors* are a pair of vectors each with the length of the number of variables in the original optimiziation problem, and denote which variables belong internally to certain clusters, and which are the *coupling variables* between which clusters. For example, an annotation such as:

 <img src="https://latex.codecogs.com/svg.latex?\Large&space; \begin{bmatrix}1 & 1 & 1 & 2 & 2 & 2 & 3 & 3 \\ 0 & 0 & 2 & 1 & 0 & 3 & 2 & 0\end{bmatrix}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}"  alt=""/> 

would define three clusters (1, 2 and 3), with:
- the first two variables being internal variables to cluster 1,
- the third and fourth variables being coupling variables between clusters 1 and 2, 
- the fifth variable being an internal variable to cluster 2,
- the sixth and seventh varaibles being coupling variables between clusters 2 and 3 and
- the eight variable being an internal variable to cluster 3.

The first and the second rows of this annotation matrix is defined as the variables <code>cluster_id_from</code> and <code>cluster_id_to</code>. 

These inputs (problem definition and annotation) are set in the <code>define_input()</code> function.

**ADMM input: **
 Besides the problem definition, ADMM-specific parameters have to be input on the  <code>runme</code>  script.

**ADMM modes: ** In this interface, four ADMM modes can be chosen from. These are categorized as follows:

1. Regular ADMM (sequential): Subproblems are created for each cluster, which are solved after one another during each iteration. To activate this mode, use <code>solver = 1</code> and <code>num_threads = 1</code> as the first and third input arguments of the <code>runme</code> function.   
2. Regular ADMM (parallel): Subproblems are created for each cluster, which are solved in parallel during each iteration. To activate this mode, use <code>solver = 1</code> and a <code>num_threads > 1</code> as the first and third input arguments of the <code>runme</code> function.
3. ADMM with bin packing algorithm (parallel): Groups of subproblems are collected in "bins". During each iteration, each subproblem within a given bin is solved after one another, and bins are solved in parallel to another. This method is developed as a countermeasure to the straggling effect, where the ADMM routine might lag because of a single, larger cluster requiring much longer to solve compared to the other clusters. To activate this mode, use <code>solver = 2</code> and a <code>num_threads > 1</code> as  the first and third input arguments of the <code>runme</code> function.
4. Asynchronous ADMM (parallel): Subproblems are created for each cluster, which are solved asynchronously. This way, each subproblem have their "local" iteration counters, and they can move onto the next iteration by receiving information from only one neighbor. Similar to the bin packing algorithm, this method is also developed as a countermeasure to the straggling effect. To activate this mode, use <code>solver = 3</code> and a <code>num_threads > 1</code> as the first and third input arguments of the <code>runme</code> function. For the details of this method, refer to the following paper:  [Asynchronous ADMM for Distributed Non-Convex Optimization in Power Systems](https://arxiv.org/abs/1710.08938)
   
|             | Regular ADMM (sequential) | Regular ADMM (parallel) | ADMM with bin packing | Asynchronous ADMM |
|-------------|:-------------------------:|:-----------------------:|:---------------------:|:-----------------:|
| solver      |             1             |            1            |           2           |         3         |
| num_threads |  higher than 1| higher than 1|  higher than 1|higher than 1|
   
**ADMM update methods: ** Moreover, five ADMM update methods can be chosen from. These are the following:
1. [Constant quadratic penalty](https://www.sciencedirect.com/science/article/pii/S0168202408700341): The quadratic penalty term is not adjusted between iterations (denoted in code as <code>Gabay_constant_rho</code>). To activate this update method, use <code>admm_option = 1</code> as the second input argument of <code>runme</code> function.   
2. [Adaptive quadratic penalty](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf): By using absolute residuals, the quadratic penalty term is adjusted between iterations to improve convergence (denoted in code as <code>Gabay_constant_rho</code>). To activate this update method, use <code>admm_option = 2</code> as the second input argument of <code>runme</code> function.  
3. [Local quadratic penalty](https://arxiv.org/abs/1506.08928): By using their respective absolute residuals, the quadratic penalty terms for each cluster are adjusted between iterations to improve convergence separately (denoted in code as <code>Gabay_constant_rho</code>). To activate this update method, use <code>admm_option = 3</code> as the second input argument of <code>runme</code> function.  
4. [Restart method](https://pdfs.semanticscholar.org/b1dc/0fa0edd9ccf77f6a3df48833f6e48fc40068.pdf?_ga=2.132469257.687946128.1616520109-1209273558.1616520109): To improve convergence, the quadratic penalty term is reset by a certain criterion (denoted in code as <code>Gabay_constant_rho</code>). To activate this update method, use <code>admm_option = 4</code> as the second input argument of <code>runme</code> function.  
5. [Relative residuals](https://arxiv.org/pdf/1704.06209.pdf): By using the relative residuals, the quadratic penalty term is adjusted between iterations to improve convergence (denoted in code as <code>Gabay_constant_rho</code>). To activate this update method, use <code>admm_option = 5</code> as the second input argument of <code>runme</code> function.  


## Requirements
This interface has been developed and tested on MATLAB R2019b. Additionally, the Optimization Toolbox and Parallel Computing Toolbox
have to be installed.

## Copyright
Copyright (C) 2021  TUM ENS

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

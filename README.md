# Parallel *O(N)* Fredkin-Koehler Solver

### Introduction

Fredkin-Koehler method is a numerical method to compute electrostatic/magnetic potential inside ferroelectric/ferromagnetic bodies. It is a hybrid method that combines finite element and boundary integral methods. The boundary integration in the original method has $O(N^2)$ scaling, where $N$ is number of unknowns in the discretized mesh. The governing equation is shown as follows, (to view the equations and this manual more clearly, please use the manual.html in **doc** folder)

$\nabla^2 \phi(\bold{x})=4 \pi \nabla \cdot \bold{M(\bold{x})} \quad$ for $\bold{x}\in R_m$,

$\nabla^2 \phi(\bold{x})=0 \qquad \qquad \quad$ for $\bold{x}\in R_e$,

where $R_m$ denotes the regions occupied by ferroelectric/ferromagnetic bodies, and $R_e$ denotes non-ferroelectric/ferromagnetic region complementary to $R_m$. One boundary condition is that the potential is zero at infinity; the other boundary condition is that the normal derivative of the potential is discontinuous at the surface of the bodies and equals to the polarization/magnetization density multiplied by normal vector at the surface. More details of the method can be found in the original paper [here](http://ieeexplore.ieee.org/document/106342/).

In a recent work, we have developed a parallel code that further improves the computational efficiency of the method. In particular, we have reduced the scaling of the boundary integration to $O(N)$ using a fast multipole method. The code is written on [libMesh](http://libmesh.github.io/) (finite element library) and [ScalFMM](https://gitlab.inria.fr/solverstack/ScalFMM) (fast multipole method). Please see the paper in the **Algorithm** section for more details of our work.


## Installation

*1. Install Moose*

Follow instructions on [GettingStarted](http://mooseframework.org/getting-started/) to install Moose Framework. If you meet any troubles during installation, please refer to [moose-users](https://groups.google.com/forum/#!forum/moose-users) forum for help.


<br>

*2. Install Parallel O(N) Fredkin-Koehler solver*

- Download the latest public version

    `git clone https://github.com/xikaij/fredkin-koehler.git`

- Compile ScalFMM

   `cd contrib`
   `./build_scalfmm`

    If you meet any trouble in compiling ScalFMM, please refer to Help and News section on ScalFMM's GitLab site [here](https://gitlab.inria.fr/solverstack/ScalFMM).

- Add path to ScalFMM into the Makefile in the top directory. Edit last two `app_INCLUDES` in the Makefile as follows:
 
    `app_INCLUDES       += -I/your/path/to/fredkin-koehler/contrib/scalfmm/Src`
    `app_INCLUDES       += -I/your/path/to/fredkin-koehler/contrib/scalfmm/Build/Src`

- Then 'make', in the top directory, to compile the code

   `make`


<br>

*3. Run example.*

An example is provided in the top directory. The input files are 

> `fredkin-koehler.i`
> `poisson.i`
> `laplace.i`

and the mesh file is `sphere_tet_approx_size0_05.e`. An unit magnetization density in the x-direction $\bold{M}=$ (1, 0, 0) set as the initial condition. Run the application as

`mpiexec -np NUM_OF_CORES ./fredkin-koehler -i fredkin-koehler.i`

Replace `NUM_OF_CORES` by the total number of cores you want to use.


## Machinery

The code uses `MultiApps` functionality in Moose. `fredkin-koehler.i` is the master app, while `poisson.i` and `laplace.i` are subapps.

The code first solves the Poisson equation in `poisson.i` on the volumetric mesh. Then the volumetric solution ($\phi_1$) is transfered from sub_app1 (`poisson.i`) to the master_app (`fredkin-koehler.i`), and is subsequently transfered from master_app (`fredkin-koehler.i`) to sub_app2 (`laplace.i`). The transfers are done using `MultiAppCopyTransfer` in Moose.

In `laplace.i`, a boundary mesh is extracted from the volumetric mesh and the boundary mesh is re-partitioned for much better load balancing during boundary integration. The boundary integral is accelerated by a fast multipole method, and results of the integration are assigned on the boundary of the volumetric mesh and serve as Dirichlet boundary condition for a subsequent Laplace solve on the volumetric mesh. The parts for extracting boundary mesh from volumetric mesh, boundary integration, and data transfers between boundary mesh and volumetric mesh is implemented in the `BoundaryIntegralFMM` object, which is derived from `GeneralUserObject` in Moose. 

The solution of the Laplace equation provides $\phi_2$, which is then transfered from sub_app2 (`laplace.i`) to the master_app (`fredkin-koehler.i`) and added upon values of $\phi_1$ in the master_app to get the final potential $\phi$. The solution transfer is done by `MultiAppAddTransfer`, which is borrowed from `MultiAppCopyTransfer`. The only difference between them is that `MultiAppCopyTransfer` uses libMesh function of solution.set(), while `MultiAppAddTransfer` uses solution.add().



## Algorithm

The method paper has details on the implementation and tests of this code: [An O(N) and parallel approach to integral problems by a kernel-independent fast multipole method: Application to polarization and magnetization of interacting particles. The Journal of Chemical Physics 145, 064307 (2016)](http://aip.scitation.org/doi/10.1063/1.4960436).



**Contributor**
-------------------------------------------

 - [Xikai Jiang](https://www.researchgate.net/profile/Xikai_Jiang), Institute for Molecular Engineering, University of Chicago. [EMAIL](xikai@uchicago.edu)


For more details on implementation and usage of this app, please email me at xikai@uchicago.edu. Comments on the app and this manual are also welcome.


**Project supervisor**
------

 - Olle G. Heinonen, Materials Science Division, Argonne National Laboratory


**License**
-------------------------------------------
The code is open-source and distributed under the GNU GPL license, and may not be used for any commercial or for-profit purposes without our permission. [Date: 01/30/2018]



## Create your own app

"Fork Fredkin-Koehler" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)



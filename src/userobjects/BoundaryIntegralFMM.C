#include "BoundaryIntegralFMM.h"

// libMesh headers
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/node.h"

// PETSc headers
#include "petscsys.h"

// ScalFMM headers
#include "Utils/FPoint.hpp"
#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Core/FFmmAlgorithmTsm.hpp"

// Global consts
const double PI_4 = 4.*M_PI;

using namespace libMesh;

template<>
InputParameters validParams<BoundaryIntegralFMM>()
{
  InputParameters params = validParams<GeneralUserObject>();

  params.set<std::string>("built_by_action") = "add_user_object";
  params.addRequiredParam<Real>("cx","x-coordinate for the center of FMM box");
  params.addRequiredParam<Real>("cy","y-coordinate for the center of FMM box");
  params.addRequiredParam<Real>("cz","z-coordinate for the center of FMM box");
  params.addRequiredParam<Real>("boxWidth","Width of the FMM cubic box");
  params.addRequiredParam<unsigned int>("TreeHeight","Octree height or level");

  return params;
}

BoundaryIntegralFMM::BoundaryIntegralFMM(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    _cx(getParam<Real>("cx")),
    _cy(getParam<Real>("cy")),
    _cz(getParam<Real>("cz")),
    _boxWidth(getParam<Real>("boxWidth")),
    _TreeHeight(getParam<unsigned int>("TreeHeight"))
{
}

void
BoundaryIntegralFMM::initialize()
{
}

void
BoundaryIntegralFMM::execute()
{
  // Get a constant reference to the mesh object
  const MeshBase & mesh_bs = _subproblem.mesh().getMesh();

  // The element dimension of boundary mesh
  const unsigned int dim = mesh_bs.mesh_dimension();

  // Get a reference to the EquationSystem
  EquationSystems & bs = _subproblem.es();

  // The boundary integral uses ExplicitSystem and
  // have only one system
  ExplicitSystem & system_bs = bs.get_system<ExplicitSystem> (0);

  std::vector<Number> global_potential(system_bs.solution->size());

  // Petsc performance log for boundary integral calculation
  // other than VecScatter
  int USER_EVENT;
  PetscLogEventRegister("BoundaryInt",0,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  // A reference to the DofMap object for this system.
  const DofMap & dof_map_bs = system_bs.get_dof_map();

  // Define dof_indices holder for phi2
  std::vector<dof_id_type> dof_indices_phi2;

  // Get system number and variable numbers
  const unsigned short int        system_number = system_bs.number();
  const unsigned short int variable_number_phi2 = 0; // phi2

  // Get a constant reference to variable phi2, get their number of components
  const Variable &                variable_phi2 = system_bs.variable(variable_number_phi2);
  const unsigned short int   variable_comp_phi2 = variable_phi2.n_components();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType bs_type = dof_map_bs.variable_type(variable_number_phi2);

  // Build a Finite Element object of the specified type. Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as an UniquePtr<FEBase>. This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> bs_face (FEBase::build(dim, bs_type));

  // Quadraure rule for surface integration with dimensionality
  // one less than the dimensionality of the element.
  QGauss qface(dim, FIFTH);

  // Tell FE object to use quadrature rule.
  bs_face->attach_quadrature_rule (&qface);

  // ScalFMM
  typedef double FReal;
  // In order to be used in a template, a constant value must be initialized
  const unsigned int ORDER = 5;
  FPoint<FReal> centerOfBox( _cx, _cy, _cz );
  const unsigned int SubTreeHeight = 1;

  // Particle, Leaf, Cell, Octree
  typedef FP2PParticleContainerIndexed<FReal>                 ContainerClass;
  typedef FTypedLeaf<FReal,ContainerClass>                    LeafClass;
  typedef FTypedChebCell<FReal,ORDER>                         CellClass;
  typedef FOctree<FReal,CellClass,ContainerClass,LeafClass>   OctreeClass;

  // Kernel D(1/r)/Dx2,D(1/r)/Dy2,D(1/r)/Dz2
  typedef FInterpMatrixKernelRx2<FReal> MatrixKernelClass1;
  typedef FInterpMatrixKernelRy2<FReal> MatrixKernelClass2;
  typedef FInterpMatrixKernelRz2<FReal> MatrixKernelClass3;

  // Three kernel classes
  typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass1,ORDER> KernelClass1;
  typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass2,ORDER> KernelClass2;
  typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass3,ORDER> KernelClass3;

  // Three kernels
  const   MatrixKernelClass1 MatrixKernel1;
  const   MatrixKernelClass2 MatrixKernel2;
  const   MatrixKernelClass3 MatrixKernel3;

  // Octrees 
  OctreeClass tree1(_TreeHeight,SubTreeHeight,_boxWidth,centerOfBox);
  OctreeClass tree2(_TreeHeight,SubTreeHeight,_boxWidth,centerOfBox);
  OctreeClass tree3(_TreeHeight,SubTreeHeight,_boxWidth,centerOfBox);

  FPoint<FReal> particlePosition;
  FSize indexPart = 0;

  // Node iterator for global mesh, because here we don't partition target points.
  MeshBase::const_node_iterator           nd = mesh_bs.nodes_begin();
  const MeshBase::const_node_iterator end_nd = mesh_bs.nodes_end();

  // Loop over all nodes - targets
  for ( ; nd != end_nd ; ++nd){   
      const Node* node_bs = *nd;

      // Target point coords
      const Real xt = (*node_bs)(0);
      const Real yt = (*node_bs)(1);
      const Real zt = (*node_bs)(2);
      particlePosition.setPosition( xt, yt, zt );

      // Insert into trees             particleType    index     physicalValue  pot forces 
      tree1.insert(particlePosition, FParticleType(1), indexPart,           0., 0., 0., 0., 0.);
      tree2.insert(particlePosition, FParticleType(1), indexPart,           0., 0., 0., 0., 0.);
      tree3.insert(particlePosition, FParticleType(1), indexPart,           0., 0., 0., 0., 0.);

      indexPart += 1;
  }

  const FSize nbTargets = indexPart;

  // Iterator el will iterate from the first to the last element on
  // this processor. Because here we partition the source points.
  const MeshBase::const_element_iterator end_el = mesh_bs.local_elements_end();
  MeshBase::const_element_iterator           el = mesh_bs.local_elements_begin();

  // Loop over all sources at quadrature points in every elements.
  // ++el requires an unnecessary temporary object.
  for ( ; el != end_el ; ++el){
      // Store a pointer to the element
      const Elem* elem_bs = *el;
  
      // The Jacobian * Quadrature Weight at the quadrature points on the face.
      const std::vector<Real>& JxW_face = bs_face->get_JxW();

      // The shape function at quadrature points
      const std::vector<std::vector<Real> >& phi = bs_face->get_phi();
  
      // The XYZ locations (in physical space) of the quadrature points on the face.
      const std::vector<Point >& qface_point = bs_face->get_xyz();
  
      // Tangent direction of xi and eta, cross product to get normal
      const std::vector<RealGradient >& qface_dxyzdxi  = bs_face->get_dxyzdxi();
      const std::vector<RealGradient >& qface_dxyzdeta = bs_face->get_dxyzdeta();
 
      // Compute the shape function values on the element face.
      bs_face->reinit(elem_bs);

      // Calculate the normal vector, we only need normal vector at one quadrature point because surface element is 2D
      // Cross product
      Real nx = qface_dxyzdxi[0](1)*qface_dxyzdeta[0](2)-qface_dxyzdxi[0](2)*qface_dxyzdeta[0](1);
      Real ny = qface_dxyzdxi[0](2)*qface_dxyzdeta[0](0)-qface_dxyzdxi[0](0)*qface_dxyzdeta[0](2);
      Real nz = qface_dxyzdxi[0](0)*qface_dxyzdeta[0](1)-qface_dxyzdxi[0](1)*qface_dxyzdeta[0](0);

      // Normalize the normal vector
      Real nunit = sqrt (nx*nx + ny*ny +nz*nz);
      nx = nx / nunit;
      ny = ny / nunit;
      nz = nz / nunit;

      // Global dof_indices for this element and variable phi2
      dof_map_bs.dof_indices(elem_bs, dof_indices_phi2, variable_number_phi2);
      // Number of dof indices for phi2 on this element
      // used to loop through all nodes to calculate phi2 on quadrature point
      const unsigned int n_phi2_dofs = dof_indices_phi2.size();

      // Loop over the face quadrature points for integration.
      for (unsigned int qp=0; qp<qface.n_points(); qp++){
          // Location of source points
          const Real x_qp = qface_point[qp](0);
          const Real y_qp = qface_point[qp](1);
          const Real z_qp = qface_point[qp](2);
 
          // Value of phi2 at quadrature point
          Real phi2_qp = 0.0;
          for (unsigned int l=0; l < n_phi2_dofs; l++){
              phi2_qp += phi[l][qp] * (*system_bs.current_local_solution)(dof_indices_phi2[l]);
          }

          Real phys = JxW_face[qp]*phi2_qp;

          particlePosition.setPosition( x_qp , y_qp , z_qp );
          // Insert into trees             particleType    index      physicalValue  pot forces
          tree1.insert(particlePosition, FParticleType(0), indexPart,       phys*nx, 0., 0., 0., 0.);
          tree2.insert(particlePosition, FParticleType(0), indexPart,       phys*ny, 0., 0., 0., 0.);
          tree3.insert(particlePosition, FParticleType(0), indexPart,       phys*nz, 0., 0., 0., 0.);

          indexPart += 1;
      }
  // End of looping over elements
  }

  const FSize nbTotal = indexPart;

  std::cout << "Particle insertion into octree done." << std::endl
            << nbTargets << " target particles, " << nbTotal-nbTargets << " source particles." << std::endl;

  // Apply kernels, here performs the compression and set M2L operators
  KernelClass1 kernel1(_TreeHeight,_boxWidth,centerOfBox,&MatrixKernel1);
  KernelClass2 kernel2(_TreeHeight,_boxWidth,centerOfBox,&MatrixKernel2);
  KernelClass3 kernel3(_TreeHeight,_boxWidth,centerOfBox,&MatrixKernel3);

  // FMM algorithm combine octree and kernel
  typedef FFmmAlgorithmTsm<OctreeClass,CellClass,ContainerClass,KernelClass1,LeafClass> FmmClass1;
  typedef FFmmAlgorithmTsm<OctreeClass,CellClass,ContainerClass,KernelClass2,LeafClass> FmmClass2;
  typedef FFmmAlgorithmTsm<OctreeClass,CellClass,ContainerClass,KernelClass3,LeafClass> FmmClass3;

  FmmClass1 algorithm1(&tree1, &kernel1);
  FmmClass2 algorithm2(&tree2, &kernel2);
  FmmClass3 algorithm3(&tree3, &kernel3);

  algorithm1.execute();
  algorithm2.execute();
  algorithm3.execute();

  // Store particle information
  struct TestParticle{
     FReal potential;
  };
  TestParticle* const particles = new TestParticle[nbTotal];

  // Get surface integral at each target point
  tree1.forEachLeaf([&](LeafClass* leaf){
      const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
          const FSize indexPartOrig = indexes[idxPart];
          particles[indexPartOrig].potential += potentials[idxPart];
          global_potential[indexPartOrig] = potentials[idxPart];
      }
  });

  tree2.forEachLeaf([&](LeafClass* leaf){
      const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
          const FSize indexPartOrig = indexes[idxPart];
          particles[indexPartOrig].potential += potentials[idxPart];
          global_potential[indexPartOrig] += potentials[idxPart];
      }
  });

  tree3.forEachLeaf([&](LeafClass* leaf){
      const FReal*const potentials = leaf->getTargets()->getPotentials();
      const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
      const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

      for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
          const FSize indexPartOrig = indexes[idxPart];
          particles[indexPartOrig].potential += potentials[idxPart];
          global_potential[indexPartOrig] += potentials[idxPart];
      }
  });

  // Get contribution from all source points among all processors.
  mesh_bs.comm().sum(global_potential);

  // Reset nd and particle index, assign values to nodes at this processor
  indexPart = 0;
  nd = mesh_bs.local_nodes_begin();
  const MeshBase::const_node_iterator end_nd_local = mesh_bs.local_nodes_end();

  for ( ; nd != end_nd_local ; ++nd){
      const Node* node_bs = *nd;
 
      // Dof_index for each node
      const dof_id_type node_dof_index_phi2 = node_bs->dof_number(system_number,
                                                                  variable_number_phi2,
                                                                  variable_comp_phi2-1);

      // Integral value, only works for smooth surface,
      // minus SIGN get (R_target-R_source) in ScalFMM
      Real bi_value = -global_potential[node_bs->id()]/PI_4
                      -1./2.*(*system_bs.solution)(node_dof_index_phi2);
 
      // Boundary integral value to solution vector
      system_bs.solution->set(node_dof_index_phi2, bi_value);
 
      indexPart += 1;
  }

  system_bs.solution->close();

  // Petsc performance log
  PetscLogEventEnd(USER_EVENT,0,0,0,0);
}

void
BoundaryIntegralFMM::finalize()
{
}

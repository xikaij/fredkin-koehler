// See LICENCE file at project root
//
// #author P. Blanchard
// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "Files/FFmaGenericLoader.hpp"


#include "Kernels/Interpolation/FInterpMatrixKernel_TensorialInteractions.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebTensorialKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
#include "Utils/FMemUtils.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

/**
 * This program runs the FMM Algorithm with the Chebyshev kernel and compares the results with a direct computation.
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Run the chebyshev FMM kernel and compare the accuracy with a direct computation.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbThreads);

    typedef double FReal;
  const char* const filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
  const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 3);
  const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
  const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);

#ifdef _OPENMP
  omp_set_num_threads(NbThreads);
  std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
  std::cout << "\n>> Sequential version.\n" << std::
#endif

    // init timer
    FTic time;

    // typedefs and infos
    typedef FInterpMatrixKernel_R_IJ<FReal> MatrixKernelClass;
    MatrixKernelClass::printInfo();

    // useful features of matrix kernel
    const unsigned int NPV  = MatrixKernelClass::NPV;
    const unsigned int NPOT = MatrixKernelClass::NPOT;
    const unsigned int NRHS = MatrixKernelClass::NRHS;
    const unsigned int NLHS = MatrixKernelClass::NLHS;

    const FReal CoreWidth = 0.;
    std::cout << "Core width: a=" << CoreWidth << std::endl;
    std::cout << std::endl;

    // Build matrix kernel
    const MatrixKernelClass MatrixKernel(CoreWidth);

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        FReal forces[3][NPOT];
        FReal physicalValue[NPV];
        FReal potential[NPOT];
    };

  // open particle file
  FFmaGenericLoader<FReal> loader(filename);
  if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

  TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
  for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
    FPoint<FReal> position;
    FReal physicalValue = 0.0;
    loader.fillParticle(&position,&physicalValue);
    // get copy
    particles[idxPart].position       = position;
    // Set physical values
    for(unsigned idxPV = 0; idxPV<NPV;++idxPV){
    //   Either copy same physical value in each component
      particles[idxPart].physicalValue[idxPV]  = physicalValue; 
    // ... or set random value
//      particles[idxPart].physicalValue[idxPV]  = physicalValue*FReal(rand())/FReal(RAND_MAX);
    }

    for(unsigned idxPot = 0; idxPot<NPOT;++idxPot){
      particles[idxPart].potential[idxPot]      = 0.0;
      particles[idxPart].forces[0][idxPot]      = 0.0;
      particles[idxPart].forces[1][idxPot]      = 0.0;
      particles[idxPart].forces[2][idxPot]      = 0.0;
    }
  }

  ////////////////////////////////////////////////////////////////////

  { // begin direct computation
    std::cout << "\nDirect Computation ... " << std::endl;
    time.tic();
    {
      for(FSize idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
        for(FSize idxOther =  idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
          FP2P::MutualParticlesKIJ(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                   particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                   particles[idxTarget].forces[0], particles[idxTarget].forces[1],
                                   particles[idxTarget].forces[2], particles[idxTarget].potential,
                                   particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                   particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                   particles[idxOther].forces[0], particles[idxOther].forces[1],
                                   particles[idxOther].forces[2], particles[idxOther].potential,
                                   &MatrixKernel);
        }
      }
    }
    time.tac();
    std::cout << "Done  " << "(@Direct Computation = "
              << time.elapsed() << "s)." << std::endl;

  } // end direct computation

  ////////////////////////////////////////////////////////////////////

  {	// begin Chebyshev kernel

    // accuracy
    const unsigned int ORDER = 5 ;
    // set box width extension
    // ... either deduce from element size
  const FReal LeafCellWidth = FReal(loader.getBoxWidth()) / FReal(FMath::pow(2.,TreeHeight-1));
  const FReal ElementSize = LeafCellWidth / FReal(3.);  
  const FReal BoxWidthExtension = 0.*ElementSize; // depends on type of element
    // ... or set to arbitrary value (0. means no extension)
//    const FReal BoxWidthExtension = FReal(0.); 

    std::cout << "LeafCellWidth=" << LeafCellWidth 
              << ", BoxWidthExtension=" << BoxWidthExtension <<std::endl;

    // stop execution if interactions are homog and box extension is required 
    if(MatrixKernelClass::Type==HOMOGENEOUS && BoxWidthExtension>0.)
      throw std::runtime_error("Extension of box width is not yet supported for homogeneous kernels! Work-around: artificially set Type to NON_HOMOGENEOUS.");

    typedef FP2PParticleContainerIndexed<FReal, NRHS,NLHS> ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
    typedef FChebCell<FReal,ORDER,NRHS,NLHS> CellClass;
    typedef FOctree<FReal,CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FChebTensorialKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
    typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
    //  typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
      std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                << " particles ..." << std::endl;
      std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
      time.tic();

      for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        // put in tree
        // PB: here we have to know NPV...
        if(NPV==1)
          tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue[0]);
        else if(NPV==3)
          tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue[0], particles[idxPart].physicalValue[1], particles[idxPart].physicalValue[2]);
        else 
          std::runtime_error("Insert correct number of physical values!");

      }

      time.tac();
      std::cout << "Done  " << "(@Creating and Inserting Particles = "
                << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
      std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;
      time.tic();
      KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel, BoxWidthExtension);
      FmmClass algorithm(&tree, &kernels);
      algorithm.execute();
      time.tac();
      std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------


    { // -----------------------------------------------------
      std::cout << "\nError computation ... " << std::endl;
      FMath::FAccurater<FReal> potentialDiff[NPOT];
      FMath::FAccurater<FReal> fx[NPOT], fy[NPOT], fz[NPOT];

      FReal checkPotential[20000][NPOT];
      FReal checkfx[20000][NPOT];

      { // Check that each particle has been summed with all other

        tree.forEachLeaf([&](LeafClass* leaf){
            for(unsigned idxPot = 0; idxPot<NPOT;++idxPot){

              const FReal*const potentials = leaf->getTargets()->getPotentials(idxPot);
              const FReal*const forcesX = leaf->getTargets()->getForcesX(idxPot);
              const FReal*const forcesY = leaf->getTargets()->getForcesY(idxPot);
              const FReal*const forcesZ = leaf->getTargets()->getForcesZ(idxPot);
              const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
              const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

              for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const FSize indexPartOrig = indexes[idxPart];

                //PB: store potential in array[nbParticles]
                checkPotential[indexPartOrig][idxPot]=potentials[idxPart];
                checkfx[indexPartOrig][idxPot]=forcesX[idxPart];

                // update accuracy
                potentialDiff[idxPot].add(particles[indexPartOrig].potential[idxPot],potentials[idxPart]);
                fx[idxPot].add(particles[indexPartOrig].forces[0][idxPot],forcesX[idxPart]);
                fy[idxPot].add(particles[indexPartOrig].forces[1][idxPot],forcesY[idxPart]);
                fz[idxPot].add(particles[indexPartOrig].forces[2][idxPot],forcesZ[idxPart]);
              }
            }// NPOT
          });
      }

      // Print for information
      std::cout << "\nRelative Inf/L2 errors: " << std::endl;
      std::cout << "  Potential: " << std::endl;
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << potentialDiff[idxPot].getRelativeInfNorm() << ", " 
                  << potentialDiff[idxPot].getRelativeL2Norm() 
                  << std::endl;
      }
      std::cout << std::endl;
      std::cout << "  Fx: " << std::endl; 
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << fx[idxPot].getRelativeInfNorm() << ", " 
                  << fx[idxPot].getRelativeL2Norm()
                  << std::endl;
      }
      std::cout  << std::endl;
      std::cout << "  Fy: " << std::endl; 
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << fy[idxPot].getRelativeInfNorm() << ", " 
                  << fy[idxPot].getRelativeL2Norm()
                  << std::endl;
      }
      std::cout  << std::endl;
      std::cout << "  Fz: " << std::endl; 
      for(unsigned idxPot = 0; idxPot<NPOT;++idxPot) {
        std::cout << "    " << idxPot << ": "
                  << fz[idxPot].getRelativeInfNorm() << ", " 
                  << fz[idxPot].getRelativeL2Norm()
                  << std::endl;
      }
      std::cout << std::endl;

    } // -----------------------------------------------------

  } // end Chebyshev kernel

  return 0;
}

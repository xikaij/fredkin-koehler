// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_MPI
// @FUSE_FFT
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "../../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"

#include "../../Src/Kernels/Uniform/FUnifCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifKernel.hpp"


#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"

#include "../../Src/Utils/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
 * This program runs the FMM Algorithm Distributed with the uniform kernel
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Test with MPI the uniform FMM and compare it to the direct computation for debugging purpose.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbThreads);

    typedef double FReal;
    const unsigned int ORDER = 5;

    typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;

//    typedef FRotationCell<FReal,ORDER> CellClass;
 //   typedef FRotationKernel<FReal,CellClass,ContainerClass,ORDER> KernelClass;


        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FUnifCell<FReal,ORDER> CellClass;
        typedef FUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;

         typedef FOctree<FReal,CellClass,ContainerClass,LeafClass> OctreeClass;
         typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

    FMpi app(argc,argv);

    const FSize nbParticles                   = FParameters::getValue(argc,argv, FParameterDefinitions::NbParticles.options, 10000000ULL);
    const unsigned int TreeHeight      = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads       = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, omp_get_max_threads());
    FTic time;

    std::cout << ">> This executable has to be used to test Proc Uniform Algorithm. \n";

    omp_set_num_threads(NbThreads);
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    // open particle file
    std::cout << "Creating : " << nbParticles << "\n" << std::endl;
    FRandomLoader<FReal> loader(nbParticles, 1.0, FPoint<FReal>(0,0,0), app.global().processId());

    OctreeClass tree(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());

    time.tic();
    TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
    memset(particles,0,(unsigned int) (sizeof(TestParticle)* loader.getNumberOfParticles()));
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position);
        particles[idxPart].physicalValue = 1.0;
    }

    FVector<TestParticle> finalParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                loader.getNumberOfParticles(),
                                                                tree.getBoxCenter(),
                                                                tree.getBoxWidth(),tree.getHeight(),
                                                                &finalParticles, &balancer);
    { // -----------------------------------------------------
        std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
        std::cout << "For a total of " << loader.getNumberOfParticles() * app.global().processCount() << " particles ..." << std::endl;
        std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
        time.tic();

        for(FSize idxPart = 0 ; idxPart < finalParticles.getSize() ; ++idxPart){
            // put in tree
            tree.insert(finalParticles[idxPart].position, idxPart, finalParticles[idxPart].physicalValue);
        }

        time.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------

    delete[] particles;
    particles = 0;

    { // -----------------------------------------------------
        std::cout << "Uniform FMM (ORDER="<< ORDER << ") ... " << std::endl;
        time.tic();
   //     KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
        // Create Matrix Kernel
          const MatrixKernelClass MatrixKernel;

        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        FmmClass algorithm(app.global(),&tree, &kernels);
        time.tac();
        std::cout << "Done  " << "(@Init = " << time.elapsed() << "s)." << std::endl;
        time.tic();
        algorithm.execute();
        time.tac();
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------


    return 0;
}

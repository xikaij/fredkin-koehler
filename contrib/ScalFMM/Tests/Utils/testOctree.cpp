// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"
#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
* In this file we show how to use octree
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Show how to use an octree (only the code is interesting)",
                         FParameterDefinitions::NbParticles);

    typedef double FReal;
    typedef FBasicParticleContainer<FReal,0,FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use the Octree\n";
    //////////////////////////////////////////////////////////////
    const FSize NbPart = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000000));
    FTic counter;

    FRandomLoader<FReal> loader(NbPart, 1, FPoint<FReal>(0.5,0.5,0.5), 1);
    OctreeClass tree(10, 3, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------
    std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
    counter.tic();

    FPoint<FReal> particlePosition;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particlePosition);
        tree.insert(particlePosition);
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
    // -----------------------------------------------------


    return 0;
}




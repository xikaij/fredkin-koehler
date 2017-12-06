// See LICENCE file at project root

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"
#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Arranger/FOctreeArranger.hpp"
#include "../../Src/Arranger/FArrangerPeriodic.hpp"
#include "../../Src/Arranger/FBasicParticleContainerIndexedMover.hpp"

#include "../../Src/Utils/FParameterNames.hpp"



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Put the particles into a tree, then change the position of some particles and update the tree.\n"
                         "This method should be used to avoid the tree reconstruction.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight);

    typedef double FReal;

    typedef FBasicCell                      CellClass;
    typedef FP2PParticleContainerIndexed<FReal>      ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    typedef FBasicParticleContainerIndexedMover<FReal, OctreeClass, ContainerClass> MoverClass;
    typedef FArrangerPeriodic<FReal, OctreeClass, ContainerClass, MoverClass> ArrangerClassPeriodic;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FSize NbPart           = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(2000000));

    FTic counter;

    srand48 ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal BoxWidth = 1.0;
    const FReal BoxCenter = 0.5;

    OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, FPoint<FReal>(BoxCenter,BoxCenter,BoxCenter));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();


    {

        FPoint<FReal> particleToFill;
        for(FSize idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particleToFill.setPosition(
                        (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)),
                        (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)),
                        (BoxWidth*FReal(drand48())) + (BoxCenter-(BoxWidth/2)));
            tree.insert(particleToFill,idxPart);
        }
    }
        counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    ArrangerClassPeriodic arrange(&tree);

    // In order to test multiple displacement of particles
    //    for(int ite=0 ; ite<10 ; ++ite){
    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    { // Create new position for each particles
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass* particles = octreeIterator.getCurrentListTargets();
            for(FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){
                particles->getWPositions()[0][idxPart] = (FReal(drand48()))*BoxWidth*4 + (BoxCenter-(BoxWidth/2));
                particles->getWPositions()[1][idxPart] = (FReal(drand48()))*BoxWidth*4 + (BoxCenter-(BoxWidth/2));
                particles->getWPositions()[2][idxPart] = (FReal(drand48()))*BoxWidth*4 + (BoxCenter-(BoxWidth/2));
            }
        } while(octreeIterator.moveRight());
    }

    counter.tac();
    std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Arrange ..." << std::endl;
    counter.tic();

    //FOctreeArranger<FReal,OctreeClass, ContainerClass, TestParticle, Converter<TestParticle> > arrange(&tree);
    arrange.rearrange();


    counter.tac();
    std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //    }
    std::cout << "Test ..." << std::endl;
    counter.tic();

    { // Check that each particle has been put into the right leaf
        long counterPart = 0;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            const MortonIndex leafIndex = octreeIterator.getCurrentGlobalIndex();

            ContainerClass* particles = octreeIterator.getCurrentListTargets();
            for(FSize idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){
                const FPoint<FReal> particlePosition( particles->getWPositions()[0][idxPart],
                                               particles->getWPositions()[1][idxPart],
                                               particles->getWPositions()[2][idxPart]);

                const MortonIndex particleIndex = tree.getMortonFromPosition( particlePosition );
                if( leafIndex != particleIndex){
                    std::cout << "Index problem, should be " << leafIndex <<
                                 " particleIndex "<< particleIndex << std::endl;
                }

            }

            counterPart += octreeIterator.getCurrentListTargets()->getNbParticles();
            if(octreeIterator.getCurrentListTargets()->getNbParticles() == 0){
                std::cout << "Problem, leaf is empty at index " << leafIndex << std::endl;
            }
        } while(octreeIterator.moveRight());

        if( counterPart != NbPart ){
            std::cout <<"Wrong particles number, should be " << NbPart << " but is " << counterPart << std::endl;
        }
    }

    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator(&tree);
        OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = NbLevels - 1;
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            // for each cells
            do{
                int countChild = 0;
                CellClass** const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                    if( child[idxChild] ){
                        countChild += 1;
                    }
                }

                if(countChild == 0){
                    std::cout << "Problem at level " << idxLevel << " cell has no child " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Test = " << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}

// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP


#include "../Utils/FAssert.hpp"
#include "FAbstractLeaf.hpp"
#include "FParticleType.hpp"


// Clang warns that args is an unused parameter, maybe it does not understand variadic templates?
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTypedLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particles
* (source XOR target) or simple system (source AND target).
*
* Particles should be typed to enable targets/sources difference.
*/
template<class FReal, class ContainerClass>
class FTypedLeaf  : public FAbstractLeaf<FReal, ContainerClass> {
    ContainerClass sources; //< The sources containers
    ContainerClass targets; //< The targets containers

public:

    /** Default destructor */
    virtual ~FTypedLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param inParticlePosition the position of the new particle
        * @param isTarget bool to know if it is a target
        * followed by other param given by the user
        */
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FParticleType type, Args ... args){
        if(type == FParticleTypeTarget) targets.push(inParticlePosition, FParticleTypeTarget, args...);
        else sources.push(inParticlePosition, FParticleTypeSource, args...);
    }

    /**
     * To add a new particle in the leaf
     * @param inParticlePosition the position of the new particle
     * @param isTarget bool to know if it is a target
     * followed by other param given by the user
     */
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, Args ... args){
        FAssert(0,"Error : cannot push a particle without specifying type (src/tgt)");
    }




    /**
    * To get all the sources in a leaf
    * @return a pointer to the list of particles that are sources
    */
    ContainerClass* getSrc() {
        return &this->sources;
    }

    /**
    * To get all the target in a leaf
    * @return a pointer to the list of particles that are targets
    */
    ContainerClass* getTargets() {
        return &this->targets;
    }

};


#endif //FTYPEDLEAF_HPP

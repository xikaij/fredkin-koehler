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
#ifndef FP2P_HPP
#define FP2P_HPP

namespace FP2P {

/**
   * @brief MutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   *
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   * @param MatrixKernel pointer to an interaction kernel evaluator
   */
template <class FReal, typename MatrixKernelClass>
inline void MutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                            FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                            const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                            FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                            const MatrixKernelClass *const MatrixKernel){

    // Compute kernel of interaction...
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[1];
    FReal dKxy[3];
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);
    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX) += dKxy[0] * coef;
    (*targetForceY) += dKxy[1] * coef;
    (*targetForceZ) += dKxy[2] * coef;
    (*targetPotential) += ( Kxy[0] * sourcePhysicalValue );

    (*sourceForceX) -= dKxy[0] * coef;
    (*sourceForceY) -= dKxy[1] * coef;
    (*sourceForceZ) -= dKxy[2] * coef;
    (*sourcePotential) += ( Kxy[0] * targetPhysicalValue );
}

/**
   * @brief NonMutualParticles (generic version)
   * P2P mutual interaction,
   * this function computes the interaction for 2 particles.
   *
   * Formulas are:
   * \f[
   * F = - q_1 * q_2 * grad K{12}
   * P_1 = q_2 * K{12} ; P_2 = q_1 * K_{12}
   * \f]
   * In details for \f$K(x,y)=1/|x-y|=1/r\f$ :
   * \f$ P_1 = \frac{ q_2 }{ r } \f$
   * \f$ P_2 = \frac{ q_1 }{ r } \f$
   * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } \f$
   */
template <class FReal, typename MatrixKernelClass>
inline void NonMutualParticles(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal sourcePhysicalValue,
                               const FReal targetX,const FReal targetY,const FReal targetZ, const FReal targetPhysicalValue,
                               FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                               const MatrixKernelClass *const MatrixKernel){

    // Compute kernel of interaction...
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[1];
    FReal dKxy[3];
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);
    FReal coef = (targetPhysicalValue * sourcePhysicalValue);

    (*targetForceX) += dKxy[0] * coef;
    (*targetForceY) += dKxy[1] * coef;
    (*targetForceZ) += dKxy[2] * coef;
    (*targetPotential) += ( Kxy[0] * sourcePhysicalValue );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tensorial Matrix Kernels: K_IJ / p_i=\sum_j K_{ij} w_j
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
   * @brief MutualParticlesKIJ
   * @param sourceX
   * @param sourceY
   * @param sourceZ
   * @param sourcePhysicalValue
   * @param sourceForceX
   * @param sourceForceY
   * @param sourceForceZ
   * @param sourcePotential
   * @param targetX
   * @param targetY
   * @param targetZ
   * @param targetPhysicalValue
   * @param targetForceX
   * @param targetForceY
   * @param targetForceZ
   * @param targetPotential
   */
template<class FReal, typename MatrixKernelClass>
inline void MutualParticlesKIJ(const FReal sourceX,const FReal sourceY,const FReal sourceZ, const FReal* sourcePhysicalValue,
                               FReal* sourceForceX, FReal* sourceForceY, FReal* sourceForceZ, FReal* sourcePotential,
                               const FReal targetX,const FReal targetY,const FReal targetZ, const FReal* targetPhysicalValue,
                               FReal* targetForceX, FReal* targetForceY, FReal* targetForceZ, FReal* targetPotential,
                               const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

    // evaluate kernel and its partial derivatives
    const FPoint<FReal> sourcePoint(sourceX,sourceY,sourceZ);
    const FPoint<FReal> targetPoint(targetX,targetY,targetZ);
    FReal Kxy[ncmp];
    FReal dKxy[ncmp][3];
    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

    for(unsigned int i = 0 ; i < 3 ; ++i){
        for(unsigned int j = 0 ; j < 3 ; ++j){

            // update component to be applied
            const int d = applyTab[i*3+j];

            // forces prefactor
            const FReal coef = -(targetPhysicalValue[j] * sourcePhysicalValue[j]);

            targetForceX[i] += dKxy[d][0] * coef;
            targetForceY[i] += dKxy[d][1] * coef;
            targetForceZ[i] += dKxy[d][2] * coef;
            targetPotential[i] += ( Kxy[d] * sourcePhysicalValue[j] );

            sourceForceX[i] -= dKxy[d][0] * coef;
            sourceForceY[i] -= dKxy[d][1] * coef;
            sourceForceZ[i] -= dKxy[d][2] * coef;
            sourcePotential[i] += ( Kxy[d] * targetPhysicalValue[j] );

        }// j
    }// i

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tensorial Matrix Kernels: K_IJ
// TODO: Implement SSE and AVX variants then move following FullMutualKIJ and FullRemoteKIJ to FP2P.hpp
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief FullMutualKIJ
 */
template <class FReal, class ContainerClass, typename MatrixKernelClass>
inline void FullMutualKIJ(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                          const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            if( inNeighbors[idxNeighbors] ){
                const FSize nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
                const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
                const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
                const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){

                    // evaluate kernel and its partial derivatives
                    const FPoint<FReal> sourcePoint(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
                    const FPoint<FReal> targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
                    FReal Kxy[ncmp];
                    FReal dKxy[ncmp][3];
                    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

                    for(unsigned int i = 0 ; i < 3 ; ++i){
                        FReal*const targetsPotentials = inTargets->getPotentials(i);
                        FReal*const targetsForcesX = inTargets->getForcesX(i);
                        FReal*const targetsForcesY = inTargets->getForcesY(i);
                        FReal*const targetsForcesZ = inTargets->getForcesZ(i);
                        FReal*const sourcesPotentials = inNeighbors[idxNeighbors]->getPotentials(i);
                        FReal*const sourcesForcesX = inNeighbors[idxNeighbors]->getForcesX(i);
                        FReal*const sourcesForcesY = inNeighbors[idxNeighbors]->getForcesY(i);
                        FReal*const sourcesForcesZ = inNeighbors[idxNeighbors]->getForcesZ(i);

                        for(unsigned int j = 0 ; j < 3 ; ++j){
                            const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
                            const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

                            // update component to be applied
                            const int d = applyTab[i*3+j];

                            // forces prefactor
                            FReal coef = -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                            targetsForcesX[idxTarget] += dKxy[d][0] * coef;
                            targetsForcesY[idxTarget] += dKxy[d][1] * coef;
                            targetsForcesZ[idxTarget] += dKxy[d][2] * coef;
                            targetsPotentials[idxTarget] += ( Kxy[d] * sourcesPhysicalValues[idxSource] );

                            sourcesForcesX[idxSource] -= dKxy[d][0] * coef;
                            sourcesForcesY[idxSource] -= dKxy[d][1] * coef;
                            sourcesForcesZ[idxSource] -= dKxy[d][2] * coef;
                            sourcesPotentials[idxSource] += Kxy[d] * targetsPhysicalValues[idxTarget];

                        }// j
                    }// i
                }
            }
        }
    }
}

template <class FReal, class ContainerClass, typename MatrixKernelClass>
inline void InnerKIJ(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        for(FSize idxSource = idxTarget + 1 ; idxSource < nbParticlesTargets ; ++idxSource){

            // evaluate kernel and its partial derivatives
            const FPoint<FReal> sourcePoint(targetsX[idxSource],targetsY[idxSource],targetsZ[idxSource]);
            const FPoint<FReal> targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
            FReal Kxy[ncmp];
            FReal dKxy[ncmp][3];
            MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

            for(unsigned int i = 0 ; i < 3 ; ++i){
                FReal*const targetsPotentials = inTargets->getPotentials(i);
                FReal*const targetsForcesX = inTargets->getForcesX(i);
                FReal*const targetsForcesY = inTargets->getForcesY(i);
                FReal*const targetsForcesZ = inTargets->getForcesZ(i);

                for(unsigned int j = 0 ; j < 3 ; ++j){
                    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);

                    // update component to be applied
                    const int d = applyTab[i*3+j];

                    // forces prefactor
                    const FReal coef = -(targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);

                    targetsForcesX[idxTarget] += dKxy[d][0] * coef;
                    targetsForcesY[idxTarget] += dKxy[d][1] * coef;
                    targetsForcesZ[idxTarget] += dKxy[d][2] * coef;
                    targetsPotentials[idxTarget] += ( Kxy[d] * targetsPhysicalValues[idxSource] );

                    targetsForcesX[idxSource] -= dKxy[d][0] * coef;
                    targetsForcesY[idxSource] -= dKxy[d][1] * coef;
                    targetsForcesZ[idxSource] -= dKxy[d][2] * coef;
                    targetsPotentials[idxSource] += Kxy[d] * targetsPhysicalValues[idxTarget];
                }// j
            }// i

        }
    }
}

/**
   * @brief FullRemoteKIJ
   */
template <class FReal, class ContainerClass, typename MatrixKernelClass>
inline void FullRemoteKIJ(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                          const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    // get information on tensorial aspect of matrix kernel
    const int ncmp = MatrixKernelClass::NCMP;
    const int applyTab[9] = {0,1,2,
                             1,3,4,
                             2,4,5};

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            if( inNeighbors[idxNeighbors] ){
                const FSize nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
                const FReal*const sourcesX = inNeighbors[idxNeighbors]->getPositions()[0];
                const FReal*const sourcesY = inNeighbors[idxNeighbors]->getPositions()[1];
                const FReal*const sourcesZ = inNeighbors[idxNeighbors]->getPositions()[2];

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){

                    // evaluate kernel and its partial derivatives
                    const FPoint<FReal> sourcePoint(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource]);
                    const FPoint<FReal> targetPoint(targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget]);
                    FReal Kxy[ncmp];
                    FReal dKxy[ncmp][3];
                    MatrixKernel->evaluateBlockAndDerivative(sourcePoint,targetPoint,Kxy,dKxy);

                    for(unsigned int i = 0 ; i < 3 ; ++i){
                        FReal*const targetsPotentials = inTargets->getPotentials(i);
                        FReal*const targetsForcesX = inTargets->getForcesX(i);
                        FReal*const targetsForcesY = inTargets->getForcesY(i);
                        FReal*const targetsForcesZ = inTargets->getForcesZ(i);

                        for(unsigned int j = 0 ; j < 3 ; ++j){
                            const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues(j);
                            const FReal*const sourcesPhysicalValues = inNeighbors[idxNeighbors]->getPhysicalValues(j);

                            // update component to be applied
                            const int d = applyTab[i*3+j];

                            // forces prefactor
                            const FReal coef = -(targetsPhysicalValues[idxTarget] * sourcesPhysicalValues[idxSource]);

                            targetsForcesX[idxTarget] += dKxy[d][0] * coef;
                            targetsForcesY[idxTarget] += dKxy[d][1] * coef;
                            targetsForcesZ[idxTarget] += dKxy[d][2] * coef;
                            targetsPotentials[idxTarget] += ( Kxy[d] * sourcesPhysicalValues[idxSource] );

                        }// j
                    }// i

                }
            }
        }
    }
}


template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                              const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsForcesX = inTargets->getForcesX();
    FReal*const targetsForcesY = inTargets->getForcesY();
    FReal*const targetsForcesZ = inTargets->getForcesZ();
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+NbFRealInComputeClass-1)/NbFRealInComputeClass;
            const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const ComputeClass*const sourcesX = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[0];
            const ComputeClass*const sourcesY = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[1];
            const ComputeClass*const sourcesZ = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[2];
            ComputeClass*const sourcesForcesX = (ComputeClass*)inNeighbors[idxNeighbors]->getForcesX();
            ComputeClass*const sourcesForcesY = (ComputeClass*)inNeighbors[idxNeighbors]->getForcesY();
            ComputeClass*const sourcesForcesZ = (ComputeClass*)inNeighbors[idxNeighbors]->getForcesZ();
            ComputeClass*const sourcesPotentials = (ComputeClass*)inNeighbors[idxNeighbors]->getPotentials();

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
                const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
                const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
                const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
                ComputeClass  tfx = FMath::Zero<ComputeClass>();
                ComputeClass  tfy = FMath::Zero<ComputeClass>();
                ComputeClass  tfz = FMath::Zero<ComputeClass>();
                ComputeClass  tpo = FMath::Zero<ComputeClass>();

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    ComputeClass Kxy[1];
                    ComputeClass dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource],
                                                             tx,ty,tz,Kxy,dKxy);
                    const ComputeClass coef = (tv * sourcesPhysicalValues[idxSource]);

                    dKxy[0] *= coef;
                    dKxy[1] *= coef;
                    dKxy[2] *= coef;

                    tfx += dKxy[0];
                    tfy += dKxy[1];
                    tfz += dKxy[2];
                    tpo += Kxy[0] * sourcesPhysicalValues[idxSource];

                    sourcesForcesX[idxSource] -= dKxy[0];
                    sourcesForcesY[idxSource] -= dKxy[1];
                    sourcesForcesZ[idxSource] -= dKxy[2];
                    sourcesPotentials[idxSource] += Kxy[0] * tv;
                }

                targetsForcesX[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfx);
                targetsForcesY[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfy);
                targetsForcesZ[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfz);
                targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
            }
        }
    }
}

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericInner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){

    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsForcesX = inTargets->getForcesX();
    FReal*const targetsForcesY = inTargets->getForcesY();
    FReal*const targetsForcesZ = inTargets->getForcesZ();
    FReal*const targetsPotentials = inTargets->getPotentials();

    {//In this part, we compute (vectorially) the interaction
        //within the target leaf.

        const FSize nbParticlesSources = (nbParticlesTargets+NbFRealInComputeClass-1)/NbFRealInComputeClass;
        const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)targetsPhysicalValues;
        const ComputeClass*const sourcesX = (const ComputeClass*)targetsX;
        const ComputeClass*const sourcesY = (const ComputeClass*)targetsY;
        const ComputeClass*const sourcesZ = (const ComputeClass*)targetsZ;
        ComputeClass*const sourcesForcesX = (ComputeClass*)targetsForcesX;
        ComputeClass*const sourcesForcesY = (ComputeClass*)targetsForcesY;
        ComputeClass*const sourcesForcesZ = (ComputeClass*)targetsForcesZ;
        ComputeClass*const sourcesPotentials = (ComputeClass*)targetsPotentials;

        for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
            const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
            const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
            const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
            const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
            ComputeClass  tfx = FMath::Zero<ComputeClass>();
            ComputeClass  tfy = FMath::Zero<ComputeClass>();
            ComputeClass  tfz = FMath::Zero<ComputeClass>();
            ComputeClass  tpo = FMath::Zero<ComputeClass>();

            for(FSize idxSource = (idxTarget+NbFRealInComputeClass)/NbFRealInComputeClass ; idxSource < nbParticlesSources ; ++idxSource){
                ComputeClass Kxy[1];
                ComputeClass dKxy[3];
                MatrixKernel->evaluateBlockAndDerivative(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource],
                                                         tx,ty,tz,Kxy,dKxy);
                const ComputeClass coef = (tv * sourcesPhysicalValues[idxSource]);

                dKxy[0] *= coef;
                dKxy[1] *= coef;
                dKxy[2] *= coef;

                tfx += dKxy[0];
                tfy += dKxy[1];
                tfz += dKxy[2];
		tpo = FMath::FMAdd(Kxy[0],sourcesPhysicalValues[idxSource],tpo);

                sourcesForcesX[idxSource] -= dKxy[0];
                sourcesForcesY[idxSource] -= dKxy[1];
                sourcesForcesZ[idxSource] -= dKxy[2];
		sourcesPotentials[idxSource] = FMath::FMAdd(Kxy[0],tv,sourcesPotentials[idxSource]);
            }

            targetsForcesX[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfx);
            targetsForcesY[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfy);
            targetsForcesZ[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfz);
            targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
        }
    }

    for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
        const FSize limitForTarget = NbFRealInComputeClass-(idxTarget%NbFRealInComputeClass);
        for(FSize idxS = 1 ; idxS < limitForTarget ; ++idxS){
            const FSize idxSource = idxTarget + idxS;
            FReal Kxy[1];
            FReal dKxy[3];
            MatrixKernel->evaluateBlockAndDerivative(targetsX[idxSource],targetsY[idxSource],targetsZ[idxSource],
                                                     targetsX[idxTarget],targetsY[idxTarget],targetsZ[idxTarget],
                                                     Kxy,dKxy);
            const FReal coef = (targetsPhysicalValues[idxTarget] * targetsPhysicalValues[idxSource]);

            dKxy[0] *= coef;
            dKxy[1] *= coef;
            dKxy[2] *= coef;

            targetsForcesX[idxTarget] += dKxy[0];
            targetsForcesY[idxTarget] += dKxy[1];
            targetsForcesZ[idxTarget] += dKxy[2];
            targetsPotentials[idxTarget] += Kxy[0] * targetsPhysicalValues[idxSource];

            targetsForcesX[idxSource] -= dKxy[0];
            targetsForcesY[idxSource] -= dKxy[1];
            targetsForcesZ[idxSource] -= dKxy[2];
            targetsPotentials[idxSource] += Kxy[0] * targetsPhysicalValues[idxTarget];
        }
    }
}

template <class FReal, class ContainerClass, class MatrixKernelClass, class ComputeClass, int NbFRealInComputeClass>
static void GenericFullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                              const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsPhysicalValues = inTargets->getPhysicalValues();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsForcesX = inTargets->getForcesX();
    FReal*const targetsForcesY = inTargets->getForcesY();
    FReal*const targetsForcesZ = inTargets->getForcesZ();
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limiteNeighbors ; ++idxNeighbors){
        if( inNeighbors[idxNeighbors] ){
            const FSize nbParticlesSources = (inNeighbors[idxNeighbors]->getNbParticles()+NbFRealInComputeClass-1)/NbFRealInComputeClass;
            const ComputeClass*const sourcesPhysicalValues = (const ComputeClass*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const ComputeClass*const sourcesX = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[0];
            const ComputeClass*const sourcesY = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[1];
            const ComputeClass*const sourcesZ = (const ComputeClass*)inNeighbors[idxNeighbors]->getPositions()[2];

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget){
                const ComputeClass tx = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsX[idxTarget]);
                const ComputeClass ty = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsY[idxTarget]);
                const ComputeClass tz = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsZ[idxTarget]);
                const ComputeClass tv = FMath::ConvertTo<ComputeClass, const FReal*>(&targetsPhysicalValues[idxTarget]);
                ComputeClass  tfx = FMath::Zero<ComputeClass>();
                ComputeClass  tfy = FMath::Zero<ComputeClass>();
                ComputeClass  tfz = FMath::Zero<ComputeClass>();
                ComputeClass  tpo = FMath::Zero<ComputeClass>();

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource){
                    ComputeClass Kxy[1];
                    ComputeClass dKxy[3];
                    MatrixKernel->evaluateBlockAndDerivative(sourcesX[idxSource],sourcesY[idxSource],sourcesZ[idxSource],
                                                             tx,ty,tz,Kxy,dKxy);
                    const ComputeClass coef = (tv * sourcesPhysicalValues[idxSource]);

                    dKxy[0] *= coef;
                    dKxy[1] *= coef;
                    dKxy[2] *= coef;

                    tfx += dKxy[0];
                    tfy += dKxy[1];
                    tfz += dKxy[2];
                    tpo += Kxy[0] * sourcesPhysicalValues[idxSource];
                }

                targetsForcesX[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfx);
                targetsForcesY[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfy);
                targetsForcesZ[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tfz);
                targetsPotentials[idxTarget] += FMath::ConvertTo<FReal, ComputeClass>(tpo);
            }
        }
    }
}

} // End namespace

template <class FReal>
struct FP2PT{
};

#if defined(SCALFMM_USE_AVX)
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, __m256d, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }


    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, __m256d, 4>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, __m256d, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, __m256, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, __m256, 8>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, __m256, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#elif defined(SCALFMM_USE_AVX2)
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, __m512d, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, __m512d, 8>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, __m512d, 8>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, __m512, 16>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }


    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, __m512, 16>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, __m512, 16>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#elif defined(SCALFMM_USE_SSE)
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, __m128d, 2>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }


    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, __m128d, 2>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, __m128d, 2>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, __m128, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, __m128, 4>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, __m128, 4>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#else
template <>
struct FP2PT<double>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<double, ContainerClass, MatrixKernelClass, double, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<double, ContainerClass, MatrixKernelClass, double, 1>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<double, ContainerClass, MatrixKernelClass, double, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};

template <>
struct FP2PT<float>{
    template <class ContainerClass, class MatrixKernelClass>
    static void FullMutual(ContainerClass* const FRestrict inTargets, ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullMutual<float, ContainerClass, MatrixKernelClass, float, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void Inner(ContainerClass* const FRestrict inTargets, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericInner<float, ContainerClass, MatrixKernelClass, float, 1>(inTargets, MatrixKernel);
    }

    template <class ContainerClass, class MatrixKernelClass>
    static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[],
                           const int limiteNeighbors, const MatrixKernelClass *const MatrixKernel){
        FP2P::GenericFullRemote<float, ContainerClass, MatrixKernelClass, float, 1>(inTargets, inNeighbors, limiteNeighbors, MatrixKernel);
    }
};
#endif

#include "FP2PMultiRhs.hpp"

#endif // FP2P_HPP

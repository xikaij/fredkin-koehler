#ifndef BOUNDARYINTEGRALFMM_H
#define BOUNDARYINTEGRALFMM_H

#include "GeneralUserObject.h"
#include "MooseMesh.h"

class BoundaryIntegralFMM;

template<>
InputParameters validParams<BoundaryIntegralFMM>();

class BoundaryIntegralFMM : public GeneralUserObject
{
public:
  BoundaryIntegralFMM(const InputParameters & parameters);

  virtual void initialize();

  virtual void execute();

  virtual void finalize();

private:
  Real _cx, _cy, _cz, _boxWidth;
  unsigned int _TreeHeight;
};

#endif

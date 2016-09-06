#ifndef BOUNDARYINTEGRALFMM_H
#define BOUNDARYINTEGRALFMM_H

#include "GeneralUserObject.h"

class BoundaryIntaegralFMM;

template<>
InputParameters validParams<BoundaryIntegralFMM>();

class BoundaryIntegralFMM : public GeneralUserObject
{
public:
  BoundaryIntegralFMM(InputParameters parameters);

  virtual void initialize();

  virtual void execute();

  virtual void finalize();

protected:

};

#endif

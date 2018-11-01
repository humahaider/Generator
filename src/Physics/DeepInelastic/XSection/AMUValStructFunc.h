//____________________________________________________________________________
/*!

\class    genie::AMUValStructFunc

\brief    Aligarh-Valencia structure function model

\ref      <include references>

\author   Huma Haider <Huma Haider <huma.haider8 \at gmail.com>
          Aligarh Muslim University

\created  October 17, 2018

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _AMU_VALENCIA_STRUCTURE_FUNCTION_MODEL_H_
#define _AMU_VALENCIA_STRUCTURE_FUNCTION_MODEL_H_

#include "Physics/DeepInelastic/XSection/QPMDISStrucFuncBase.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/PartonDistributions/PDFModelI.h"

namespace genie {

class AMUValStructFunc : public QPMDISStrucFuncBase {

public:
  AMUValStructFunc();
  AMUValStructFunc(string config);
  virtual ~AMUValStructFunc();

  // overload QPMDISStrucFuncBase::Calculate() to override
  // the standard calculation of F1-F6
  virtual void Calculate (const Interaction * interaction) const;

  // overload Algorithm::Configure() to read the config. registry
  // at the algorithm initialization and set private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

protected:

  void Init         (void);

  // override part of the DISStructureFuncModel implementation
  // to compute all the corrections applied by the Bodek-Yang model.
  // double ScalingVar (const Interaction * i) const;
  // void   KFactors   (const Interaction * i, double & kuv,
  //                    double & kdv, double & kus, double & kds) const;

  // Bodek-Yang model-specific parameters

  //  double fA;     ///< better scaling var parameter A
  //  double fB;     ///< better scaling var parameter B
  //  double fCsU;   ///< U-sea K factor parameter
  //  double fCsD;   ///< D-sea K factor parameter
  //  double fCv1U;  ///< U-val K factor parameter
  //  double fCv2U;  ///< U-val K factor parameter
  //  double fCv1D;  ///< D-val K factor parameter
  //  double fCv2D;  ///< D-val K factor parameter
};

}         // genie namespace

#endif    // _AMU_VALENCIA_STRUCTURE_FUNCTION_MODEL_H_

//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author:  Huma Haider <Huma Haider <huma.haider8 \at gmail.com>
          Aligarh Muslim University
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/DeepInelastic/XSection/AMUValStructFunc.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
AMUValStructFunc::AMUValStructFunc() :
QPMDISStrucFuncBase("genie::AMUValStructFunc")
{
  this->Init();
}
//____________________________________________________________________________
AMUValStructFunc::AMUValStructFunc(string config):
QPMDISStrucFuncBase("genie::AMUValStructFunc", config)
{
  this->Init();
}
//____________________________________________________________________________
AMUValStructFunc::~AMUValStructFunc()
{

}
//____________________________________________________________________________
void AMUValStructFunc::Calculate(const Interaction * interaction) const
{
  // Reset mutable members
  fF1 = 0;
  fF2 = 0;
  fF3 = 0;
  fF4 = 0;
  fF5 = 0;
  fF6 = 0;

  // Get process info & perform various checks
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const InitialState & init_state = interaction->InitState();
  const Target & tgt = init_state.Tgt();

  int  nuc_pdgc    = tgt.HitNucPdg();
  int  probe_pdgc  = init_state.ProbePdg();
  bool is_p        = pdg::IsProton       ( nuc_pdgc    );
  bool is_n        = pdg::IsNeutron      ( nuc_pdgc    );
  bool is_nu       = pdg::IsNeutrino     ( probe_pdgc  );
  bool is_nubar    = pdg::IsAntiNeutrino ( probe_pdgc  );
  bool is_lepton   = pdg::IsLepton       ( probe_pdgc  );
  bool is_dm       = pdg::IsDarkMatter   ( probe_pdgc  );
  bool is_CC       = proc_info.IsWeakCC();
  bool is_NC       = proc_info.IsWeakNC();
  bool is_EM       = proc_info.IsEM();

  if ( !is_lepton && !is_dm ) return;
  if ( !is_p && !is_n       ) return;
  if ( tgt.N() == 0 && is_n ) return;
  if ( tgt.Z() == 0 && is_p ) return;

  // Flags switching on/off quark contributions so that this algorithm can be
  // used for both l + N -> l' + X, and l + q -> l' + q' level calculations

  double switch_uv    = 1.;
  double switch_us    = 1.;
  double switch_ubar  = 1.;
  double switch_dv    = 1.;
  double switch_ds    = 1.;
  double switch_dbar  = 1.;
  double switch_s     = 1.;
  double switch_sbar  = 1.;
  double switch_c     = 1.;
  double switch_cbar  = 1.;

  if(tgt.HitQrkIsSet()) {

     switch_uv    = 0.;
     switch_us    = 0.;
     switch_ubar  = 0.;
     switch_dv    = 0.;
     switch_ds    = 0.;
     switch_dbar  = 0.;
     switch_s     = 0.;
     switch_sbar  = 0.;
     switch_c     = 0.;
     switch_cbar  = 0.;

     int  qpdg = tgt.HitQrkPdg();
     bool sea  = tgt.HitSeaQrk();

     bool is_u    = pdg::IsUQuark     (qpdg);
     bool is_ubar = pdg::IsAntiUQuark (qpdg);
     bool is_d    = pdg::IsDQuark     (qpdg);
     bool is_dbar = pdg::IsAntiDQuark (qpdg);
     bool is_s    = pdg::IsSQuark     (qpdg);
     bool is_sbar = pdg::IsAntiSQuark (qpdg);
     bool is_c    = pdg::IsCQuark     (qpdg);
     bool is_cbar = pdg::IsAntiCQuark (qpdg);

     if      (!sea && is_u   ) { switch_uv   = 1; }
     else if ( sea && is_u   ) { switch_us   = 1; }
     else if ( sea && is_ubar) { switch_ubar = 1; }
     else if (!sea && is_d   ) { switch_dv   = 1; }
     else if ( sea && is_d   ) { switch_ds   = 1; }
     else if ( sea && is_dbar) { switch_dbar = 1; }
     else if ( sea && is_s   ) { switch_s    = 1; }
     else if ( sea && is_sbar) { switch_sbar = 1; }
     else if ( sea && is_c   ) { switch_c    = 1; }
     else if ( sea && is_cbar) { switch_cbar = 1; }
     else return;

     // make sure user inputs make sense
    if(is_nu    && is_CC && is_u   ) return;
    if(is_nu    && is_CC && is_c   ) return;
    if(is_nu    && is_CC && is_dbar) return;
    if(is_nu    && is_CC && is_sbar) return;
    if(is_nubar && is_CC && is_ubar) return;
    if(is_nubar && is_CC && is_cbar) return;
    if(is_nubar && is_CC && is_d   ) return;
    if(is_nubar && is_CC && is_s   ) return;
  }

  // Compute PDFs [both at (scaling-var,Q2) and (slow-rescaling-var,Q2)
  // Applying all PDF K-factors abd scaling variable corrections

  this -> CalcPDFs (interaction);

  //
  // Compute structure functions for the EM, NC and CC cases
  //

  double F2val=0, xF3val=0;

  // ***  NEUTRAL CURRENT

  // Include DM in NC
  if(is_NC) {

    if(!is_nu && !is_nubar) return;

    double GL   = (is_nu) ? ( 0.5 - (2./3.)*fSin2thw) : (     - (2./3.)*fSin2thw); // clu
    double GR   = (is_nu) ? (     - (2./3.)*fSin2thw) : ( 0.5 - (2./3.)*fSin2thw); // cru
    double GLp  = (is_nu) ? (-0.5 + (1./3.)*fSin2thw) : (       (1./3.)*fSin2thw); // cld
    double GRp  = (is_nu) ? (       (1./3.)*fSin2thw) : (-0.5 + (1./3.)*fSin2thw); // crd

    double gvu  = GL  + GR;
    double gau  = GL  - GR;
    double gvd  = GLp + GRp;
    double gad  = GLp - GRp;
    double gvu2 = TMath::Power(gvu, 2.);
    double gau2 = TMath::Power(gau, 2.);
    double gvd2 = TMath::Power(gvd, 2.);
    double gad2 = TMath::Power(gad, 2.);

    double q2   = (switch_uv   * fuv + switch_us   * fus + switch_c    * fc)  * (gvu2+gau2) +
                  (switch_dv   * fdv + switch_ds   * fds + switch_s    * fs)  * (gvd2+gad2);
    double q3   = (switch_uv   * fuv + switch_us   * fus + switch_c    * fc)  * (2*gvu*gau) +
                  (switch_dv   * fdv + switch_ds   * fds + switch_s    * fs)  * (2*gvd*gad);

    double qb2  = (switch_ubar * fus + switch_cbar * fc)  * (gvu2+gau2) +
                  (switch_dbar * fds + switch_sbar * fs)  * (gvd2+gad2);
    double qb3  = (switch_ubar * fus + switch_cbar * fc)  * (2*gvu*gau) +
                  (switch_dbar * fds + switch_sbar * fs)  * (2*gvd*gad);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("AMUValDIS", pINFO) << "f2 : q = " << q2 << ", bar{q} = " << qb2;
    LOG("AMUValDIS", pINFO) << "xf3: q = " << q3 << ", bar{q} = " << qb3;
#endif

    F2val  = q2+qb2;
    xF3val = q3-qb3;
  }

  // ***  CHARGED CURRENT

  if(is_CC) {
    double q=0, qbar=0;

    if (is_nu) {
      q    = ( switch_dv * fdv   + switch_ds * fds   ) * fVud2 +
             ( switch_s  * fs                        ) * fVus2 +
             ( switch_dv * fdv_c + switch_ds * fds_c ) * fVcd2 +
             ( switch_s  * fs_c                      ) * fVcs2;

      qbar = ( switch_ubar * fus  ) * fVud2 +
             ( switch_ubar * fus  ) * fVus2 +
             ( switch_cbar * fc_c ) * fVcd2 +
             ( switch_cbar * fc_c ) * fVcs2;
    }
    else
    if (is_nubar) {
      q    = ( switch_uv * fuv + switch_us * fus    ) * fVud2 +
             ( switch_uv * fuv + switch_us * fus    ) * fVus2 +
             ( switch_c  * fc_c                     ) * fVcd2 +
             ( switch_c  * fc_c                     ) * fVcs2;

      qbar = ( switch_dbar * fds_c ) * fVcd2 +
             ( switch_dbar * fds   ) * fVud2 +
             ( switch_sbar * fs    ) * fVus2 +
             ( switch_sbar * fs_c  ) * fVcs2;
    }
    else {
      return;
    }

    F2val  = 2*(q+qbar);
    xF3val = 2*(q-qbar);
  }

  // ***  ELECTROMAGNETIC

  if(is_EM) {

    if(!pdg::IsChargedLepton(probe_pdgc)) return;

    double sq23 = TMath::Power(2./3., 2.);
    double sq13 = TMath::Power(1./3., 2.);

    double qu   = sq23 * ( switch_uv   * fuv + switch_us * fus );
    double qd   = sq13 * ( switch_dv   * fdv + switch_ds * fds );
    double qs   = sq13 * ( switch_s    * fs  );
    double qbu  = sq23 * ( switch_ubar * fus );
    double qbd  = sq13 * ( switch_dbar * fds );
    double qbs  = sq13 * ( switch_sbar * fs  );

    double q    = qu  + qd  + qs;
    double qbar = qbu + qbd + qbs;

    F2val  = q + qbar;;
    xF3val = 0.;

  }

  double Q2val = this->Q2        (interaction);
  double x     = this->ScalingVar(interaction);
  double f     = this->NuclMod   (interaction); // nuclear modification
  double r     = this->R         (interaction); // R ~ FL

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("AMUValDIS", pDEBUG) << "Nucl. mod   = " << f;
  LOG("AMUValDIS", pDEBUG) << "R(=FL/2xF1) = " << r;
#endif

  //It was confirmed by A.Bodek that the modified scaling variable
  //should just be used to compute the strucure functions F2 and xF3,
  //but that the usual Bjorken x should be used for the relations
  //between the structure functions.
  //For the same reason remove the freezing of Q2 at 0.8 for those relations,
  //although it has not been explicitly asked to A.Bodek if it should be done.
  const Kinematics & kinematics = interaction->Kine();
  double bjx = kinematics.x();

  double a = TMath::Power(bjx,2.) / TMath::Max(Q2val, fLowQ2CutoffF1F2);
  double c = (1. + 4. * kNucleonMass2 * a) / (1.+r);

  fF3 = f * xF3val/bjx;
  fF2 = f * F2val;
  fF1 = fF2 * 0.5*c/bjx;
  fF5 = fF2/bjx;           // Albright-Jarlskog relation
  fF4 = 0.;                // Nucl.Phys.B 84, 467 (1975)

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("AMUValDIS", pDEBUG)
     << "F1-F5 = "
     << fF1 << ", " << fF2 << ", " << fF3 << ", " << fF4 << ", " << fF5;
#endif
}
//____________________________________________________________________________
void AMUValStructFunc::Configure(const Registry & config)
{
// Overload Algorithm::Configure() to read the config. registry and set
// private data members.
// QPMDISStrucFuncBase::Configure() creates the owned PDF object that gets
// configured with the specified PDFModelI

  QPMDISStrucFuncBase::Configure(config);

  // <customize as needed>
}
//____________________________________________________________________________
void AMUValStructFunc::Configure(string param_set)
{
  QPMDISStrucFuncBase::Configure(param_set);

  // <customize as needed>
}
//____________________________________________________________________________
void AMUValStructFunc::Init(void)
{
  // <initialize your variables here>
}
//____________________________________________________________________________

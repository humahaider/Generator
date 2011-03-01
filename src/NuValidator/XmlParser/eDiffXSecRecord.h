//_____________________________________________________________________________
/*!

\class    genie::nuvld::eDiffXSecRecord

\brief    XML record structure for eN differential cross section data points

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003
*/
//_____________________________________________________________________________

#ifndef _E_DIFF_XSEC_RECORD_H_
#define _E_DIFF_XSEC_RECORD_H_

#include <iostream>

#include "RecordStructure.h"

using std::ostream;

namespace genie {
namespace nuvld {
  
class eDiffXSecRecord : public RecordStructure {

public:
  static eDiffXSecRecord * Instance(void);

  friend ostream & operator<<(ostream & stream, const eDiffXSecRecord & rec_xsec); 

private:

  eDiffXSecRecord(); 
  eDiffXSecRecord(const eDiffXSecRecord & rec_xsec);
  virtual ~eDiffXSecRecord();
  
  static eDiffXSecRecord * _instance;
};

} // nuvld namespace
} // genie namespace

#endif // _E_DIFF_XSEC_RECORD_H_
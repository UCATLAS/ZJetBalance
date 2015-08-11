#ifndef ZJetBalanceMiniTreeAnaSkeleton_H
#define ZJetBalanceMiniTreeAnaSkeleton_H

#include <ZJetBalance/ZJetBalanceMiniTreeAnaBase.h>

//algorithm wrapper
#include <TH1F.h>
#include <TH2F.h>

class ZJetBalanceMiniTreeAnaSkeleton : public ZJetBalanceMiniTreeAnaBase
{
 public:
  // this is a standard constructor
  ZJetBalanceMiniTreeAnaSkeleton();
  
  virtual EL::StatusCode configure (); 
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(ZJetBalanceMiniTreeAnaSkeleton, 1);
  
 private:
  // declaration for your histograms
  TH1F* m_h_ZM; //!

};  
#endif

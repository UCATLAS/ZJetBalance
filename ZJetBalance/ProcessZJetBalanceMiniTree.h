#ifndef ZJetBalance_ProcessZJetBalanceMiniTree_H
#define ZJetBalance_ProcessZJetBalanceMiniTree_H

#include <EventLoop/StatusCode.h>
#include <EventLoop/Algorithm.h>
#include <EventLoop/Worker.h>

//algorithm wrapper
#include <xAODAnaHelpers/Algorithm.h>

// ROOT include(s):
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <vector>

class ProcessZJetBalanceMiniTree : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
 public:
  // float cutValue;
  int m_eventCounter;  //!
  
  std::string m_name;  //!
  
 private:
  bool m_debug; //! set verbose mode
  
  
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  // this is a standard constructor
  ProcessZJetBalanceMiniTree ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  
  // these are the functions not inherited from Algorithm
  virtual EL::StatusCode configure ();
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(ProcessZJetBalanceMiniTree, 1);
};

#endif

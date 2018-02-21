///////////////////////////////////////////////////////////////////////
// Author: P T Surukuchi
// Date: June 07, 2016
// Built on top of macros created by K. Gilje
// TODO:
//   Implement error level
// In the file, search for TODOs and EDITs:
//   TODO means this feature/detail has to be implemented
//   EDIT means this feature/detail/bug has to be edited
////////////////////////////////////////////////////////////////////////

#ifndef TREACTOR_HH
#define TREACTOR_HH

//Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//Include ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"
#include "TRandom3.h"

//Include Osc_Sens related headers
#include "TMacroInterface.hh"

/// Reactor class
/**
 *  Information of the reactor including dimensions.
 Use this as a base-class for deriving other reactors like HFIR.
 May need to make the data-members private and setup getters for encapsulation or const , but for now it works.
 Could create a TReactorSpectrum class separate and use it as a data-memeber here.
 */
class TReactor
{
public:
  
  /// Default constructor
  TReactor();
  
  /// Destructor
  virtual ~TReactor(){}
  
  
  ///////////////////////////////////////////////////////////////////////
  //public attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// The HFIR reactor core radius (in meters)
  double fReactorRadius = 0.2;
  
  /// The HFIR reactor core height (in meters)
  double fReactorHeight = 0.5;
  
  /// The HFIR reactor core power (MW turned into GW)
  double fReactorPower = 85.0 / 1000.0;
  
  /// Number of Events to create
  int fNEvents;
  
  /// A tree containing reactor event points
  TTree* fReactorEventTree;
  
  ///////////////////////////////////////////////////////////////////////
  //public member functions
  ///////////////////////////////////////////////////////////////////////
  /// Pull a powermap-weighted random point from the reactor
  /*
   Any subclasses have to redefine their own method, hence the method has to be declared pure virtual in case there are any derived classes
   */
  virtual TVector3 GetRandReactorPosition();
  
  /// Create a tree of events in the Reactor.
  void CreateReactorEventTree();
  
  /// This function resets all the tree variables to default in the Reactor event tree.
  void ResetReactorEventTree();
  
  /// Create Reactor events and fill a tree with values.
  void FillReactorEventTree();
  
  /// Access o reading the reactor event tree, filling the specified variables
  void ReadReactorEventTree(double&, double&, double&);
  
  /// Access to reading the Reactor event tree
  void ReadReactorEventTree(TString);
  
  /// Write the Reactor Event Tree
  void WriteReactorEventTree(TFile&);
  
  /// print information of the reactor
  virtual void Print();
  
private:
  
  ///////////////////////////////////////////////////////////////////////
  //private attributes
  ///////////////////////////////////////////////////////////////////////
  
  /// The x position of a reactor event in global detector coordinates
  double reactionGlobalX = 0.;
  
  /// The y position of a reactor event in global detector coordinates
  double reactionGlobalY = 0.;
  
  /// The z position of a reactor event in global detector coordinates
  double reactionGlobalZ = 0.;
  
  /// Check whether detector event tree is created
  bool isReactorEventTreeCreated = false;
  
  ///////////////////////////////////////////////////////////////////////
  //private member functions
  ///////////////////////////////////////////////////////////////////////
  /// Function to read inputs from macro file
  virtual void ReadMacroInputs();
  
};

#endif

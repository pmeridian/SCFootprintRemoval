// -*- C++ -*-
//
// Package:    SuperClusterFootprintRemoval
// Class:      SuperClusterFootprintRemoval
// 
/**\class SuperClusterFootprintRemoval SuperClusterFootprintRemoval.cc PFIsolation/SuperClusterFootprintRemoval/plugins/SuperClusterFootprintRemoval.cc

 Description: Implements the algo for removal of PF clusters from the SC footprint

 Implementation:
     Runs on AOD in 4_2. Electron MVA cut for 4_2.
*/
//
// Original Author:  Marco Peruzzi,32 4-C16,+41227676829,
//         Created:  Sat Sep 29 17:58:21 CEST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TGeoTube.h"
#include "TGeoPara.h"
#include "TVector3.h"
#include "TMath.h"

#ifndef __SUPERCLUSTERFOOTPRINTREMOVAL__HH__
#define __SUPERCLUSTERFOOTPRINTREMOVAL__HH__


//
// constants, enums and typedefs
//

typedef struct {
  int nxtals;
  std::vector<TVector3> xtalposition;
  std::vector<float> xtaletawidth;
  std::vector<float> xtalphiwidth;
  std::vector<std::vector<TVector3> > xtalcorners;
} sc_xtal_information;

typedef struct {
  float dR;
  float dEta;
  float dPhi;
} angular_distances_struct;


//
// class declaration
//

class SuperClusterFootprintRemoval {

public:

  SuperClusterFootprintRemoval(const edm::Event&, const edm::EventSetup&);
  ~SuperClusterFootprintRemoval();

  // get the vector of the indices of those ***photon*** PF candidates (in the pfCandidates collection) that are inside the SC footprint or are duplicata of the RECO object with SuperCluster sc
  std::vector<int> GetPFCandInFootprint(reco::SuperClusterRef sc);

  // get the angular distance between the photon pf candidate pfCandidates[pfindex] and the SuperCluster sc
  angular_distances_struct GetPFCandDistanceFromSC(reco::SuperClusterRef sc, int pfindex);

private:

  TVector3 PropagatePFCandToEcal(int pfcandindex, float position, bool isbarrel);
  sc_xtal_information GetSCXtalInfo(reco::SuperClusterRef sc);
  std::vector<int> GetMatchedPFCandidates(reco::SuperClusterRef sc);

  CaloSubdetectorGeometry *barrelGeometry;
  CaloSubdetectorGeometry *endcapGeometry;
  TGeoPara eegeom;

  edm::Handle<reco::GsfElectronCollection> electronHandle;
  edm::Handle<reco::PFCandidateCollection> pfCandidates;  

  float global_linkbyrechit_enlargement;
  
};

#endif


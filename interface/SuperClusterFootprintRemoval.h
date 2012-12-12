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
// $Id: SuperClusterFootprintRemoval.h,v 1.2 2012/11/27 13:00:45 peruzzi Exp $
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
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TGeoTube.h"
#include "TGeoPara.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRotation.h"
#include "TString.h"

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

//    The constructor can take the following arguments from configuration (default value):
//
//    isolation_cone_size_forSCremoval (0.4) : isolation cone used for PF Iso calculation
//    tag_pfCandidates_forSCremoval ("particleFlow") : collection of PF candidates to use
//    tag_Vertices_forSCremoval ("offlinePrimaryVertices") : collection of vertices to use
//    rechit_link_enlargement_forSCremoval (0.25) : enlargement of linear dimension of xtals for rechit matching

  SuperClusterFootprintRemoval(const edm::Event& iEvent, const edm::ParameterSet& iConfig, const edm::EventSetup& iSetup);
  ~SuperClusterFootprintRemoval();

  // Calculate the PF isolation around the object consisting of the SuperCluster sc. Component can be "neutral","charged" or "photon".
  // The vertexforchargediso parameter should be passed for charged isolation ONLY, and tells with respect to which vertex (within the vertices collection) we should cut on dxy (0.1 cm) and dz (0.2 cm) of the track associated to charged PF candidate. Ideally, this should be the vertex from which the object consisting of the SuperCluster sc comes from.
  float PFIsolation(TString component, reco::SuperClusterRef sc, int vertexforchargediso=-999);

  // Get the vector of the indices of those PF candidates (neutrals,charged and photons, in the pfCandidates collection) that are inside the SC footprint or are duplicata of the RECO object with SuperCluster sc
  std::vector<int> GetPFCandInFootprint(reco::SuperClusterRef sc);

  // Get the angular distance between the pf candidate pfCandidates[pfindex] and the SuperCluster sc
  angular_distances_struct GetPFCandHitDistanceFromSC(reco::SuperClusterRef sc, int pfindex);

private:

  TVector3 PropagatePFCandToEcal(int pfcandindex, float position, bool isbarrel);
  sc_xtal_information GetSCXtalInfo(reco::SuperClusterRef sc);
  std::vector<int> GetMatchedPFCandidates(reco::SuperClusterRef sc);

  CaloSubdetectorGeometry *barrelGeometry;
  CaloSubdetectorGeometry *endcapGeometry;
  TGeoPara eegeom;
  MagneticField *magField;

  edm::Handle<reco::PFCandidateCollection> pfCandidates;  
  edm::Handle<reco::VertexCollection> vertexHandle;

  double global_linkbyrechit_enlargement;
  double global_isolation_cone_size;

  int FindPFCandType(int id);
  
};

#endif


// -*- C++ -*-
//
// Package:    SuperClusterFootprintRemoval
// Class:      SuperClusterFootprintRemoval
// 
/**\class SuperClusterFootprintRemoval SuperClusterFootprintRemoval.cc PFIsolation/SuperClusterFootprintRemoval/plugins/SuperClusterFootprintRemoval.cc

 Description: Implements the algo for removal of PF clusters from the SC footprint

 Implementation:
     Runs on AOD.
*/
//
// Original Author:  Marco Peruzzi,32 4-C16,+41227676829,
//         Created:  Sat Sep 29 17:58:21 CEST 2012
// $Id: SuperClusterFootprintRemoval.cc,v 1.6 2012/12/15 18:10:07 peruzzi Exp $
//
//


#ifndef __SUPERCLUSTERFOOTPRINTREMOVAL__CC__
#define __SUPERCLUSTERFOOTPRINTREMOVAL__CC__

#include "../interface/SuperClusterFootprintRemoval.h"

//
// constructors and destructor
//
SuperClusterFootprintRemoval::SuperClusterFootprintRemoval(const edm::Event& iEvent, const edm::ParameterSet& iConfig, const edm::EventSetup& iSetup)
{
   //now do what ever initialization is needed

  eegeom = TGeoPara(1,1,1,0,0,0);

  global_isolation_cone_size = iConfig.getUntrackedParameter<double>("isolation_cone_size_forSCremoval",0.4);
  global_linkbyrechit_enlargement = iConfig.getUntrackedParameter<double>("rechit_link_enlargement_forSCremoval",0.25);

  edm::ESHandle<CaloGeometry> geometry ;
  iSetup.get<CaloGeometryRecord>().get(geometry);
  barrelGeometry = (CaloSubdetectorGeometry*)(geometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel));
  endcapGeometry = (CaloSubdetectorGeometry*)(geometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap));

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  magField = (MagneticField*)(magneticField.product());

  //PFcandidates
  edm::InputTag pfCandidateTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_pfCandidates_forSCremoval",edm::InputTag("particleFlow"));
  iEvent.getByLabel(pfCandidateTag, pfCandidates);

  //Vertices
  edm::InputTag fVertexTag = iConfig.getUntrackedParameter<edm::InputTag>("tag_Vertices_forSCremoval",edm::InputTag("offlinePrimaryVertices"));
  iEvent.getByLabel(fVertexTag, vertexHandle);

}


SuperClusterFootprintRemoval::~SuperClusterFootprintRemoval()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


TVector3 SuperClusterFootprintRemoval::PropagatePFCandToEcal(int pfcandindex, float position, bool isbarrel){
  // WARNING: this propagates until EE+ or EE- at the given TMath::Abs(position.z()) for isbarrel=0, depending on where the candidate is pointing.



  int i = pfcandindex;
  int type = FindPFCandType((*pfCandidates)[i].pdgId());

  if (type>2) {
    std::cout << "Asking propagation for lepton, not implemented" << std::endl;
    return TVector3(0,0,1e10);
  }

  TVector3 pfvertex((*pfCandidates)[i].vx(),(*pfCandidates)[i].vy(),(*pfCandidates)[i].vz());
  TVector3 pfmomentum((*pfCandidates)[i].px(),(*pfCandidates)[i].py(),(*pfCandidates)[i].pz());
  pfmomentum = pfmomentum.Unit();
  TVector3 ecalpfhit(0,0,0);
  if (isbarrel){
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    double dist = TGeoTube::DistFromInsideS(p,d,0,position,1e+10);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }
  else { // EE
    double dim[6]={1e+10,1e+10,fabs(position),0,0,0};
    double p[3] = {pfvertex.x(),pfvertex.y(),pfvertex.z()};
    double d[3] = {pfmomentum.x(),pfmomentum.y(),pfmomentum.z()};
    eegeom.SetDimensions(dim);
    double dist = eegeom.DistFromInside(p,d);
    ecalpfhit = pfvertex + dist*pfmomentum;
  }


  if (type==1){ // APPROXIMATE helix propagation

    float field = magField->inTesla(GlobalPoint(0.,0.,0.)).z();
    float charge = (*pfCandidates)[i].charge();
    float curvature = ((*pfCandidates)[i].pt())/(0.3*field*charge)*100; // curvature radius in cm
    float final_radius = ecalpfhit.Perp();
    float addphi = -TMath::ASin(final_radius/curvature/2);

    TRotation r;
    r.RotateZ(addphi);
    ecalpfhit *= r;

  }


  return ecalpfhit;

}

sc_xtal_information SuperClusterFootprintRemoval::GetSCXtalInfo(reco::SuperClusterRef sc){

  sc_xtal_information out;

  bool isbarrel = (sc->seed()->seed().subdetId()==EcalBarrel);

  std::vector<DetId> cristalli;         
  for (reco::CaloCluster_iterator bc=sc->clustersBegin(); bc!=sc->clustersEnd(); ++bc){
    const std::vector< std::pair<DetId, float> > & seedrechits = (*bc)->hitsAndFractions();
    for (unsigned int i=0; i<seedrechits.size(); i++) cristalli.push_back(seedrechits[i].first);
    sort(cristalli.begin(),cristalli.end());
    std::vector<DetId>::iterator it;
    it = unique(cristalli.begin(),cristalli.end());
    cristalli.resize(it-cristalli.begin());
  }
          
  unsigned int i=0;
  for (i=0; i<cristalli.size(); i++){

    if (cristalli.at(i).subdetId()!=EcalBarrel && cristalli.at(i).subdetId()!=EcalEndcap) continue;

    CaloCellGeometry *cellGeometry = NULL;
    if (cristalli.at(i).subdetId()!=sc->seed()->seed().subdetId()) {
      std::cout << "Problem with subdetId() in SuperClusterFootprintRemoval" << std::endl;
      continue;
    } 

    EBDetId ebDetId;
    EEDetId eeDetId;
    if (isbarrel) {
      ebDetId = cristalli.at(i); 
      cellGeometry = (CaloCellGeometry*)(barrelGeometry->getGeometry(ebDetId));
    }
    else {
      eeDetId = cristalli.at(i);
      cellGeometry = (CaloCellGeometry*)(endcapGeometry->getGeometry(eeDetId));
    }
    TVector3 xtal_position(cellGeometry->getPosition().x(),cellGeometry->getPosition().y(),cellGeometry->getPosition().z());

    float dphi; float deta;
    if (isbarrel){
    dphi=(dynamic_cast<const EcalBarrelGeometry*>(barrelGeometry))->deltaPhi(ebDetId);
    deta=(dynamic_cast<const EcalBarrelGeometry*>(barrelGeometry))->deltaEta(ebDetId);
    }
    else {
    dphi=(dynamic_cast<const EcalEndcapGeometry*>(endcapGeometry))->deltaPhi(eeDetId);
    deta=(dynamic_cast<const EcalEndcapGeometry*>(endcapGeometry))->deltaEta(eeDetId);
    }
    
    const CaloCellGeometry::CornersVec& cellCorners (cellGeometry->getCorners());
    std::vector<TVector3> corners;
    for (int k=0; k<4; k++){
      TVector3 thiscorner((float)(cellCorners[k].x()),(float)(cellCorners[k].y()),(float)(cellCorners[k].z()));
      corners.push_back(thiscorner);
    }
    out.xtalposition.push_back(xtal_position);
    out.xtaletawidth.push_back(deta);
    out.xtalphiwidth.push_back(dphi);
    out.xtalcorners.push_back(corners);
  }
  out.nxtals=i;

  assert (out.nxtals==(int)out.xtalposition.size());

  return out;

}

std::vector<int> SuperClusterFootprintRemoval::GetMatchedPFCandidates(reco::SuperClusterRef sc){

  std::vector<int> out;

  for (unsigned int i=0; i<pfCandidates->size(); i++) {
    if ((*pfCandidates)[i].pdgId()==22) {
      if ((*pfCandidates)[i].mva_nothing_gamma()>0){
	if( (*pfCandidates)[i].superClusterRef()==sc) {
	  out.push_back(i);
	}
      }
    }
  }

  return out;

}

std::vector<int> SuperClusterFootprintRemoval::GetPFCandInFootprint(reco::SuperClusterRef sc){

  bool isbarrel = (sc->seed()->seed().subdetId()==EcalBarrel);

  sc_xtal_information infos = GetSCXtalInfo(sc);
  std::vector<int> matchedpfcand = GetMatchedPFCandidates(sc);

  std::vector<int> result;

  for (unsigned int i=0; i<pfCandidates->size(); i++){

    int type = FindPFCandType((*pfCandidates)[i].pdgId());
    if (type>2) continue;

    bool inside=false;

    for (unsigned int j=0; j<matchedpfcand.size(); j++) if ((int)i==matchedpfcand.at(j)) inside=true;

    for (int j=0; j<infos.nxtals; j++){
      
      TVector3 xtal_position = infos.xtalposition[j];
      TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? xtal_position.Perp() : xtal_position.z(), isbarrel);

      if (isbarrel){
	float xtalEtaWidth = infos.xtaletawidth[j]*(1.0+global_linkbyrechit_enlargement);
	float xtalPhiWidth = infos.xtalphiwidth[j]*(1.0+global_linkbyrechit_enlargement);
	float deta = ecalpfhit.Eta()-xtal_position.Eta();
	float dphi = reco::deltaPhi(ecalpfhit.Phi(),xtal_position.Phi());
	if (fabs(deta)<xtalEtaWidth/2 && fabs(dphi)<xtalPhiWidth/2) inside=true;
      }
      else { // EE
	if (ecalpfhit.z()*xtal_position.z()>0){
	  TVector3 xtal_corners[4];
	  for (int k=0; k<4; k++) xtal_corners[k] = infos.xtalcorners[j][k];
	  float hitx = ecalpfhit.x();
	  float hity = ecalpfhit.y();
	  float polx[5];
	  float poly[5];
	  for (int k=0; k<4; k++) polx[k] = xtal_corners[k].x();
	  for (int k=0; k<4; k++) poly[k] = xtal_corners[k].y();
	  polx[4]=polx[0]; poly[4]=poly[0]; // closed polygon
	  float centerx = (polx[0]+polx[1]+polx[2]+polx[3])/4;
	  float centery = (poly[0]+poly[1]+poly[2]+poly[3])/4;
	  hitx = centerx + (hitx-centerx)/(1.0+global_linkbyrechit_enlargement);
	  hity = centery + (hity-centery)/(1.0+global_linkbyrechit_enlargement);
	  if (TMath::IsInside(hitx,hity,5,polx,poly)) inside=true;
	}
      }

    }

    if (inside) result.push_back(i);

  }

  return result;

}


angular_distances_struct SuperClusterFootprintRemoval::GetPFCandHitDistanceFromSC(reco::SuperClusterRef sc, int pfindex){

  int i = pfindex;
  int type = FindPFCandType((*pfCandidates)[i].pdgId());

  if (type>2) {
    std::cout << "propagation not implemented for lepton objects!!!" << std::endl;
    angular_distances_struct out;
    out.dR=999; out.dEta=999; out.dPhi=999;
    return out;
  }

  bool isbarrel = (sc->seed()->seed().subdetId()==EcalBarrel);

  TVector3 sc_position = TVector3(sc->x(),sc->y(),sc->z());

  TVector3 pfvertex((*pfCandidates)[i].vx(),(*pfCandidates)[i].vy(),(*pfCandidates)[i].vz());
  TVector3 pfmomentum((*pfCandidates)[i].px(),(*pfCandidates)[i].py(),(*pfCandidates)[i].pz());
  pfmomentum = pfmomentum.Unit();

  TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? sc_position.Perp() : sc_position.z(),isbarrel);

  angular_distances_struct out;
  out.dR = reco::deltaR(sc_position.Eta(),sc_position.Phi(),ecalpfhit.Eta(),ecalpfhit.Phi());
  out.dEta = ecalpfhit.Eta()-sc_position.Eta();
  out.dPhi = reco::deltaPhi(ecalpfhit.Phi(),sc_position.Phi());

  return out;

}

int SuperClusterFootprintRemoval::FindPFCandType(int id){

  int type = -1;

  if (id==111 || id==130 || id==310 || id==2112) type=0; //neutral hadrons
  if (fabs(id)==211 || fabs(id)==321 || id==999211 || fabs(id)==2212) type=1; //charged hadrons
  if (id==22) type=2; //photons
  if (fabs(id)==11) type=3; //electrons
  if (fabs(id)==13) type=4; //muons

  return type;
}

float SuperClusterFootprintRemoval::PFIsolation(TString component, reco::SuperClusterRef sc, int vertexforchargediso){

  if (component!="charged") {
    if (vertexforchargediso!=-999) std::cout << "WARNING: Why are you specifying a vertex for neutral or photon isolation? This will be ignored." << std::endl;
    vertexforchargediso=-1;
  }

  if (component=="charged") {
      if (vertexforchargediso==-999) {std::cout << "WARNING: You did not specify a vertex for the charged component of PF isolation. Deactivating vertex cuts (vertexforchargediso=-1) by default." << std::endl; vertexforchargediso=-1;}
      if (vertexforchargediso > (int)(vertexHandle->size())-1 || vertexforchargediso<-1) {std::cout << "ERROR: Invalid vertexforchargediso specified. Returning 999." << std::endl; return 999;}
    }

  float result = 0;
  std::vector<int> removed = GetPFCandInFootprint(sc);

  int thistype=-999;
  if (component=="neutral") thistype=0;
  if (component=="charged") thistype=1;
  if (component=="photon")  thistype=2;

  if (thistype==-999) {std::cout << "ERROR: Incorrect PF isolation component selected. Returning 999." << std::endl; return 999;}

  for (unsigned int i=0; i<pfCandidates->size(); i++){

    int type = FindPFCandType((*pfCandidates)[i].pdgId());
    if (type!=thistype) continue;

    bool toremove = false;
    for (unsigned int j=0; j<removed.size(); j++) if ((int)i==removed.at(j)) toremove = true;
    if (toremove) continue;

    angular_distances_struct distance = GetPFCandHitDistanceFromSC(sc,i);
    if (distance.dR>global_isolation_cone_size) continue;

    if (type==1 && vertexforchargediso>-1){
      TVector3 pfvertex((*pfCandidates)[i].vx(),(*pfCandidates)[i].vy(),(*pfCandidates)[i].vz());
      TVector3 vtxmom((*pfCandidates)[i].trackRef()->px(),(*pfCandidates)[i].trackRef()->py(),(*pfCandidates)[i].trackRef()->pz());
      TVector3 phovtx((*vertexHandle)[vertexforchargediso].x(),(*vertexHandle)[vertexforchargediso].y(),(*vertexHandle)[vertexforchargediso].z());
      float dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
      float dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
      dxy=fabs(dxy);
      dz=fabs(dz);
      if (dz>0.2) continue;
      if (dxy>0.1) continue;
    }

    result+=(*pfCandidates)[i].pt();

  }

  return result;

}

#endif

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
// $Id: SuperClusterFootprintRemoval.cc,v 1.10 2013/05/15 12:57:16 peruzzi Exp $
//
//


#ifndef __SUPERCLUSTERFOOTPRINTREMOVAL__CC__
#define __SUPERCLUSTERFOOTPRINTREMOVAL__CC__

#include "../interface/SuperClusterFootprintRemoval.h"

//
// constructors and destructor
//
SuperClusterFootprintRemoval::SuperClusterFootprintRemoval(const edm::Event& iEvent, const edm::EventSetup& iSetup, double conesize, edm::ParameterSet iConfig)
{
   //now do what ever initialization is needed

  eegeom = TGeoPara(1,1,1,0,0,0);

  global_isolation_cone_size = conesize;
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

  //Photons
  iEvent.getByLabel("photons", photonHandle);

  //Jets
  iEvent.getByLabel("ak5PFJets", jetHandle);

  randomgen = new TRandom3(0);

}


SuperClusterFootprintRemoval::~SuperClusterFootprintRemoval()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete randomgen;

}


//
// member functions
//


TVector3 SuperClusterFootprintRemoval::PropagatePFCandToEcal(int pfcandindex, float position, bool isbarrel){
  // WARNING: this propagates until EE+ or EE- at the given TMath::Abs(position.z()) for isbarrel=0, depending on where the candidate is pointing.

  if (!((*pfCandidates)[pfcandindex].pt()>0)) {
    std::cout << "Warning: called propagation to ECAL for object with negative or zero pt. Returning TVector3(0,0,1e10)." << std::endl;
    return TVector3(0,0,1e10);
  }

  int i = pfcandindex;
  int type = FindPFCandType((*pfCandidates)[i].pdgId());

  if (type>2) {
    std::cout << "Asking propagation for lepton, not implemented. Returning TVector3(0,0,1e10)." << std::endl;
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
    bool isbarrel = (cristalli.at(i).subdetId()==EcalBarrel);

    CaloCellGeometry *cellGeometry = NULL;
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

std::vector<int> SuperClusterFootprintRemoval::GetPFCandInFootprint(reco::SuperClusterRef sc, float rotation_phi){

  bool isbarrel = (fabs(sc->eta())<1.5);

  sc_xtal_information infos = GetSCXtalInfo(sc);
  std::vector<int> matchedpfcand = GetMatchedPFCandidates(sc);

  std::vector<int> result;

  for (unsigned int i=0; i<pfCandidates->size(); i++){

    if (!((*pfCandidates)[i].pt()>0)) continue;

    int type = FindPFCandType((*pfCandidates)[i].pdgId());
    if (type>2) continue;

    bool inside=false;

    for (unsigned int j=0; j<matchedpfcand.size(); j++) if ((int)i==matchedpfcand.at(j)) inside=true;

    for (int j=0; j<infos.nxtals; j++){
      
      TVector3 xtal_position = infos.xtalposition[j];
      if (rotation_phi!=0) {
	TRotation r; r.RotateZ(rotation_phi);
	xtal_position *= r;
      }

      TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? xtal_position.Perp() : xtal_position.z(), isbarrel);

      if (ecalpfhit.Perp()==0) continue;

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


angular_distances_struct SuperClusterFootprintRemoval::GetPFCandHitDistanceFromSC(reco::SuperClusterRef sc, int pfindex, float rotation_phi){

  int i = pfindex;
  int type = FindPFCandType((*pfCandidates)[i].pdgId());

  if (type>2) {
    std::cout << "propagation not implemented for lepton objects!!!" << std::endl;
    angular_distances_struct out;
    out.dR=999; out.dEta=999; out.dPhi=999;
    return out;
  }

  bool isbarrel = (fabs(sc->eta())<1.5);

  TVector3 sc_position = TVector3(sc->x(),sc->y(),sc->z());
  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    sc_position *= r;
  }

  TVector3 pfvertex((*pfCandidates)[i].vx(),(*pfCandidates)[i].vy(),(*pfCandidates)[i].vz());
  TVector3 pfmomentum((*pfCandidates)[i].px(),(*pfCandidates)[i].py(),(*pfCandidates)[i].pz());
  pfmomentum = pfmomentum.Unit();

  TVector3 ecalpfhit = PropagatePFCandToEcal(i,isbarrel ? sc_position.Perp() : sc_position.z(),isbarrel);

  if (ecalpfhit.Perp()==0){
    std::cout << "GetPFCandHitDistanceFromSC: impact position found in the origin of the transverse plane. Returning error state." << std::endl;
    angular_distances_struct out;
    out.dR = 999;
    out.dEta = 999;
    out.dPhi = 999;
    return out;
  }

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

PFIsolation_struct SuperClusterFootprintRemoval::PFIsolation(reco::SuperClusterRef sc, int vertexforchargediso, float rotation_phi)
{
  PFIsolation_struct out;
  out.chargediso=999;
  out.chargediso_primvtx=999;
  out.neutraliso=999;
  out.photoniso=999;
  
  if (vertexforchargediso==-999) {std::cout << "WARNING: You did not specify a vertex for the charged component of PF isolation. Deactivating vertex cuts (vertexforchargediso=-1) by default." << std::endl; vertexforchargediso=-1;}
  if (vertexforchargediso > (int)(vertexHandle->size())-1 || vertexforchargediso<-1) {std::cout << "ERROR: Invalid vertexforchargediso specified. Returning 999." << std::endl; return out;}
  
  
  std::vector<int> removed = GetPFCandInFootprint(sc,rotation_phi);
  out.chargediso=0;
  out.chargediso_primvtx=0;
  out.neutraliso=0;
  out.photoniso=0;

  for (unsigned int i=0; i<pfCandidates->size(); i++){

    if (!((*pfCandidates)[i].pt()>0)) continue;

    int type = FindPFCandType((*pfCandidates)[i].pdgId());
    if (type>2 ) continue;

    angular_distances_struct distance = GetPFCandHitDistanceFromSC(sc,i,rotation_phi);
    if (distance.dR>global_isolation_cone_size) continue;

    bool toremove = false;
    for (unsigned int j=0; j<removed.size(); j++) if ((int)i==removed.at(j)) toremove = true;
    if (toremove) continue;

    out.chargediso=0;
    out.chargediso_primvtx=999;
    out.neutraliso=0;
    out.photoniso=0;
    if (type==1 ){
	  TVector3 pfvertex((*pfCandidates)[i].vx(),(*pfCandidates)[i].vy(),(*pfCandidates)[i].vz());
	  TVector3 vtxmom((*pfCandidates)[i].trackRef()->px(),(*pfCandidates)[i].trackRef()->py(),(*pfCandidates)[i].trackRef()->pz());
	  TVector3 phovtx((*vertexHandle)[0].x(),(*vertexHandle)[0].y(),(*vertexHandle)[0].z());
	  float dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
	  float dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
	  dxy=fabs(dxy);
	  dz=fabs(dz);
	  
	  if (dz>0.2) continue;
	  if (dxy>0.1) continue;
	  out.chargediso_primvtx+=(*pfCandidates)[i].pt();

	  if (vertexforchargediso>-1)
	    {
	      
	      TVector3 phovtx((*vertexHandle)[vertexforchargediso].x(),(*vertexHandle)[vertexforchargediso].y(),(*vertexHandle)[vertexforchargediso].z());
	      float dxy = ( -(pfvertex.x()-phovtx.x())*vtxmom.y() +(pfvertex.y()-phovtx.y())*vtxmom.x() ) / vtxmom.Perp();
	      float dz = (pfvertex.z()-phovtx.z()) - ( (pfvertex.x()-phovtx.x())*vtxmom.x() + (pfvertex.y()-phovtx.y())*vtxmom.y() ) / vtxmom.Perp() * vtxmom.z() / vtxmom.Perp();
	      dxy=fabs(dxy);
	      dz=fabs(dz);
	      
	      if (dz>0.2) continue;
	      if (dxy>0.1) continue;
	      out.chargediso+=(*pfCandidates)[i].pt();
	    }
    }
    else if (type==0 ){
      out.neutraliso=0;
      out.neutraliso+=(*pfCandidates)[i].pt();
    }
    else if (type==2 ){
      out.photoniso=0;
      out.photoniso+=(*pfCandidates)[i].pt();
    }

  }
  return out;
}

bool SuperClusterFootprintRemoval::FindCloseJetsAndPhotons(reco::SuperClusterRef sc, float rotation_phi){

  TVector3 photon_position = TVector3(sc->x(),sc->y(),sc->z());
  if (rotation_phi!=0) {
    TRotation r; r.RotateZ(rotation_phi);
    photon_position *= r;
  }
  double eta = photon_position.Eta();
  double phi = photon_position.Phi();
  
  const float mindR = 0.8;
  bool found=false;

  for (reco::PFJetCollection::const_iterator jet=jetHandle->begin(); jet!=jetHandle->end(); jet++){
    if (jet->pt()<20) continue;
    float dR = reco::deltaR(eta,phi,jet->eta(),jet->phi());
    if (dR<mindR) found=true;
  }

  for (reco::PhotonCollection::const_iterator pho=photonHandle->begin(); pho!=photonHandle->end(); pho++){
    if (pho->pt()<10) continue;
    float dR = reco::deltaR(eta,phi,pho->eta(),pho->phi());
    if (dR<mindR) found=true;
  }

  return found;

}

PFIsolation_RandomCone_struct SuperClusterFootprintRemoval::RandomConeIsolation(reco::SuperClusterRef sc, int vertexforchargediso){

  PFIsolation_RandomCone_struct out;
  out.iso.chargediso=999;
  out.iso.chargediso_primvtx=999;
  out.iso.neutraliso=999;
  out.iso.photoniso=999;
  out.randomCone.randomcone_eta=999;
  out.randomCone.randomcone_phi=999;
  out.randomCone.randomcone_isok=true;

  const double pi = TMath::Pi();

  double rotation_phi = pi/2;

  bool isok = !(FindCloseJetsAndPhotons(sc,rotation_phi));
  if (!isok) {
    rotation_phi = -pi/2;
    isok=!(FindCloseJetsAndPhotons(sc,rotation_phi));
  }

  int count=0;
  while (!isok && count<20) {
    rotation_phi = randomgen->Uniform(0.8,2*pi-0.8);
    isok=!(FindCloseJetsAndPhotons(sc,rotation_phi));
    count++;
  }

  if (count==20){
    //    std::cout << "It was not possible to find a suitable direction for the random cone in this event. This is not a problem."  << std::endl;
    out.randomCone.randomcone_isok=false;
    return out;
  };

  
  out.iso=PFIsolation(sc,vertexforchargediso,rotation_phi);

  out.randomCone.randomcone_eta=TVector3(sc->x(),sc->y(),sc->z()).Eta();
  float newphi = TVector3(sc->x(),sc->y(),sc->z()).Phi()+rotation_phi;
  while (newphi>pi) newphi-=2*pi;
  while (newphi<-pi) newphi+=2*pi;
  out.randomCone.randomcone_phi=newphi;

  return out;  
}

#endif

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "fastjet/Selector.hh" //.......... Background Sutraction event by event
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"//.......... Background Sutraction event by event
//#include "include/tools/Subtractor.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"

#include "TPDGCode.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "THn.h"

#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

static const int debug = 0;
static const int do_bkg = 0;
static const int charged_jets = 1;
static const int charged_constituents = 1;
static const float pt_cut_rho = 1; // GeV, like CMS; only on particles for rho
static const float max_r_rho = 0.5; // has to be at least equal to largest jet radius


// cuts which i don't think are necessary in this model
static const float LH_CUT = 0.;
static const float CO_CUT = 0.;

int is_stable(const HepMC::GenParticle *part) {
    // copied from AliStack::IsStable()	 
    int pdg = abs(part->pdg_id());
    if(pdg>1000000000)return kTRUE;

    const Int_t kNstable = 18;
    Int_t i;


    Int_t pdgStable[kNstable] = {
        kGamma,             // Photon
        kElectron,          // Electron
        kMuonPlus,          // Muon 
        kPiPlus,            // Pion
        kKPlus,             // Kaon
        kK0Short,           // K0s
        kK0Long,            // K0l
        kProton,            // Proton 
        kNeutron,           // Neutron
        kLambda0,           // Lambda_0
        kSigmaMinus,        // Sigma Minus
        kSigmaPlus,         // Sigma Plus
        3312,               // Xsi Minus 
        3322,               // Xsi 
        3334,               // Omega
        kNuE,               // Electron Neutrino 
        kNuMu,              // Muon Neutrino
        kNuTau              // Tau Neutrino
    };

    Bool_t isStable = kFALSE;
    for (i = 0; i < kNstable; i++) {
        if (pdg == abs(pdgStable[i])) {
            isStable = kTRUE;
            break;
        }
    }

    return isStable;
}

int is_stable_charged(const HepMC::GenParticle *part) {
    // copied from AliStack::IsStable()	 
    int pdg = abs(part->pdg_id());
    if(pdg>1000000000)return kTRUE;

    const Int_t kNstableCharged = 9;
    Int_t i;


    Int_t pdgStableCharged[kNstableCharged] = {
        kElectron,          // Electron
        kMuonPlus,          // Muon 
        kPiPlus,            // Pion
        kKPlus,             // Kaon
        kProton,            // Proton 
        kSigmaMinus,        // Sigma Minus
        kSigmaPlus,         // Sigma Plus
        3312,               // Xsi Minus 
        3334                // Omega
    };

    Bool_t isStable = kFALSE;
    for (i = 0; i < kNstableCharged; i++) {
        if (pdg == abs(pdgStableCharged[i])) {
            isStable = kTRUE;
            break;
        }
    }

    return isStable;
}

int is_charged(const HepMC::GenParticle *part) {
    int abs_kf = abs(part->pdg_id());

    if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
        return 1;
    else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16)
        cout << " Unexpected particle: kf=" << abs_kf << endl;
    return 0;
}

float dphi(float phi1, float phi2) {
    float dphi=phi1-phi2;
    float pi = 3.14159;
    if (dphi < 0)
        dphi+=2*pi;
    if (dphi > 2.*pi)
        dphi-=2*pi;
    /*
       if (dphi < -TMath::Pi())
       dphi+=2*TMath::Pi();
       if (dphi > TMath::Pi())
       dphi-=2*TMath::Pi();
       */
    return dphi;
}


int main(int argc, char **argv) {
    //
    // Takes two arguments: infile (HEPMC format) outfile (ROOT format)
    //

    if (argc < 2) {
        cerr << "Need two arguments: infile outfile" << endl << "infile is HEPMC ascii format; outfile will be root format" << endl;
        return 1;
    }

    char *inname = argv[1];
    // specify an input file
    HepMC::IO_GenEvent ascii_in(inname,std::ios::in);

    // Make histos

    char *outname = argv[2];
    cout << "Input: " << inname << ", output " << outname << endl;

    TFile fout(outname,"RECREATE");
    static const int nR = 1; // radii 0.2, 0.3, 0.4, 0.5

    TH2F *hPtJetEta[nR];

    TH2F *hPtJetEtaBgSub[nR];//.......... Background Sutraction event by event

    TH2F *hPtLeadJetEta[nR];
    TH2F *hRhoNTrackLeadJet[nR];
    TH2F *hRhoPtLeadJet[nR];
    TH2F *hRhoNTrackLeadJetBG[nR];
    TH2F *hRhoPtLeadJetBG[nR];
    TH2F *hPtLeadJetPtTrack[nR];
    TH2F *hPtLeadJetZ[nR];
    TH2F *hPtLeadJetXi[nR];
    THnF *hPtEtaLeadJetDetaDphiP[nR];
    THnF *hPtEtaLeadJetDetaDphiN[nR];  

    TH2F *hPtLeadJetRho[nR];

    const int nEtaBins = 120;
    const float maxEtaHist = 6;

    TH2D* dPhiTracks = new TH2D("dPhiTracks", "dPhiTracks", 100, -1.*TMath::Pi(), TMath::Pi(), 100, 0, 100);
    TH2D* dPhiJets = new TH2D("dPhiJets", "dPhiJets", 8, 0 , TMath::TwoPi(), 150, 0, 150);
    dPhiJets->Sumw2();
    dPhiTracks->Sumw2();


    TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
    hNEvent->Sumw2();
    TH1F *hPtLead = new TH1F("hPtLead","leading hadron pt;p_{T}",100,0,100);
    hPtLead->Sumw2();
    TH1D *hPtAllTrack = new TH1D("hPtAllTrack","Distribution of projected Pt;p_{T}",600,-150,150);
    hPtAllTrack->Sumw2();
    TH2D *hPtPartEta = new TH2D("hPtPartEta","particle pt,eta;p_{T};#eta",100,0,100,nEtaBins,-maxEtaHist,maxEtaHist);
    hPtPartEta->Sumw2();
    TH2D *hPtChTrackEta = new TH2D("hPtChTrackEta","charged particle pt,eta;p_{T};#eta",100,0,100,nEtaBins,-maxEtaHist,maxEtaHist);
    hPtChTrackEta->Sumw2();

    TH3F *hPtAllTrackLeadJetEta = new TH3F("hPtAllTrackLeadJetEta","Pt AllTrack vs LeadJet vs Eta;p_{T,Lead_Jet}(GeV/c);p_{T,All_Track}(GeV/c);#eta",150,0,150,300,-150,150,nEtaBins,-maxEtaHist,maxEtaHist);
    hPtAllTrackLeadJetEta->Sumw2();

    for (Int_t iR = 0; iR < nR; iR++)
    {
        Float_t jetR = 0.2 + 0.1*iR;
        char hname[256];
        char htitle[256];
        sprintf(hname,"hPtJetEta_R%02d",iR+2);
        sprintf(htitle,"jet spectrum R=%.1f;p_{T,Jet} (GeV/c);#eta",jetR);
        hPtJetEta[iR] = new TH2F(hname,htitle,150,0,150,nEtaBins,-maxEtaHist,maxEtaHist);
        hPtJetEta[iR]->Sumw2();

        sprintf(hname,"hPtJetEtaBgSub_R%02d",iR+2);
        sprintf(htitle,"jet spectrum Background Subtraction R=%.1f;p_{T,Jet} (GeV/c);#eta",jetR);
        hPtJetEtaBgSub[iR] = new TH2F(hname,htitle,150,0,150,nEtaBins,-maxEtaHist,maxEtaHist);
        hPtJetEtaBgSub[iR]->Sumw2();

        sprintf(hname,"hPtLeadJetEta_R%02d",iR+2);
        sprintf(htitle,"lead_jet spectrum R=%.1f;p_{T,Lead_jet} (GeV/c);#eta",jetR);
        hPtLeadJetEta[iR] = new TH2F(hname,htitle,150,0,150,nEtaBins,-maxEtaHist,maxEtaHist);
        hPtLeadJetEta[iR]->Sumw2();

        sprintf(hname,"hPtLeadJetPtTrack_R%02d",iR+2);
        sprintf(htitle,"track p_{T} in jets R=%.1f;p_{T,Lead_jet} (GeV/c);p_{T,track}",jetR);
        hPtLeadJetPtTrack[iR] = new TH2F(hname,htitle,150,0,150,150,0,150);
        hPtLeadJetPtTrack[iR]->Sumw2();

        sprintf(hname,"hPtLeadJetZ_R%02d",iR+2);
        sprintf(htitle,"track Z in jets R=%.1f;p_{T,Lead_jet} (GeV/c);z",jetR);
        hPtLeadJetZ[iR] = new TH2F(hname,htitle,150,0,150,50,0,1);
        hPtLeadJetZ[iR]->Sumw2();
        sprintf(hname,"hPtLeadJetXi_R%02d",iR+2);
        sprintf(htitle,"track #xi in jets R=%.1f;p_{T,Lead_jet} (GeV/c);#xi",jetR);
        hPtLeadJetXi[iR] = new TH2F(hname,htitle,150,0,150,70,0,7);
        hPtLeadJetXi[iR]->Sumw2();

        sprintf(hname,"hRhoNTrackLeadJet_R%02d",iR+2);
        sprintf(htitle," Radial Particle Distribution R=%.1f;p_{T,Jet} (GeV/c);r",jetR);
        hRhoNTrackLeadJet[iR] = new TH2F(hname,htitle,150,0,150,25,0,max_r_rho);
        hRhoNTrackLeadJet[iR]->Sumw2();

        sprintf(hname,"hRadialPtLeadJet_R%02d",iR+2);
        sprintf(htitle,"Radial pt_Distribution R=%.1f;1/p_{T,lead jet}d#Sigma p_{T}/dr;r",jetR);

        hRhoPtLeadJet[iR] = new TH2F(hname,htitle,150,0,150,25,0,max_r_rho);
        hRhoPtLeadJet[iR]->Sumw2();

        sprintf(hname,"hRhoNTrackLeadJetBG_R%02d",iR+2);
        sprintf(htitle," Radial Particle Distribution BG R=%.1f;p_{T,Jet} (GeV/c);r",jetR);
        hRhoNTrackLeadJetBG[iR] = new TH2F(hname,htitle,150,0,150,25,0,max_r_rho);
        hRhoNTrackLeadJetBG[iR]->Sumw2();

        sprintf(hname,"hRadialPtLeadJetBG_R%02d",iR+2);
        sprintf(htitle,"Radial pt_Distribution BG R=%.1f;1/p_{T,lead jet}d#Sigma p_{T}/dr;r",jetR);

        hRhoPtLeadJetBG[iR] = new TH2F(hname,htitle,150,0,150,25,0,max_r_rho);
        hRhoPtLeadJetBG[iR]->Sumw2();

        sprintf(hname,"Pt_distribution_R%02d",iR+2);
        sprintf(htitle,"Pt_distribution R=%.1f;p_{T,Lead_Jet} (GeV/c); #eta_{Lead_jet}; #Delta#phi; #Delta#eta",jetR);
        int nbins1[]={30,12,64,30};
        double xmin1[]={0,-3,-TMath::Pi(),-1.5};
        double xmax1[]={150,3,TMath::Pi(),1.5};
        hPtEtaLeadJetDetaDphiP[iR]= new THnF(hname,htitle,4,nbins1,xmin1,xmax1);
        hPtEtaLeadJetDetaDphiP[iR]->GetAxis(3)->SetTitle("#Delta#eta"); // parsing of title string does not seem to work
        hPtEtaLeadJetDetaDphiP[iR]->Sumw2();

        sprintf(hname,"Particle_distribution_R%02d",iR+2);
        sprintf(htitle,"Particle Distribution R=%.1f;p_{T,Lead_Jet} (GeV/c); #eta_{Lead_jet}; #Delta#phi; #Delta#eta",jetR);
        int nbins2[]={30,12,64,30};
        double xmin2[]={0,-3,-TMath::Pi(),-1.5};
        double xmax2[]={150,3,TMath::Pi(),1.5};
        hPtEtaLeadJetDetaDphiN[iR]= new THnF(hname,htitle,4,nbins2,xmin2,xmax2);
        hPtEtaLeadJetDetaDphiN[iR]->GetAxis(3)->SetTitle("#Delta#eta"); // parsing of title string does not seem to work
        hPtEtaLeadJetDetaDphiN[iR]->Sumw2();

        hPtLeadJetRho[iR] = new TH2F(Form("hPtLeadJetRho_R%02d",iR+2),Form("#rho_{BG} vs ptlead jet R=%.1f;p_{T,lead jet}^{unsub} (GeV/c);#rho_{BG} (GeV)"),150,0,150,50,0,50); 
    }

    // get the first event
    HepMC::GenEvent* evt = ascii_in.read_next_event();

    // loop until we run out of events
    while ( evt )
    {
        // analyze the event
        if (debug)
            cout << "Event " << endl;
        // also get the heavy-ion info
        HepMC::HeavyIon* ion = evt->heavy_ion();
        double ep_angle = ion->event_plane_angle();
        if (!evt) 
            cerr << "Input file not found " << inname << endl;


        hNEvent->Fill(0.5,evt->weights()[0]); // count events

        // from example_UsingIterators.cc

        float pt_lead = -1;
        float phi_lead = -100;
        float max_eta_jet = .7;  // ALICE track acceptance cut
        float max_eta_track = .9;
        float min_pt =10;

        int index = 0;
        std::vector <fastjet::PseudoJet> fjInputs;
        for ( HepMC::GenEvent::particle_iterator pit = evt->particles_begin();
                pit != evt->particles_end(); ++pit )
        {
            const HepMC::GenParticle *p = *pit;
            // REDMER: probably change line belwo to only use final state charged particles ? 
            // i.e. if is_charged(p) 
            if (!p->end_vertex() && p->status()==1 && (!charged_jets || is_charged(p)) && fabs(p->momentum().eta()) < max_eta_track);
            { // final state charged particle
                //if (p->status()==1 && is_stable_charged(p) && fabs(p->momentum().eta()) < max_eta) { // charged primary (NB heavy flavour decay products probably not included)
                if (p->momentum().perp() < CO_CUT) continue;  // constituent cut
                hPtPartEta->Fill(p->momentum().perp(), p->momentum().eta(),evt->weights()[0]);
                if (is_charged(p))
                    // so here now just charged tracks
                    hPtChTrackEta->Fill(p->momentum().perp(), p->momentum().eta(),evt->weights()[0]);

                if ( fabs(p->momentum().eta()) < max_eta_track )
                    // also falling in the eta range
                {
                    // bookkeep dphi prof to see what the ep orientatin is
                    // x-axis is dphi, y-axis is pt
                    dPhiTracks->Fill(dphi(p->momentum().phi()-TMath::Pi(), ep_angle), p->momentum().perp(), evt->weights()[0]);


                    if (p->momentum().perp() > pt_lead)
                        // and exceeding leading track pt
                    {
                        pt_lead = p->momentum().perp();
                        phi_lead = p->momentum().phi();
                    }
                    double mom = sqrt(p->momentum().x()*p->momentum().x() + 
                            p->momentum().y()*p->momentum().y() +
                            p->momentum().z()*p->momentum().z());

                    // so here all tracks which fall into the correct eta / charge range 
                    // are pushed to fjInputs as a pseudojet
                    //
                    // if pt-lead is never exceeded then it is still -1
                    fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),mom);
                    jInp.set_user_index(index);
                    fjInputs.push_back(jInp);
                    index++;
                }
            } 

            }   
            if(debug) cout << " Accepted " << index << " track for this event " << endl;

            // Do jet finding
            // Need R =0.2 and R=0.4 later on...
            //
            //
            // i don't understand this syntax
            fastjet::GhostedAreaSpec ghostSpec(max_eta_track,1,0.01);
            fastjet::Strategy               strategy = fastjet::Best;
            fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
            fastjet::AreaType areaType =   fastjet::active_area;
            fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);

            // but this is now theleading pt or the event, not the jet ????
            hPtLead->Fill(pt_lead,evt->weights()[0]);

            float pt_all_track = 0;
            // keepin mind what to do with this variable
            float pt_lead_jet = -1;
            float pt_lead_jet_raw = -1;
            float rho_lead_jet = -1;
            float eta_lead_jet=0;


            for (int iR =0; iR < nR; iR++)
            {
                float jetR = 0.2+0.1*iR;
                // modification 12/05/2014

                pt_lead_jet = -1;
                float phi_lead_jet= 0;
                pt_all_track=0;
                // end modification 12/05/2014 ... 
                fastjet::RangeDefinition range(-max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);

                fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, jetR,recombScheme, strategy);
                fastjet::ClusterSequenceArea clustSeqCh(fjInputs, jetDefCh,areaDef);

                // here we get the jets
                vector <fastjet::PseudoJet> inclusiveJetsCh = clustSeqCh.inclusive_jets();

                fastjet::JetMedianBackgroundEstimator bge;  //.......... Background Sutraction event by event   
                fastjet::ClusterSequenceArea *clustSeqBG = 0;

                // REDMER shouldn't be necessary in this case .. i guess ? 
                /*
                   if (do_bkg) {
                   fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, jetR, recombScheme, strategy);
                   fastjet::Selector  BGSelector = fastjet::SelectorStrip(2*jetR);  //.......... Background Sutraction event by event
                //	BGSelector.set_reference(inclusiveJetsCh[iJet]);
                //	fastjet::ClusterSequenceArea clustSeqBG(fjInputs, jetDefBG,areaDef);

                // BGSelector.set_jets(inclusiveJetsCh);//...............
                clustSeqBG = new fastjet::ClusterSequenceArea(fjInputs, jetDefBG,areaDef);//............
                vector <fastjet::PseudoJet> BGJets = clustSeqBG->inclusive_jets();

                bge.set_selector(BGSelector);
                bge.set_jets(BGJets);
                }*/
                // modification 12/05/2014
                eta_lead_jet=0; 
                float jet_eta= 0;
                float jet_phi= 0;
                float jet_pt= 0;
                float jet_pt_bgsub= 0;

                // end modification 12/05/2014 ...    

                for (unsigned int iJet = 0; iJet < inclusiveJetsCh.size(); iJet++)
                {
                    if (!range.is_in_range(inclusiveJetsCh[iJet]))   
                        continue;

                    // get the jet contituents
                    std::vector<fastjet::PseudoJet> constituents = clustSeqCh.constituents(inclusiveJetsCh[iJet]);
                    double maxpt(0);
                    for(int i(0); i < constituents.size(); i++) {
                        if(constituents[i].perp() > maxpt) maxpt = constituents[i].perp();
                    }
                    if(maxpt < LH_CUT) continue;


                    jet_pt = inclusiveJetsCh[iJet].perp();
                    jet_eta = inclusiveJetsCh[iJet].eta();
                    float dphi_jh = dphi(inclusiveJetsCh[iJet].phi(),phi_lead);

                    float rho = 0;

                    if (do_bkg) {
                        rho= bge.rho(inclusiveJetsCh[iJet]);
                        float jet_area = clustSeqCh.area(inclusiveJetsCh[iJet]);
                        float bgsub = rho*jet_area;
                        jet_pt_bgsub = jet_pt - bgsub;     //.......... Background Sutraction event by event
                    }
                    else 
                        jet_pt_bgsub = jet_pt;

                    // REDMER fill dphi, pt histo
                    dPhiJets->Fill(dphi(inclusiveJetsCh[iJet].phi(),ep_angle), jet_pt,evt->weights()[0]);



                    // modification 12/05/2014
                    jet_phi = inclusiveJetsCh[iJet].phi();
                    // end modification 12/05/2014 ...
                    // modification 16/05/2014 : Leadjet>15GeV...
                    if (jet_pt_bgsub> pt_lead_jet)
                    {
                        pt_lead_jet_raw = jet_pt;
                        pt_lead_jet = jet_pt_bgsub;
                        rho_lead_jet = rho;
                        phi_lead_jet = inclusiveJetsCh[iJet].phi();

                        // modification 12/05/2014
                        eta_lead_jet = inclusiveJetsCh[iJet].eta();
                        // end modification 12/05/2014

                    }

                    // fill histos
                    hPtJetEta[iR]->Fill(jet_pt,jet_eta,evt->weights()[0]);
                    hPtJetEtaBgSub[iR]->Fill(jet_pt_bgsub,jet_eta,evt->weights()[0]);

                    //hPtLeadPtJetDPhi[iR]->Fill(pt_lead,jet_pt,dphi_jh,evt->weights()[0]);

                }  // end of jet loop
                // modification 16/05/2014  // moved outside loop 06/05/2014
                if (pt_lead_jet > 0)
                {
                    hPtLeadJetEta[iR]->Fill(pt_lead_jet,eta_lead_jet,evt->weights()[0]);
                    hPtLeadJetRho[iR]->Fill(pt_lead_jet_raw,rho_lead_jet,evt->weights()[0]);
                    // end modification 16/05/2014

                    // compte the projection of the track Pt onto the leading jet

                    // float pt_all_track = 0;

                    //modification 12/05/2014 
                    // float phi_track = 0;
                    // float eta_track = 0;
                    // end modification 12/05/2014

                    /*
                       for ( HepMC::GenEvent::particle_iterator pit = evt->particles_begin();
                       pit != evt->particles_end(); ++pit )
                       {
                       const HepMC::GenParticle *p = *pit;

                       if (!p->end_vertex() && p->status()==1) {
                       float pt_track = p->momentum().perp();
                       float phi_track = p->momentum().phi();
                       pt_all_track += pt_track * cos( phi_track - phi_lead_jet );
                       if (!charged_constituents || is_charged(p))
                       { 

                    // modification 12/05/2014
                    float eta_track = p->momentum().eta();

                    float dphi_track_leadjet = dphi(phi_track , phi_lead_jet);
                    float deta_track_leadjet = (eta_track - eta_lead_jet);

                    // select bg at 90 degrees
                    if (fabs(fabs(dphi_track_leadjet)- 0.5*TMath::Pi()) < max_r_rho &&
                    fabs(deta_track_leadjet) < max_r_rho && 
                    pt_track > pt_cut_rho)
                    {
                    float R_track_bgjet = TMath::Sqrt((fabs(dphi_track_leadjet)-0.5*TMath::Pi())*(fabs(dphi_track_leadjet)-0.5*TMath::Pi()) + deta_track_leadjet*deta_track_leadjet);
                    if (R_track_bgjet < max_r_rho) {
                    hRhoNTrackLeadJetBG[iR]->Fill(pt_lead_jet,R_track_bgjet,evt->weights()[0]);

                    hRhoPtLeadJetBG[iR]->Fill(pt_lead_jet,R_track_bgjet,(evt->weights()[0] * pt_track/pt_lead_jet));
                    }
                    }
                    // select jet neighbourhood;
                    // Maybe loop over constituents instead?		  
                    if ( fabs(dphi_track_leadjet)<max_r_rho && 
                    fabs(deta_track_leadjet)<max_r_rho ) {
                    float R_track_leadjet = 99;
                    if (pt_track > pt_cut_rho)
                    {
                    R_track_leadjet = TMath::Sqrt(dphi_track_leadjet*dphi_track_leadjet + deta_track_leadjet*deta_track_leadjet);
                    hRhoNTrackLeadJet[iR]->Fill(pt_lead_jet,R_track_leadjet,evt->weights()[0]);

                    hRhoPtLeadJet[iR]->Fill(pt_lead_jet,R_track_leadjet,(evt->weights()[0] * pt_track/pt_lead_jet));
                    }
                    if (R_track_leadjet >90 &&  // calc radius for low pt tracks
                    fabs(dphi_track_leadjet) < jetR && 
                    fabs(deta_track_leadjet) < jetR) {
                    R_track_leadjet = TMath::Sqrt(dphi_track_leadjet*dphi_track_leadjet + deta_track_leadjet*deta_track_leadjet);
                    }

                    if (R_track_leadjet < jetR) {
                    hPtLeadJetPtTrack[iR]->Fill(pt_lead_jet, pt_track, evt->weights()[0]);
                    hPtLeadJetZ[iR]->Fill(pt_lead_jet, pt_track/pt_lead_jet, evt->weights()[0]);
                    hPtLeadJetXi[iR]->Fill(pt_lead_jet, log(pt_lead_jet/pt_track), evt->weights()[0]);
                    }
                    }		

                    double filling1[4] = {pt_lead_jet, eta_lead_jet, dphi_track_leadjet, deta_track_leadjet};
                    hPtEtaLeadJetDetaDphiP[iR]->Fill(filling1,(evt->weights()[0] * pt_track));

                    double filling2[4] = {pt_lead_jet, eta_lead_jet, dphi_track_leadjet, deta_track_leadjet};
                    hPtEtaLeadJetDetaDphiN[iR]->Fill(filling2,evt->weights()[0]);

                    }
                    }
                    }*/
                }
                // charged primary (NB heavy flavour decay products probably not included)
                if (clustSeqBG) {
                    delete clustSeqBG;
                    clustSeqBG = 0;
                }
            }  // R loop

            // Use pt_all_track from largest R (R=0.5)
            hPtAllTrack->Fill(pt_all_track,evt->weights()[0]);

            // modification 12/05/2014 ....modification 16/05/2014 (3D pt-lead_jet,pt-all_track,eta)

            hPtAllTrackLeadJetEta->Fill(pt_lead_jet,pt_all_track,eta_lead_jet,evt->weights()[0]);

            // end modification 12/05/2014


            if (debug)
                cout << "Leading particle pt " << pt_lead << " phi " << phi_lead << endl; 
            // delete the created event from memory
            delete evt;
            // read the next event
            ascii_in >> evt;    
        }

        fout.Write();

        for (int iR =0; iR < nR; iR++) {
            hPtEtaLeadJetDetaDphiP[iR]->Write();
            hPtEtaLeadJetDetaDphiN[iR]->Write();
        }

        fout.Close();

        return 0;
    } 

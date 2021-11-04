//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 14 17:03:33 2020 by ROOT version 6.20/04
// from TTree tau3x1/
// found on file: 1x1_taupair.root
//////////////////////////////////////////////////////////

#ifndef DataMC_h
#define DataMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DataMC {
    public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    
    // Fixed size dimensions of array or collections stored in the TTree if any.
    
    // Declaration of leaf types
    Double_t        nPhotons_negThrust;
    Double_t        nPhotons_posThrust;
    Double_t        nPi0s_negThrust;
    Double_t        nPi0s_posThrust;
    Double_t        thrust;
    Double_t        visibleEnergyOfEventCMS;
    Double_t        tauPlusMCMode;
    Double_t        tauMinusMCMode;
    Double_t        tauPlusMCProng;
    Double_t        tauMinusMCProng;
    Double_t        track_negThrust_E_CMS;
    Double_t        track_negThrust_px_CMS;
    Double_t        track_negThrust_py_CMS;
    Double_t        track_negThrust_pz_CMS;
    Double_t        track_negThrust_EoverP;
    Double_t        track_posThrust_E_CMS;
    Double_t        track_posThrust_px_CMS;
    Double_t        track_posThrust_py_CMS;
    Double_t        track_posThrust_pz_CMS;
    Double_t        track_posThrust_EoverP;
    Double_t        Ymin2;
    Double_t        Ymax2;
    
    // List of branches
    TBranch        *b_nPhotons_negThrust;   //!
    TBranch        *b_nPhotons_posThrust;   //!
    TBranch        *b_nPi0s_negThrust;   //!
    TBranch        *b_nPi0s_posThrust;   //!
    TBranch        *b_thrust;   //!
    TBranch        *b_visibleEnergyOfEventCMS;   //!
    TBranch        *b_tauPlusMCMode;   //!
    TBranch        *b_tauMinusMCMode;   //!
    TBranch        *b_tauPlusMCProng;   //!
    TBranch        *b_tauMinusMCProng;   //!
    TBranch        *b_track_negThrust_E_CMS;   //!
    TBranch        *b_track_negThrust_px_CMS;   //!
    TBranch        *b_track_negThrust_py_CMS;   //!
    TBranch        *b_track_negThrust_pz_CMS;   //!
    TBranch        *b_track_negThrust_EoverP;   //!
    TBranch        *b_track_posThrust_E_CMS;   //!
    TBranch        *b_track_posThrust_px_CMS;   //!
    TBranch        *b_track_posThrust_py_CMS;   //!
    TBranch        *b_track_posThrust_pz_CMS;   //!
    TBranch        *b_track_posThrust_EoverP;   //!
    TBranch        *b_Ymin2;   //!
    TBranch        *b_Ymax2;   //!
    
    DataMC(TTree *tree=0);
    virtual ~DataMC();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    //virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef DataMC_cxx
DataMC::DataMC(TTree *tree) : fChain(0)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("1x1_taupair.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("1x1_taupair.root");
        }
        f->GetObject("tau3x1",tree);
        
    }
    Init(tree);
}

DataMC::~DataMC()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t DataMC::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t DataMC::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void DataMC::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("nPhotons_negThrust", &nPhotons_negThrust, &b_nPhotons_negThrust);
    fChain->SetBranchAddress("nPhotons_posThrust", &nPhotons_posThrust, &b_nPhotons_posThrust);
    fChain->SetBranchAddress("nPi0s_negThrust", &nPi0s_negThrust, &b_nPi0s_negThrust);
    fChain->SetBranchAddress("nPi0s_posThrust", &nPi0s_posThrust, &b_nPi0s_posThrust);
    fChain->SetBranchAddress("thrust", &thrust, &b_thrust);
    fChain->SetBranchAddress("visibleEnergyOfEventCMS", &visibleEnergyOfEventCMS, &b_visibleEnergyOfEventCMS);
    fChain->SetBranchAddress("tauPlusMCMode", &tauPlusMCMode, &b_tauPlusMCMode);
    fChain->SetBranchAddress("tauMinusMCMode", &tauMinusMCMode, &b_tauMinusMCMode);
    fChain->SetBranchAddress("tauPlusMCProng", &tauPlusMCProng, &b_tauPlusMCProng);
    fChain->SetBranchAddress("tauMinusMCProng", &tauMinusMCProng, &b_tauMinusMCProng);
    fChain->SetBranchAddress("track_negThrust_E_CMS", &track_negThrust_E_CMS, &b_track_negThrust_E_CMS);
    fChain->SetBranchAddress("track_negThrust_px_CMS", &track_negThrust_px_CMS, &b_track_negThrust_px_CMS);
    fChain->SetBranchAddress("track_negThrust_py_CMS", &track_negThrust_py_CMS, &b_track_negThrust_py_CMS);
    fChain->SetBranchAddress("track_negThrust_pz_CMS", &track_negThrust_pz_CMS, &b_track_negThrust_pz_CMS);
    fChain->SetBranchAddress("track_negThrust_EoverP", &track_negThrust_EoverP, &b_track_negThrust_EoverP);
    fChain->SetBranchAddress("track_posThrust_E_CMS", &track_posThrust_E_CMS, &b_track_posThrust_E_CMS);
    fChain->SetBranchAddress("track_posThrust_px_CMS", &track_posThrust_px_CMS, &b_track_posThrust_px_CMS);
    fChain->SetBranchAddress("track_posThrust_py_CMS", &track_posThrust_py_CMS, &b_track_posThrust_py_CMS);
    fChain->SetBranchAddress("track_posThrust_pz_CMS", &track_posThrust_pz_CMS, &b_track_posThrust_pz_CMS);
    fChain->SetBranchAddress("track_posThrust_EoverP", &track_posThrust_EoverP, &b_track_posThrust_EoverP);
    fChain->SetBranchAddress("Ymin2", &Ymin2, &b_Ymin2);
    fChain->SetBranchAddress("Ymax2", &Ymax2, &b_Ymax2);
    Notify();
}

Bool_t DataMC::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}

void DataMC::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t DataMC::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
//#endif // #ifdef DataMC_cxx


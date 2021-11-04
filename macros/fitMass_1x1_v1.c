////////////////////////////////////////////////////
//
// E.C.B 2020
////////////////////////////////////////////////////
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

#include "DataMC.h"

using namespace RooFit;


void Data_filter(DataMC *t, RooDataSet *data, RooRealVar x, Int_t mode, Int_t type);
void GetSandB(DataMC *t, Double_t &S, Double_t &B, Double_t bdt);
void fit();
void PlotSignificance();
void SetReader();

TGraph *gr;

TMVA::Reader *reader; //= new TMVA::Reader( "!Color:!Silent" );
Float_t v[12]; //for reader



// Zoom
Int_t bins = 50;


Double_t Ymin_min = 1.755;
//Double_t Ymin_max = 1.81;

//Double_t Ymin_min = 1.73;
Double_t Ymin_max = 1.85;


Double_t Ymax_min = 1.7;
Double_t Ymax_max = 1.85;

//Double_t Ymax_min = 0.0;
//Double_t Ymax_max = 5.0;


/*
 
 Double_t Ymin_min = 0.5;
 Double_t Ymin_max = 2.5;
 
 Double_t Ymax_min = 0;
 Double_t Ymax_max = 7.;
 */

void fitMass_1x1_v1()
{
    SetReader();
    PlotSignificance();
    fit();
    
    return;
}

void SetReader()
{
    
    reader = new TMVA::Reader( "!Color:Silent" );
    reader->AddVariable( "v1", &v[0] );
    reader->AddVariable( "v2", &v[1] );
    reader->AddVariable( "v3", &v[2] );
    reader->AddVariable( "v4", &v[3] );
    reader->AddVariable( "v5", &v[4] );
    reader->AddVariable( "v6", &v[5] );
    reader->AddVariable( "v7", &v[6] );
    reader->AddVariable( "v8", &v[7] );
    reader->AddVariable( "v9", &v[8] );
    reader->AddVariable( "v10", &v[9] );
    //reader->AddVariable( "v11", &v[10] );
    //reader->AddVariable( "v12", &v[11] );
    
    reader->BookMVA("BDT", "dataset/weights/_MyRBDT.weights.xml");
    return;
}



void PlotSignificance()
{
    
    //Adding data in a TChain
    TChain *chData = new TChain("tau3x1");
    chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_taupair.root");
    chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_qqbar_bkg.root");
    chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_mumu1_bkg.root");
    chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_mumu2_bkg.root");
    chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_ee_bkg.root");
    //chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_eemumu_bkg.root");
    chData->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_eeee_bkg.root");
    
    
    TTree *data = (TTree*) chData;
    DataMC *t = new DataMC(data);
    
    Double_t bdt_min = -0.;
    Double_t bdt_max = 0.4;
    Double_t step = 0.05;
    Double_t bdtval = bdt_min;
    
    Double_t Ns = 0;
    Double_t Nb = 0;
    Double_t x[30];
    Double_t y[30];
    Int_t i = 0;
    
    while(bdtval<bdt_max)
    {
        GetSandB(t,Ns,Nb,bdtval);
        //cout<<Ns<<"  "<<Nb<<endl;
        if((Ns+Nb)!=0)
        {
            x[i] = bdtval;
            //y[i] = Ns/TMath::Sqrt(Ns + Nb);
            y[i] = TMath::Sqrt(Ns + Nb) - TMath::Sqrt(Nb);
            i++;
            bdtval += step;
        }
        else bdtval=2;
    }
    Int_t n = i;
    
    
    gr = new TGraph(n,x,y);
    gr->Draw("AC*");
    
    //delete chData;
    //delete data;
    
}



void fit()
{
    
    
    //TChain *chData = new TChain("tree");
    //chData->Add("../data/*.root");
    
    //Adding data in a TChain
    TChain *chData_taupair = new TChain("tau3x1");
    chData_taupair->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_taupair.root");
    TTree *taudata = (TTree*) chData_taupair;
    
    
    TChain *chData_bkg = new TChain("tau3x1");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_qqbar_bkg.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_mumu1_bkg.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_mumu2_bkg.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_ee_bkg.root");
    //chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_eemumu_bkg.root");
    chData_bkg->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_eeee_bkg.root");
    
    
    TTree *treeS = (TTree*) chData_taupair;
    DataMC *tS = new DataMC(treeS);
    
    TTree *treeB = (TTree*) chData_bkg;
    DataMC *tB = new DataMC(treeB);
    
    RooRealVar Ymin("Ymin", "Y_{min}", Ymin_min, Ymin_max, "GeV/c^{2}");
    RooDataSet *dataYmin = new RooDataSet("dataYmin","dataYmin",RooArgSet(Ymin));
    RooDataSet *dataYminP = new RooDataSet("dataYminP","dataYminP",RooArgSet(Ymin));
    RooDataSet *dataYminB = new RooDataSet("dataYminB","dataYminB",RooArgSet(Ymin));
    
    RooRealVar Ymax("Ymax", "Y_{max}", Ymax_min, Ymax_max, "GeV/c^{2}");
    RooDataSet *dataYmax = new RooDataSet("dataYmax","dataYmax",RooArgSet(Ymax));
    RooDataSet *dataYmaxP = new RooDataSet("dataYmaxP","dataYmaxP",RooArgSet(Ymax));
    RooDataSet *dataYmaxB = new RooDataSet("dataYmaxB","dataYmaxB",RooArgSet(Ymax));
    
    RooCategory tagCat("tagCat","tagging Category");
    tagCat.defineType("signal");
    tagCat.defineType("taupair");
    tagCat.defineType("bkg");
    
    
    tagCat.setLabel("signal");
    
    Data_filter(tS,dataYmin,Ymin,1,1);
    dataYmin->addColumn(tagCat);
    
    //Data_filter(tS,dataYmax,Ymax,2,1);
    //dataYmax->addColumn(tagCat);
    
    
    tagCat.setLabel("taupair");
    
    Data_filter(tS,dataYminP,Ymin,1,0);
    dataYminP->addColumn(tagCat);
    
    //Data_filter(tS,dataYmaxP,Ymax,2,0);
    //dataYmaxP->addColumn(tagCat);
    
    
    //dataMmin->Print("v");
    //dataYmin->Print("v");
    //dataYmax->Print("v");
    
    
    tagCat.setLabel("bkg");
    
    Data_filter(tB,dataYminB,Ymin,1,0);
    dataYminB->addColumn(tagCat);
    
    //Data_filter(tB,dataYmaxB,Ymax,2,0);
    //dataYmaxB->addColumn(tagCat);
    
    //dataMminB->Print("v");
    //dataYminB->Print("v");
    //dataYmaxB->Print("v");
    
    
    dataYmin->append(*dataYminP);
    //dataYmax->append(*dataYmaxP);
    
    dataYmin->append(*dataYminB);
    //dataYmax->append(*dataYmaxB);
    
    
    //RooRealVar P11("P11","P11", 1.778, 0, 1000.0);
    //RooRealVar P21("P21","P21", 0.01, -1000, 1000);
    RooRealVar P11("P11","P11", 1.778, 1.76, 1.8);
    RooRealVar P21("P21","P21", 0.01, -1, 1);
    RooRealVar P31("P31","P31", -0.42/3.15, -100, 100);
    RooRealVar P41("P41","P41", -0.10/3.15, -50.0, 50.0);
    RooRealVar P51("P51","P51", -1.0/3.15, -50, 50);
    RooRealVar P61("P61","P61", 0, -100, 100);
    RooRealVar P71("P71","P71", 0, -100, 100);
    RooRealVar P81("P81","P81", 0, -100, 100);
    RooRealVar P91("P91","P91", 0, -100, 100);
    
    RooRealVar P12("P12","P12", 1.778, 0, 1000.0);
    RooRealVar P22("P22","P22", 0.01, -1000, 1000);
    RooRealVar P32("P32","P32", -0.42/3.15, -100, 100);
    RooRealVar P42("P42","P42", -0.10/3.15, -50.0, 50.0);
    RooRealVar P52("P52","P52", -1.0/3.15, -50, 50);
    RooRealVar P62("P62","P62", 0, -100, 100);
    RooRealVar P72("P72","P72", 0, -100, 100);
    RooRealVar P82("P82","P82", 0, -100, 100);
    
    //P11.setVal(1.778);
    //P21.setVal(0.01);
    //P31.setVal(-0.01);
    
    
    
    P12.setVal(1.775);
    P22.setVal(-0.006);
    P32.setVal(-0.175);
    P42.setVal(-0.04);
    P52.setVal(-0.5);
    P62.setVal(0.008);
    P72.setVal(0.01);
    P82.setVal(0.01);
    
    //P81.setVal(0);
    //P81.setConstant(1);
    
    //TF1 flow("flow","[0]*TMath::Erfc((x-[2])/(TMath::Sqrt(2)*[1]))+[3]*x+[4]",xvalmin,xvalmax);
    //TF1 fup("fup","[0]*TMath::Erf((x-[2])/(TMath::Sqrt(2)*[1]))+[3]*x+[4]",yvalmin,yvalmax);
    
    //RooGenericPdf MassModelY("PsMassPDFY", "(@3 + @4*@0 + @6*@0*@0 + @7*@0*@0*@0 + @8*@0*@0*@0*@0) * atan((@0 - @1)/@2) + @5 * @0 + 1",
    //               RooArgSet(Ymin,P11,P21,P31,P41,P51,P61,P71,P81));
    
    //RooGenericPdf MassModelY("PsMassPDFY", "(@3 + @4*@0 + @6*@0*@0) * TMath::Erfc((@0 - @1)/@2) + @7*@0*@0* + @5 * @0 + 1",
    RooGenericPdf MassModelY("PsMassPDFY", "(@3 + @4*@0 + @6*@0*@0 + @7*@0*@0*@0) * TMath::Erfc((@0 - @1)/@2) + @8*@0*@0 + @5 * @0 + 1",
                             RooArgSet(Ymin,P11,P21,P31,P41,P51,P61,P71,P81));
    
    RooGenericPdf MassModelY2("PsMassPDFY2", "(@3 + @4*@0 + @6*@0*@0 + @7*@0*@0*@0 + @8*@0*@0*@0*@0) * atan((@0 - @1)/@2) + @5 * @0 + 1",
                              RooArgSet(Ymax,P12,P22,P32,P42,P52,P62,P72,P82));
    
    
    //RooFitResult *fitmas = MassModel.fitTo(*dataMmin, Save(kTRUE), Strategy(1), NumCPU(4));
    
    RooFitResult *fitmasY = MassModelY.fitTo(*dataYmin, Save(kTRUE), Strategy(1), NumCPU(2));
    
    //RooFitResult *fitmasY2 = MassModelY2.fitTo(*dataYmax, Save(kTRUE), Strategy(1), NumCPU(4));
    
    //fitmas->Print("v");
    fitmasY->Print("v");
    //fitmasY2->Print("v");
    

    
    
    //TCanvas *cMmin = new TCanvas("cMmin","cMmin",800,800);
    //RooPlot* thrFrame = Mmin.frame(Mmin_min,Mmin_max,bins);
    //dataMmin->plotOn(thrFrame, Name("PseudoMass"), DataError(RooAbsData::SumW2));
    //dataMmin->plotOn(thrFrame, Name("Mmin"), DataError(RooAbsData::SumW2),Cut("tagCat==tagCat::signal"));
    //dataMmin->plotOn(thrFrame, Name("Mmin"), DataError(RooAbsData::SumW2),Cut("tagCat==tagCat::taupair"));
    //dataMmin->plotOn(thrFrame, Name("Mmin"), DataError(RooAbsData::SumW2),Cut("tagCat==tagCat::bkg"));
    //MassModel.plotOn(thrFrame, Name("model"));
    //thrFrame->GetXaxis()->SetTitleSize(0);
    //thrFrame->Draw();
    
    TCanvas *cYmin = new TCanvas("cYmin","cYmin",800,800);
    RooPlot* thrFrameY = Ymin.frame(Ymin_min,Ymin_max,bins);
    dataYmin->plotOn(thrFrameY, Name("Ymin"), DataError(RooAbsData::SumW2));
    dataYmin->plotOn(thrFrameY, Name("Ymin"), DataError(RooAbsData::SumW2),MarkerColor(kBlue),Cut("tagCat==tagCat::signal"));
    dataYmin->plotOn(thrFrameY, Name("Ymin"), DataError(RooAbsData::SumW2), MarkerColor(kRed),Cut("tagCat==tagCat::taupair"));
    dataYmin->plotOn(thrFrameY, Name("Ymin"), DataError(RooAbsData::SumW2),MarkerColor(24),Cut("tagCat==tagCat::bkg"));
    MassModelY.plotOn(thrFrameY, Name("model"));
    MassModelY.paramOn(thrFrameY,Layout(0.6,0.9,0.9));
    thrFrameY->GetXaxis()->SetTitleSize(0);
    thrFrameY->Draw();
    /*
     TCanvas *cYmax = new TCanvas("cYmax","cYmax",800,800);
     RooPlot* thrFrameYmax = Ymax.frame(Ymax_min,Ymax_max,bins);
     dataYmax->plotOn(thrFrameYmax, Name("Ymax"), DataError(RooAbsData::SumW2));
     dataYmax->plotOn(thrFrameYmax, Name("Ymax"), DataError(RooAbsData::SumW2),MarkerColor(kBlue),Cut("tagCat==tagCat::signal"));
     dataYmax->plotOn(thrFrameYmax, Name("Ymax"), DataError(RooAbsData::SumW2),MarkerColor(kRed),Cut("tagCat==tagCat::taupair"));
     dataYmax->plotOn(thrFrameYmax, Name("Ymax"), DataError(RooAbsData::SumW2),MarkerColor(24),Cut("tagCat==tagCat::bkg"));
     //MassModelY2.plotOn(thrFrameYmax, Name("model"));
     //MassModelY2.paramOn(thrFrameYmax,Layout(0.65,0.95,0.95));
     thrFrameYmax->GetXaxis()->SetTitleSize(0);
     thrFrameYmax->Draw();
     */
    
    
    
}

void Data_filter(DataMC *t, RooDataSet *data, RooRealVar x, Int_t mode, Int_t type)
{
    Double_t xmin,xmax;
    if(mode == 1)
    {
        xmin = Ymin_min;
        xmax = Ymin_max;
    }
    if(mode == 2)
    {
        xmin = Ymax_min;
        xmax = Ymax_max;
    }
    /*
     TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
     Float_t v[12];
     
     reader->AddVariable( "v1", &v[0] );
     reader->AddVariable( "v2", &v[1] );
     reader->AddVariable( "v3", &v[2] );
     reader->AddVariable( "v4", &v[3] );
     reader->AddVariable( "v5", &v[4] );
     reader->AddVariable( "v6", &v[5] );
     reader->AddVariable( "v7", &v[6] );
     reader->AddVariable( "v8", &v[7] );
     reader->AddVariable( "v9", &v[8] );
     reader->AddVariable( "v10", &v[9] );
     //reader->AddVariable( "v11", &v[10] );
     //reader->AddVariable( "v12", &v[11] );
     
     reader->BookMVA("BDT", "dataset/weights/_BDT.weights.xml");
     */
    Long64_t nentries = t->fChain->GetEntries();
    cout<<"No. events in data: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    
    Bool_t IsSignal = false;
    
    for (Long64_t jentry=0; jentry!=nentries;jentry++)
    {
        Long64_t ientry = t->LoadTree(jentry);
        if (ientry < 0) break;
        nb = t->fChain->GetEntry(jentry);
        nbytes += nb;
        if(jentry%(nentries/10)==0) cout<<(100.0*jentry)/nentries<<" %"<<endl;
        
        if(t->visibleEnergyOfEventCMS<2.0) continue;
        if(t->visibleEnergyOfEventCMS>10.0) continue;
        if(t->thrust>0.99) continue;
        if(t->nPi0s_negThrust>0) continue;
        if(t->nPi0s_posThrust>0) continue;
        //if(t->nPhotons_negThrust>0) continue;
        //if(t->nPhotons_posThrust>0) continue;
        
        //if(t->track_negThrust_EoverP > 0.8 ||
        // t->track_posThrust_EoverP > 0.8 ) continue;
        
        
        
        if(t->tauMinusMCMode==3 && t->tauPlusMCMode==3) IsSignal=true;
        else IsSignal=false;
        
        v[0] = t->thrust;
        v[1] = t->visibleEnergyOfEventCMS;
        v[2] = t->track_negThrust_EoverP;
        v[3] = t->track_posThrust_EoverP;
        v[4] = t->nPhotons_negThrust;
        v[5] = t->nPhotons_posThrust;
        
        /*
         v[0] = t->thrust;
         v[1] = t->visibleEnergyOfEventCMS;
         v[2] = t->track_negThrust_EoverP;
         v[3] = t->track_posThrust_EoverP;
         v[4] = t->nPi0s_negThrust;
         v[5] = t->nPi0s_posThrust;
         v[6] = t->nPhotons_negThrust;
         v[7] = t->nPhotons_posThrust;
         v[8] = TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS);
         v[9] = TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS);
         v[10] = TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS + t->track_negThrust_pz_CMS*t->track_negThrust_pz_CMS);
         v[11] = TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS + t->track_posThrust_pz_CMS*t->track_posThrust_pz_CMS);
         */
        
        Double_t pt1 = TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS);
        Double_t pt2 =  TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS);
        
        Double_t p1 =  TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS + t->track_negThrust_pz_CMS*t->track_negThrust_pz_CMS);
        Double_t p2 = TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS + t->track_posThrust_pz_CMS*t->track_posThrust_pz_CMS);
        
        if(pt1>pt2)
        {
            v[6] = pt1;
            v[7] = pt2;
        }
        else
        {
            v[6] = pt2;
            v[7] = pt1;
        }
        
        if(p1>p2)
        {
            v[8] = pt1;
            v[9] = pt2;
        }
        else
        {
            v[8] = pt2;
            v[9] = pt1;
        }
        
        
        
        Double_t val = reader->EvaluateMVA("BDT" );
        //cout<< " BDT : "<<val<<endl;
        
        if(val<0.2) continue;
        if(type==1 && IsSignal==false) continue;
        
        Double_t xval = 0;
        if(mode==1) xval = TMath::Sqrt(t->Ymin2);
        if(mode==2) xval = TMath::Sqrt(t->Ymax2);
        if(mode!=1 && mode!=2)
        {
            cout<<" Mode no valid "<<endl;
            continue;
        }
        
        if(xval > xmin && xval < xmax)
        {
            x = xval;
            data->add(x);
        }
        
        //if(!IsSignal) cout<<t->tauMinusMCMode<<"  "<<t->tauPlusMCMode<<endl;
    }
    
    
    //delete reader;
    return;
}

void GetSandB(DataMC *t, Double_t &S, Double_t &B, Double_t bdt)
{
    
    
    Long64_t nentries = t->fChain->GetEntries();
    cout<<"No. events in data: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    
    Bool_t IsSignal = false;
    Double_t nSig = 0;
    Double_t nBkg = 0;
    
    for (Long64_t jentry=0; jentry!=nentries;jentry++)
    {
        Long64_t ientry = t->LoadTree(jentry);
        if (ientry < 0) break;
        nb = t->fChain->GetEntry(jentry);
        nbytes += nb;
        //if(jentry%(nentries/10)==0) cout<<(100.0*jentry)/nentries<<" %"<<endl;
        
        if(t->visibleEnergyOfEventCMS<2.0) continue;
        if(t->visibleEnergyOfEventCMS>10.0) continue;
        if(t->thrust>0.99) continue;
        if(t->nPi0s_negThrust>0) continue;
        if(t->nPi0s_posThrust>0) continue;
        //if(t->nPhotons_negThrust>0) continue;
        //if(t->nPhotons_posThrust>0) continue;
        
        if(t->tauMinusMCMode==3 && t->tauPlusMCMode==3) IsSignal=true;
        else IsSignal=false;
        
        v[0] = t->thrust;
        v[1] = t->visibleEnergyOfEventCMS;
        v[2] = t->track_negThrust_EoverP;
        v[3] = t->track_posThrust_EoverP;
        v[4] = t->nPhotons_negThrust;
        v[5] = t->nPhotons_posThrust;
        
        /*
         v[0] = t->thrust;
         v[1] = t->visibleEnergyOfEventCMS;
         v[2] = t->track_negThrust_EoverP;
         v[3] = t->track_posThrust_EoverP;
         v[4] = t->nPi0s_negThrust;
         v[5] = t->nPi0s_posThrust;
         v[6] = t->nPhotons_negThrust;
         v[7] = t->nPhotons_posThrust;
         v[8] = TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS);
         v[9] = TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS);
         v[10] = TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS + t->track_negThrust_pz_CMS*t->track_negThrust_pz_CMS);
         v[11] = TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS + t->track_posThrust_pz_CMS*t->track_posThrust_pz_CMS);
         */
        
        Double_t pt1 = TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS);
        Double_t pt2 =  TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS);
        
        Double_t p1 =  TMath::Sqrt(t->track_negThrust_px_CMS*t->track_negThrust_px_CMS + t->track_negThrust_py_CMS*t->track_negThrust_py_CMS + t->track_negThrust_pz_CMS*t->track_negThrust_pz_CMS);
        Double_t p2 = TMath::Sqrt(t->track_posThrust_px_CMS*t->track_posThrust_px_CMS + t->track_posThrust_py_CMS*t->track_posThrust_py_CMS + t->track_posThrust_pz_CMS*t->track_posThrust_pz_CMS);
        
        if(pt1>pt2)
        {
            v[6] = pt1;
            v[7] = pt2;
        }
        else
        {
            v[6] = pt2;
            v[7] = pt1;
        }
        
        if(p1>p2)
        {
            v[8] = pt1;
            v[9] = pt2;
        }
        else
        {
            v[8] = pt2;
            v[9] = pt1;
        }
        
        
        
        Double_t val = reader->EvaluateMVA("BDT" );
        //cout<< " BDT : "<<val<<endl;
        
        if(val<bdt) continue;
        
        if(IsSignal==true) nSig++;
        else nBkg++;
        
    }
    
    
    cout<<bdt<<":  S = "<<nSig<<"  :  B = "<<nBkg<<endl;
    S = nSig;
    B = nBkg;
    
    //delete reader;
    return;
    
    
    
}


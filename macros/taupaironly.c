#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TH1.h"


#include "DataMC.h"
using namespace RooFit;

void Data_filtro(DataMC *t, RooDataSet *data, RooRealVar x, Int_t mode, Int_t type);
void plot();

Int_t bins = 100;
Double_t Ymin_min = 1.6;
Double_t Ymin_max = 1.85;

Double_t Ymax_min = 1.6;
Double_t Ymax_max = 1.85;


void taupaironly()
{
    plot();
}

void plot()
{
    
    TChain *chData_taupair = new TChain("tau3x1");
    chData_taupair->Add("~mac/Documents/Johan/Belle2/Tau/newmethod/taumass/data/1x1_taupair.root");
    
    TTree *taudata = (TTree*) chData_taupair;
    
    TTree *treeS = (TTree*) chData_taupair;
    DataMC *tS = new DataMC(treeS);
    
    RooRealVar Ymin("Ymin", "Y_{min}", Ymin_min, Ymin_max, "GeV/c^{2}");
    RooDataSet *dataYmin = new RooDataSet("dataYmin","dataYmin",RooArgSet(Ymin));
    
    //RooRealVar Ymax("Ymax", "Y_{max}", Ymax_min, Ymax_max, "GeV/c^{2}");
    //RooDataSet *dataYmax = new RooDataSet("dataYmax","dataYmax",RooArgSet(Ymax));
    
    RooCategory tagCat("tagCat","tagging Category");
    tagCat.defineType("taupair");

    tagCat.setLabel("taupair");
    Data_filtro(tS,dataYmin,Ymin,1,1);
    dataYmin->addColumn(tagCat);
    
    
    
    //Data_filtro(tS,dataYmax,Ymax,2,1);
    //dataYmax->addColumn(tagCat);
    
    
    
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
    
    //RooGenericPdf MassModelY2("PsMassPDFY2", "(@3 + @4*@0 + @6*@0*@0 + @7*@0*@0*@0 + @8*@0*@0*@0*@0) * atan((@0 - @1)/@2) + @5 * @0 + 1",
                              //RooArgSet(Ymax,P12,P22,P32,P42,P52,P62,P72,P82));
    
    
    //RooFitResult *fitmas = MassModel.fitTo(*dataMmin, Save(kTRUE), Strategy(1), NumCPU(4));
    
    RooFitResult *fitmasY = MassModelY.fitTo(*dataYmin, Save(kTRUE), Strategy(1), NumCPU(2));
    
    //RooFitResult *fitmasY2 = MassModelY2.fitTo(*dataYmax, Save(kTRUE), Strategy(1), NumCPU(4));
    
    //fitmas->Print("v");
    fitmasY->Print("v");
    //fitmasY2->Print("v");
    
    TCanvas *cYmin = new TCanvas("cYmin","cYmin",800,800);
    RooPlot* thrFrameY = Ymin.frame(Ymin_min,Ymin_max,bins);
    dataYmin->plotOn(thrFrameY, Name("Ymin"), DataError(RooAbsData::SumW2));
    dataYmin->plotOn(thrFrameY, Name("Ymin"), DataError(RooAbsData::SumW2),MarkerColor(kBlue),Cut("tagCat==tagCat::taupair"));
    MassModelY.plotOn(thrFrameY, Name("model"));
    MassModelY.paramOn(thrFrameY,Layout(0.6,0.9,0.9));
    thrFrameY->GetXaxis()->SetTitleSize(0);
    thrFrameY->Draw();

    
}




void Data_filtro(DataMC *t, RooDataSet *data, RooRealVar x, Int_t mode, Int_t type)
{
    Double_t xmin,xmax;
    if(mode == 1)
    {
        xmin = Ymin_min;
        xmax = Ymin_max;
    }
    
    Long64_t nentries = t->fChain->GetEntries();
    cout<<"No. events in data: "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;
    
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
        //if(t->nPi0s_negThrust==3) continue;
        if(t->nPi0s_negThrust>0) continue;
        if(t->nPi0s_posThrust>0) continue;
        //if(t->nPi0s_posThrust==3) continue;
        
        
        if(t->tauMinusMCMode!=3 || t->tauPlusMCMode!=3) continue;
        
        //if(t->tauMinusMCProng==1 && t->tauPlusMCProng==1) continue;
        
        Double_t px1,py1,pz1,px2,py2,pz2,E1,E2,Ecm;
        Double_t z_pi1, z_pi2, a_pi1x, a_pi1y, a_pi1z, a_pi2x, a_pi2y, a_pi2z;
        Double_t punto_a_b,A_0,B_0,C_0,D_0,M_tau_Edge,M_nu_Edge,M_tau_Max;
        
        Ecm = 10.5794;//Gev
        px1 = t->track_posThrust_px_CMS;
        py1 = t->track_posThrust_py_CMS;
        pz1 = t->track_posThrust_pz_CMS;
        E1 = t->track_posThrust_E_CMS;
        
        px2 = t->track_negThrust_px_CMS;
        py2 = t->track_negThrust_py_CMS;
        pz2 = t->track_negThrust_pz_CMS;
        E2 = t->track_negThrust_E_CMS;
        
        //Normaliced energies (pion neutrino)
        z_pi1 = (E1)/(Ecm);
        z_pi2 = (E2)/(Ecm);
        
        //normaliced moments in components
        a_pi1x = px1/Ecm;
        a_pi1y = py1/Ecm;
        a_pi1z = pz1/Ecm;
        a_pi2x = px2/Ecm;
        a_pi2y = py2/Ecm;//
        a_pi2z = pz2/Ecm;
        
        
        // create the moments vectors
        TVector3 a_pi1, a_pi2, sum_a_b,cruz_a_b,comp1,comp2;
        a_pi1.SetXYZ(a_pi1x,a_pi1y,a_pi1z);
        a_pi2.SetXYZ(a_pi2x,a_pi2y,a_pi2z);
        
        // define some variables previous for add vectors, scalar and vectorial products
        
        sum_a_b = a_pi1+a_pi2;
        punto_a_b = a_pi1.Dot(a_pi2);
        cruz_a_b = a_pi1.Cross(a_pi2);
        comp1 = ((z_pi2*z_pi2) - z_pi2)*a_pi1 + ((z_pi1*z_pi1) - z_pi1)*a_pi2;
        
        // Definiendo los parametroa A0, B0, C0 y D0
        A_0 = sum_a_b.Mag2();//magnitud al cuadrado
        B_0 = (2.0*a_pi1.Mag2())*((z_pi2*z_pi2) - z_pi2) + (2.0*a_pi2.Mag2())*((z_pi1*z_pi1) - z_pi1) + (2.0*punto_a_b)*((z_pi1*z_pi1) + (z_pi2*z_pi2) - z_pi1 - z_pi2 - (sum_a_b.Mag2()));
        C_0 = 4.0*cruz_a_b.Mag2();
        D_0 = a_pi1.Mag2()*a_pi2.Mag2()*sum_a_b.Mag2() - cruz_a_b.Mag2() +(2.0*a_pi1.Mag2()*a_pi2.Mag2())*(z_pi1 + z_pi2 - (z_pi1*z_pi1) - (z_pi2*z_pi2)) + comp1.Mag2() - (2.0*punto_a_b)*(a_pi1.Mag2()*((z_pi2*z_pi2) - z_pi2) + a_pi2.Mag2()*((z_pi1*z_pi1) - z_pi1));
        // masa límite para el tau y el neutrino
        Double_t mu_tau_Edge2 = (4.0*(B_0*B_0) + 3.0*(C_0*C_0) - 16.0*(A_0*D_0) - 8.0*(B_0*C_0))/(16.0*(A_0*C_0));
        Double_t mu_nu_Edge2 = (4.0*(B_0*B_0) - (C_0*C_0) - 16.0*(A_0*D_0))/(16*(A_0*C_0));
        Double_t mu_tau2_Edge2 = ((TMath::Sqrt((B_0*B_0) - 4.0*(A_0*D_0))) - B_0)/(2.0*(A_0));
        
        Double_t mu_tau_Max2 = ((B_0 - C_0)*(B_0 - C_0))/(4.0*(A_0*C_0)) - ((D_0)/(C_0));
        if (mu_nu_Edge2>=0)
        {
            M_tau_Edge = TMath::Sqrt(mu_tau_Edge2)*Ecm;
            M_nu_Edge = TMath::Sqrt(mu_nu_Edge2)*Ecm;
            M_tau_Max = TMath::Sqrt(mu_tau_Max2)*Ecm;
        }
        else
        {
            M_tau_Edge = M_tau_Edge = TMath::Sqrt(mu_tau2_Edge2)*Ecm;
            M_nu_Edge = 0;
        }
        
        Double_t xval = 0;
        if(mode==1) xval = M_tau_Edge;
        //if(mode==2) xval = M_tau_Edge;
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
    }
    return;
    
    
}

#include <iostream>
#include <fstream>

using namespace std;

TH1 * StrucFact(TH1 * hist,double &val, int & binMax)
{
    hist->Scale(1./hist->GetEntries());

    TH1 *hm =0;

    hm = hist->FFT(hm, "MAG R2C M");
    hm->Multiply(hm);
    if(hm->GetYaxis()->GetNbins()>1)
    {
        hm->Scale(1.0/hm->GetBinContent(1,1));
        hm->SetBinContent(1,1,0);
    }
    else
    {
        hm->Scale(1.0/hm->GetBinContent(1));
        hm->SetBinContent(1,0);
    }
    binMax=hm->GetMaximumBin();
    val=hm->GetBinContent(binMax);
    return hm;

}
TH1 * CalStruFact(TH1 * hist)
{
    TH1 * xyProj=nullptr;
    if(hist->GetZaxis()->GetNbins()>1)
    {
          xyProj=((TH3D*)hist)->Project3D("xy");
    }
    else
    {
        xyProj=(TH1 *)hist->Clone();
    }

    ofstream  theMaxXY("theMaxXY",std::ofstream::out | std::ofstream::app);
    double val;
    int binMax;
    TH1 * h1=StrucFact(xyProj,val,binMax);
    delete xyProj;
    theMaxXY<<val<<" "<<binMax<<endl;
    theMaxXY.close();

    return h1;

}

void StructFactor ()
{
    TVirtualFFT::SetTransform(0);
    TFile *MyFile = new TFile("RootFile.root","UPDATE");
    TVectorD *v=nullptr;
gROOT->SetBatch(kTRUE);
        TH1 * Ave =nullptr;


    size_t stp=0;


    v = (TVectorD*)gDirectory->Get("v");
    cout<<"v="<<(*v)[0]<<endl;
    TH1* hist = nullptr;
    bool var=1;

    TCanvas* c1 = new TCanvas("c1", "c1", 1200,1000);
    gStyle->SetOptStat(0);


    for(int i=1;i<=(((*v)[0]<1000)?((*v)[0]):(999));i++)
        {


            hist=(TH1* )gDirectory->Get(("pos" + to_string(i)).c_str());
            if(hist!=nullptr)
            {

                if(var)
                {
                    Ave=CalStruFact(hist);

                    var=0;

                }
                else
                {
                   TH1 * varr=CalStruFact(hist);
                   Ave->Add(varr);

                   delete varr;
                   delete gROOT->FindObject("out_MAG R2C M");


                }

            }
            delete hist;
            stp=i;

        }


    Ave->GetXaxis()->SetTitle("X");
    Ave->GetYaxis()->SetTitle("Y");

    Ave->GetYaxis()->CenterTitle(true);

    Ave->GetXaxis()->CenterTitle(true);
    Ave->Draw();
    cout<<"Full NEntries="<<Ave->GetEntries()<<endl;

    c1->Print("StructFactor.png");




}

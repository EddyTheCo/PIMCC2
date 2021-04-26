#include<iostream>
#include <iomanip>
#include "site.hh"
#include"constants.hh"
#include"lattice.hh"
#include"block.hh"
#include<stack>
#include <cstring>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>


const bool restart=ReadFromInput<string>(10)=="restart";
#ifdef USEROOT
#include<TFile.h>
#include<TH1.h>
#include<TH2D.h>
#include<TH3D.h>


TFile *lattice::RootFile = (restart)?new TFile("RootFile.root","UPDATE"):new TFile("RootFile.root","RECREATE");
TH1* lattice::hpos = nullptr;
TH2D* lattice::Greens = nullptr;
TH1D* lattice::PCFUp = nullptr;
#endif

using namespace std;
using namespace Constants;


const size_t SAMPLING=ReadFromInput<size_t>(20);
const size_t NPartiIni=(restart)?ReadFromInput<size_t>(1,".restart.conf"):ReadFromInput<size_t>(1);
ofstream  lattice::thesweep("sweep",(restart)?std::ofstream::out | std::ofstream::app:std::ofstream::out),lattice::theratios("ratios",(restart)?std::ofstream::out | std::ofstream::app:std::ofstream::out);
size_t Warmup=ReadFromInput<int>(21);
bool isGrandCanonical=ReadFromInput<string>(11)=="GrandCanonical";

const size_t lattice::NRep=ReadFromInput<size_t>(5);

const size_t lattice::NSweeps=ReadFromInput <size_t> (8);




array<vector<Site>,10000>*  const lattice::grid=new array<vector<Site>,10000>;
lattice::lattice()
{


    for(size_t i=0;i<NTimeSlices;i++)
    {
            grid->at(i).reserve(10000);

            for(size_t j = 0; j<NPartiIni; j++)
            {
                if(restart)
                {
                    double var;
                    vector<double> x;
                    for(size_t k=0;k<d;k++)
                    {
                        (* position::inFile)>>var;
                        x.push_back(var);

                    }





                     grid->at(i).push_back(Site(j,i,position(x),true));
                }
                else
                {


                    if(i==0)
                    {
                        vector<double> x;
                        for(size_t k=0;k<d;k++)
                        {
                            x.push_back(Constants::giveRanD(position::L.x.at(k))-position::L.x.at(k)/2);
                        }
                        grid->at(i).push_back(Site(j,i,position(x),true));
                    }
                    else
                    {
                        grid->at(i).push_back(Site(j,i,grid->at(0).at(j).pos,true));
                    }



                }

            }


    }
#ifdef USEROOT
 if (RootFile->IsZombie()) {
      cout << "Error opening file" << endl;
  }
    if(!RootFile->IsOpen())
	{
		RootFile=new TFile("RootFile.root","RECREATE");
	}
    initialize_histos();
#endif




}
#ifdef USEROOT
TVectorD *lattice::v=nullptr;
 void lattice::initialize_histos(void)
{

    if(restart)
    {
        Greens = (TH2D*)gDirectory->Get("Greens");
        PCFUp = (TH1D*)gDirectory->Get("PCFUp");


        if(!Greens)Greens=new TH2D("Greens","",5000,0,sqrt((position::L).norm()),NTimeSlices,-0.5,NTimeSlices-0.5);
        if(!PCFUp)PCFUp=new TH1D("PCFUp","",5000,0,sqrt((position::L).norm()));
        v = (TVectorD*)gDirectory->Get("v");
        if(!v)
        {
            TVectorD v(1);
            v[0]=0.;
            v.Write("v");
        }

    }
    else
    {
        Greens=new TH2D("Greens","",5000,0,sqrt((position::L).norm()),NTimeSlices,-0.5,NTimeSlices-0.5);
        PCFUp=new TH1D("PCFUp","",5000,0,sqrt((position::L).norm()));

        TVectorD v(1);
        v[0]=0.;
        v.Write("v");
    }


    v = (TVectorD*)gDirectory->Get("v");

    if(d==3)
    hpos=new TH3D("pos","",100,-position::L.x.at(0)/2,position::L.x.at(0)/2,100,-position::L.x.at(1)/2,position::L.x.at(1)/2,100,-position::L.x.at(2)/2,position::L.x.at(2)/2);
    if(d==2)
    hpos=new TH2D("pos","",500,-position::L.x.at(0)/2,position::L.x.at(0)/2,500,-position::L.x.at(1)/2,position::L.x.at(1)/2);
    if(d==1)
    hpos=new TH1D("pos","",500,-position::L.x.at(0)/2,position::L.x.at(0)/2);

}
#endif




void lattice::setup()const
{

    grid->at(0).at(0).setLattice(grid);


    for(size_t i=0;i<NTimeSlices;i++)
    {
            for(size_t j = 0; j<NPartiIni; j++)
            {
                 grid->at(i).at(j).setPtr();
            }
    }
#ifndef WARMUP
    for(size_t i=0;i<NTimeSlices;i++)
    {

            for(size_t j = 0; j<NPartiIni; j++)
            {

                 grid->at(i).at(j).totalEnergy();
            }


    }
#endif

}


void lattice::move()const
{


#ifndef WARMUP
if(!restart)thesweep<< left << setw(12) <<"KEnergy"<< left << setw(12) <<"PEnergy"<< left << setw(12) <<"TEnergy"<< left << setw(12) <<"WormLenght"<< left << setw(12) <<"SuperFlDens"<< right << setw(12) <<"NParti"<</* right << setw(12) <<"S(k,w=0)"<<*/endl;
#endif

theratios<< left << setw(12)<<"Ropen"<<left << setw(12)<<"RClose"<<left << setw(12)<<"Rmove"<<left << setw(12)<<"Rswap"<<left << setw(12)<<"RInsert"<<left << setw(12)<<"RWiggle"<<left << setw(12)<<"RShift"<<left << setw(12)<<"RRmove"<<endl;
    for(size_t step=0;step<NRep;step++)
    {

        const auto myBlock=block(NSweeps
#ifdef USEROOT
                                 ,Greens
#endif
                                 );
cout<<"finished block "<<step<<
#ifdef WARMUP
      " NPati="<<Site::getNparti()<<" "<< left << setw(12) <<getOpenRatio()<< left << setw(12) <<getCloseRatio()<< left << setw(12) <<getMoveRatio()<< left << setw(12)
   <<getSwapRatio()<< left << setw(12) <<getInsertRatio()<< left << setw(12) <<getWiggleRatio()<< left << setw(12) <<getShiftRatio()<< left << setw(12) <<getRemoveRatio()<<
#endif
      endl;
            PrintConfiguration(
#ifdef USEROOT
                        (*v)[0]+1
#endif
  );

#ifndef WARMUP
thesweep << left << setw(12) << myBlock.getKineticEnergy()<<" "<< left << setw(12) <<myBlock.SumOfPotential/(NTimeSlices)
         <<" "<< left << setw(12) << myBlock.getKineticEnergy()+myBlock.SumOfPotential/(NTimeSlices)
          <<" "<<left << setw(12)    << myBlock.Wormlenght/(NTimeSlices-1)<<" "<<left << setw(12)    << myBlock.getSuperfluidDensity()<<" "<<right << setw(12)    <<myBlock.NumberOfParticles<</*" "<<right << setw(12)    <<StructureFactor<<*/endl;



       theratios<< left << setw(12) <<getOpenRatio()<< left << setw(12) <<getCloseRatio()<< left << setw(12) <<getMoveRatio()<< left << setw(12)
                <<getSwapRatio()<< left << setw(12) <<getInsertRatio()<< left << setw(12) <<getWiggleRatio()<< left << setw(12) <<getShiftRatio()<< left << setw(12) <<getRemoveRatio()<<endl;
#endif
    }

#ifdef USEROOT
    RootFile->Close();
#endif
    thesweep.close();
    theratios.close();
}

void lattice::PrintConfiguration (
        #ifdef USEROOT
        const size_t step=0
        #else
        void
        #endif
        ) const
{

   ofstream RestartConf(".restartVAR.conf");
   ofstream RestartPtrConf(".restartPtrVAR.conf");
#ifdef SAVECONF
   static size_t step=0;
   std::ofstream Data;
     Data.open ("conf.dat", std::ofstream::out | std::ofstream::app);
   static int var=1;
   if(var)
   {
       Data<<"#Beta="<<beta<<" Nparticles="<<Site::getNparti()<<" NtimesLices="<<NTimeSlices<<endl;
       var=0;
   }
#endif
   RestartConf<<Site::getNparti()<<" #Particles"<<endl;
   RestartConf.precision(12);


   for(size_t i=0;i<NTimeSlices;i++)
   {

           for(size_t j = 0; j<Site::getNparti(); j++)
           {

               const auto var=grid->at(i).at(j);
               RestartConf<<var.pos;
#ifdef SAVECONF
               if(!(step%SAMPLING)&&step!=0)Data<<var.pos;
#endif

               RestartPtrConf<<var.left->ParticleOnBead<<" "<<var.left->TimeSliceOnBead<<" "<<var.right->ParticleOnBead<<" "<<var.right->TimeSliceOnBead<<" ";
#ifdef USEROOT

                if(d==3)
                ((TH3D*)hpos)->Fill(grid->at(i).at(j).pos.perio(0),grid->at(i).at(j).pos.perio(1),grid->at(i).at(j).pos.perio(2));
                if(d==2)
                ((TH2D*)hpos)->Fill(grid->at(i).at(j).pos.perio(0),grid->at(i).at(j).pos.perio(1));
                if(d==1)
                ((TH1D*)hpos)->Fill(grid->at(i).at(j).pos.perio(0));

                for(size_t k = 0; k<Site::getNparti(); k++)
                {
                    if(j!=k)
                    {
                        const double dis=sqrt((grid->at(i).at(j).pos-grid->at(i).at(k).pos).norm());
                        PCFUp->Fill(dis);
                    }
                }


#endif
            }
#ifdef SAVECONF
           if(!(step%SAMPLING)&&step!=0)Data<<endl;
#endif
   }
#ifdef SAVECONF
   if(!(step%SAMPLING)&&step!=0)Data<<endl;
   Data.close();
   step++;
#endif
#ifdef USEROOT

    (*v)[0] = step;



       (*v).Write("v",TObject::kOverwrite);

        hpos->Write(("pos" + to_string(step%1000)).c_str(),TObject::kOverwrite);
        hpos->Reset("ICESM");

   gDirectory->Write("", TObject::kOverwrite);

#endif


Constants::saveRandom();
   RestartConf.close();
   RestartPtrConf.close();
   rename(".restartVAR.conf", ".restart.conf");
   rename(".restartPtrVAR.conf", ".restartPtr.conf");







   

}


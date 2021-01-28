#include"block.hh"
#include<iostream>
#include<fstream>
#include"lattice.hh"
#include<vector>
#include"constants.hh"
#ifdef USEROOT
#include<TH2D.h>
#endif
using namespace std;

#include"site.hh"




block::block(array<vector<Site>,10000>* particles, const size_t &NTimeSlices, const size_t &NSweeps
             #ifdef USEROOT
             , TH2D * const Greens
             #endif
             )
{



double TSumOfdisplacement=0,TSumOfPotential=0,TNumberOfParticles=0;
size_t TWormlenght=0,step=0,measureCounter=0,measureCounter1=0;
double TWinding=0;
Site* const start=&(particles->at(0).at(0));

#ifdef WARMUP
size_t h=0;
size_t War=Warmup;
static size_t CorrectNpart=0;
#endif
while(step<NSweeps)
{



//start->printLattice();
//cout<<" Nparti="<<start->NParti_<<" Npartini="<<NPartiIni<<endl;

#ifdef WARMUP
    if(!(h%1000)&&War&&isGrandCanonical)
    {
        if(start->NParti_<War)
        {
            Site::mu+=1;
        }
        if(start->NParti_>War)
        {
            Site::mu-=1;
        }
        if(start->NParti_==War)
        {
            CorrectNpart++;
            if(CorrectNpart>100)
            {
                War=0;
                isGrandCanonical=0;
            }
        }

        cout<<"mu="<<Site::mu<<" eta="<<Site::eta<<" Npar="<<start->NParti_<<" C="<<CorrectNpart<<endl;


    }
    h++;

#endif
         if(start->ThereIsAWorm)
        {
#ifdef WARMUP
                if(!Warmup)
#endif
                {
                    measureCounter1++;
                TWormlenght+=start->NInactiveLinks();
#ifdef USEROOT
                Greens->Fill(sqrt((start->Rbead->pos-start->Lbead->pos).norm()),abs(1.*start->Rbead->TimeSliceOnBead-1.*start->Lbead->TimeSliceOnBead));
#endif
                }
            switch ((!isGrandCanonical)?giveRanI(2):giveRanI(3)) {
            case 0:
            {
               //cout<<"closing worm"<<endl;
                    start->NCloseP++;
                    if(!(start->cantClose(MBar)))
                    {
                    if(start->Lbead->CloseWorm(0))
                    {
                        step++;

                    }
                    }
                break;
            }
            case 1:
            {
               //cout<<"MoveWorm"<<endl;
                start->NMoveP++;
                    start->MoveWorm();
                     break;
            }
            case 2:
            {
               //cout<<"swap"<<endl;
                start->NSwapP++;
               if(start->NParti_>1)
               start->PrepareSwap();

               break;
            }
            case 3:
            {

                //cout<<"removeWorm"<<endl;
                start->removeWorm();
                break;
            }



            }

        }
        else
        {
             if(!Warmup)
             {
                TSumOfdisplacement+=start->TEnergy;
                TSumOfPotential+=start->TPotential;
                TNumberOfParticles+=start->NParti_;
                (d>2)?TWinding+=start->TWinding.normxy():TWinding+=start->TWinding.norm();
                measureCounter++;
             }


             switch ((isGrandCanonical)?giveRanI(3):giveRanI(2)) {
             case 0:
             {

               if(start->NParti_)
               {
                  // cout<<"OpenWorm"<<endl;
                    start->NOpenP++;
                   const size_t posiTimes=giveRanI(NTimeSlices-1) ;
                   const size_t posiParti=giveRanI(particles->at(posiTimes).size()-1);
                   const size_t var2=  giveRanI(MBar-2);
                   Site* const Ranbead=&(particles->at(posiTimes).at(posiParti));
                   start->ThereIsAWorm= Ranbead->OpenWorm(var2,var2+1,0,Ranbead->pos);
               }

                break;
             }
             case 1:
            {
                 if(start->NParti_)
                 {
                    start->NWiggleP++;
                   //cout<<"wiggle"<<endl;

                     const size_t posiTimes=giveRanI(NTimeSlices-1) ; //Choose a random time slice
                     const size_t posiParti=giveRanI(particles->at(posiTimes).size()-1); //Choose the particle

                        start->Lbead=&(particles->at(posiTimes).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)
                        const size_t var2= giveRanI(MBar-3)+1;

                        start->Rbead=start->Lbead->searchBead(true,var2);
                        start->Lbead->oldpos=start->Lbead->pos;

                        start->Lbead->Wiggle(0);

                        start->Lbead=nullptr;
                        start->Rbead=nullptr;

                }
                break;
            }
              case 3:
             {
                 //cout<<"insertworminclose "<<Constants::giveRanD(1)<<" "<<Constants::giveRanDNormal(0,1)<<endl;
                 start->insertWorm();
                 break;
             }
             case 2:
            {
                 start->NShiftP++;
                 if(start->NParti_)
                 {

//                 cout<<"ShiftParticle"<<endl;


                     const size_t posiParti=giveRanI(start->NParti_-1); //Choose the particle

                        start->Lbead=&(particles->at(0).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)


                        vector<double> varVec;

                        varVec.push_back(-position::L.TheX()/2.0+giveRanD(position::L.TheX()));
                        if(d>1)varVec.push_back(-position::L.TheY()/2.0+giveRanD(position::L.TheY()));
                        if(d>2)varVec.push_back(-position::L.TheZ()/2.0+giveRanD(position::L.TheZ()));
                        const position p=position(varVec);


                        start->Lbead->shiftParticle(0,p);

                        start->Lbead=nullptr;


                }
                break;
            }

             }


        }

    }

#ifdef WARMUP
if(!Warmup)
#endif
{
        SumofDisplacement=TSumOfdisplacement/measureCounter;
        SumOfPotential=TSumOfPotential/measureCounter;
        NumberOfParticles=TNumberOfParticles/measureCounter;
        Wormlenght=1.*TWormlenght/measureCounter1;
        SumofWinding=TWinding/measureCounter;
}
#ifdef WARMUP
if(!War&&Warmup)
{
    Warmup=0;
}
#endif
}




block::~block()
{

}

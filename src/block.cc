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




block::block(const size_t &NSweeps
             #ifdef USEROOT
             , TH2D * const Greens
             #endif
             )
{



double TSumOfdisplacement=0,TSumOfPotential=0,TNumberOfParticles=0;
size_t TWormlenght=0,step=0,measureCounter=0,measureCounter1=0;
double TWinding=0;


#ifdef WARMUP
size_t h=0;
size_t War=Warmup;
static size_t CorrectNpart=0;
#endif
while(step<NSweeps)
{





#ifdef WARMUP
    if(!(h%1000)&&War&&isGrandCanonical)
    {
        if(Site::getNparti()<War)
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
         if(Site::ThereIsAWorm)
        {
#ifdef WARMUP
                if(!Warmup)
#endif
                {
                    measureCounter1++;
                TWormlenght+=Site::NInactiveLinks();
#ifdef USEROOT
                Greens->Fill(sqrt((Site::Rbead->pos-Site::Lbead->pos).norm()),abs(1.*Site::Rbead->TimeSliceOnBead-1.*Site::Lbead->TimeSliceOnBead));
#endif
                }
            switch ((!isGrandCanonical)?giveRanI(2):giveRanI(3)) {
            case 0:
            {
               //cout<<"closing worm"<<endl;
                    Site::NCloseP++;
                    if(!(Site::cantClose(MBar)))
                    {
                    if(Site::Lbead->CloseWorm(0))
                    {
                        step++;

                    }
                    }
                break;
            }
            case 1:
            {
               //cout<<"MoveWorm"<<endl;
                Site::NMoveP++;
                    Site::MoveWorm();
                     break;
            }
            case 2:
            {
               //cout<<"swap"<<endl;
                Site::NSwapP++;
               if(Site::getNparti()>1)
               Site::PrepareSwap();

               break;
            }
            case 3:
            {

                //cout<<"removeWorm"<<endl;
                Site::removeWorm();
                break;
            }



            }

        }
        else
        {
             if(!Warmup)
             {
                TSumOfdisplacement+=Site::TEnergy;
                TSumOfPotential+=Site::TPotential;
                TNumberOfParticles+=Site::getNparti();
                (d>2)?TWinding+=Site::TWinding.normxy():TWinding+=Site::TWinding.norm();
                measureCounter++;
             }


             switch ((isGrandCanonical)?giveRanI(2):giveRanI(1)) {
             case 0:
             {

               if(Site::getNparti())
               {
                  // cout<<"OpenWorm"<<endl;
                    Site::NOpenP++;
                   const size_t posiTimes=giveRanI(NTimeSlices-1) ;
                   const size_t posiParti=giveRanI(Site::getNparti() -1);
                   const size_t var2=  giveRanI(MBar-2);
                   Site* const Ranbead=&(Site::theParticles->at(posiTimes).at(posiParti));
                   Site::ThereIsAWorm= Ranbead->OpenWorm(var2,var2+1,0,Ranbead->pos);
               }

                break;
             }
             case 1:
            {
                 if(Site::getNparti())
                 {
                    Site::NWiggleP++;
                   //cout<<"wiggle"<<endl;

                     const size_t posiTimes=giveRanI(NTimeSlices-1) ; //Choose a random time slice
                     const size_t posiParti=giveRanI(Site::getNparti()-1); //Choose the particle

                        Site::Lbead=&(Site::theParticles->at(posiTimes).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)
                        const size_t var2= giveRanI(MBar-3)+1;

                        Site::Rbead=Site::Lbead->searchBead(true,var2);
                        Site::Lbead->oldpos=Site::Lbead->pos;

                        Site::Lbead->Wiggle(0);

                        Site::Lbead=nullptr;
                        Site::Rbead=nullptr;

                }
                break;
            }
              case 2:
             {
                 //cout<<"insertworminclose "<<Constants::giveRanD(1)<<" "<<Constants::giveRanDNormal(0,1)<<endl;
                 Site::insertWorm();
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

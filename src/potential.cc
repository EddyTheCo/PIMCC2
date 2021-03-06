#include"potential.hh"
#include<iostream>
#include"position.hh"
#include"site.hh"
const position potential::k=ReadFromInput<position>(12);
const double potential::CLieb=ReadFromInput<double>(18);
const double potential::CDip=ReadFromInput<double>(18);
const double potential::V_0Softcore=ReadFromInput<double>(18);
const double potential::SigmaBonin=ReadFromInput<double>(18);
const bool potential::harmonic=ReadFromInput<string>(16)=="harmonic";
const bool potential::freeExt=ReadFromInput<string>(16)=="free";
const bool potential::freeInt=ReadFromInput<string>(17)=="free";
const bool potential::infinteWell=ReadFromInput<string>(16)=="infinteWell";
const bool potential::LiebLin=ReadFromInput<string>(17)=="LiebLiniger";
const bool potential::Dipola=ReadFromInput<string>(17)=="Dipolar";
const bool potential::Softcor=ReadFromInput<string>(17)=="Softcore";
const bool potential::Bonin=ReadFromInput<string>(17)=="Bonin";
using namespace std;


void potential::U(const Site* const bead,double& dU,position& graddU ) const
{
    if(!freeExt)ExternPot(bead,dU,graddU);
    if(!freeInt)PairInteraction(bead,dU,graddU);
}

void potential::ExternPot(const Site * const bead, double& dU, position & graddU)const
{

if(harmonic)
{
    position var=position(0.);
    for(size_t i=0;i<d;i++)
    {
        dU+=k.x.at(i)*bead->pos.x.at(i)*bead->pos.x.at(i)/2;
        var.x.at(i)=k.x.at(i)*bead->pos.x.at(i);
    }
    graddU=graddU+var;
    return;
}
if(infinteWell)
{
    if(bead->pos>position::L/2.)
    {
        dU+= 1000000000000000;
        graddU=graddU+position(10000000000000000.0);
    }
    return;
}
//cout<<"##ERROR error in the potential"<<endl;
}

void potential::PairInteraction(const Site *const bead,double& dU,position & graddU)const
{

   const Site * ptr=bead->up;
while(ptr!=bead)
{

   if(ptr->active&&ptr!=ptr->theZeta)
   {
       if(Dipola)
       {
           Dipolar(dU,graddU,bead,ptr);
           ptr=ptr->up;
           continue;
       }
       if(LiebLin)
       {
           //cout<<"lieb"<<endl;
           LiebLini(dU,graddU,bead,ptr);
           ptr=ptr->up;
           continue;
       }
       if(Softcor)
       {
           //cout<<"Softcore"<<endl;
           Softcore(dU,graddU,bead,ptr);
           ptr=ptr->up;
           continue;
       }
       if(Bonin)
       {
           BoninPot(dU,graddU,bead,ptr);
           ptr=ptr->up;
           continue;
       }

   }

   ptr=ptr->up;

}
return;
}
inline void potential::LiebLini(double & dU,position& graddU,const Site* const bead,const Site* const ptr)const
{

    if(sqrt((bead->pos-ptr->pos).norm())<=1.)
    {
        dU+= CLieb;
    }

}
inline void potential::Dipolar(double & dU, position& graddU, const Site * const bead, const Site * const ptr)const
{

   dU+= 3./pow((bead->pos-ptr->pos).norm(),1.5)+CDip/pow((bead->pos-ptr->pos).norm(),6);

   graddU=graddU+(bead->pos-ptr->pos)*(-9./pow((bead->pos-ptr->pos).norm(),2.5)) -(bead->pos-ptr->pos)*(12.*CDip/pow((bead->pos-ptr->pos).norm(),7));
}
inline void potential::Softcore(double & dU, position& graddU, const Site * const bead, const Site * const ptr)const
{
    const position dif=bead->pos-ptr->pos;
    const double difNorm=dif.norm();
   dU+= V_0Softcore/(1.+pow(difNorm,3.));
   graddU=graddU+(dif)*(-6.*V_0Softcore*pow(difNorm,2.)/(pow(1.+pow(difNorm,3.),2.)));
}
inline void potential::BoninPot(double &dU, position &graddU, const Site *const bead, const Site *const ptr) const
{
    const position dif=bead->pos-ptr->pos;
    const double difNorm=dif.norm();
dU+= (1.-3.*(dif.x.at(2)*dif.x.at(2))/(difNorm))/(pow(difNorm,1.5))+pow(SigmaBonin*SigmaBonin/difNorm,6.);

   vector<double> UsrGrad;
   vector<double> varVec;
   varVec.push_back(0.);
   varVec.push_back(0.);
   varVec.push_back(-6*dif.x.at(2)/pow(difNorm,2.5));
   position var(varVec);

   graddU=graddU+ dif*(-12*pow(SigmaBonin*SigmaBonin/difNorm,6.))/difNorm + dif*(-3)/pow(difNorm,2.5) + dif*dif.x.at(2)*dif.x.at(2)*15/pow(difNorm,3.5) + var;
}

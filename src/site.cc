#include<iostream>
#include "site.hh"
#include"constants.hh"
#include"position.hh"
#include"potential.hh"
#include"lattice.hh"
#include<stack>
#include<map>
using namespace std;
using namespace Constants;


const bool Site::fourthOrder=ReadFromInput<string>(19)=="true";
size_t Site::NpartixNT=0;

const double landa=ReadFromInput<double>(4);
const double tao=ReadFromInput<double>(2);
const double beta=ReadFromInput<double>(3);

double Site::mu=ReadFromInput<double>(9);
double Site::eta=ReadFromInput<double>(13);



const double variance=2*landa*tao;
ifstream Site::RestartPtrConf(".restartPtr.conf");
Site* Site::theZeta=nullptr;
const size_t NTimeSlices=beta/tao;
const size_t MBar=(ReadFromInput<size_t>(15)<1||ReadFromInput<size_t>(15)>NTimeSlices)?NTimeSlices:ReadFromInput<size_t>(15);

#ifndef WARMUP
double Site::TEnergy=0,Site::TPotential=0;
double Site::TEnergyVar=0,Site::TPotentialVar=0;

position Site::TWinding=position(0.);
position Site::TWindingVar=position(0.);
#endif

const potential Site::ThePotential=potential();


bool Site::ThereIsAWorm=false;
array<vector<Site>,10000>* Site::theParticles=nullptr;

Site* Site::Lbead=nullptr;
Site* Site::Rbead=nullptr;

size_t Site::NClose=1,Site::NOpen=1,Site::NMove=1,Site::NWiggle=1,Site::NWiggleP=1,Site::NShift=1,Site::NShiftP=1,Site::NSwap=1,Site::NCloseP=1,Site::NOpenP=1,Site::NMoveP=1,Site::NSwapP=1,Site::NInsertP=1,Site::NInsert=1,Site::NRemoP=1,Site::NRemo=1;

Site::Site(const size_t i, const size_t j, const position var, const bool act):active(act),pos(var),oldpos(pos),TimeSliceOnBead(j),ParticleOnBead(i),left(nullptr),right(nullptr),up(nullptr),down(nullptr)
{


NpartixNT++;
}



Site::~Site(){
    NpartixNT--;
}
bool Site::OpenWorm(const size_t step, const size_t ab, double dU, const position& start )
{

    dU+=-mu*tao;
    if(step>0)
    {
        double U=0;
        right->ChangeInU(true,dU,U);

        if(right->OpenWorm(step-1,ab,dU,start))
        {
            Lbead=this;
            right->active=false;
#ifndef WARMUP
            const auto Dist=pos-right->pos;
            TEnergy-=Dist.norm();
            TWinding=TWinding+Dist;
            TPotential-=U;
#endif
            return true;
        }
        return false;

    }
    else
    {


        if(eta*getNparti()/(propagator(start,right->pos,ab,-dU)*position::volumen)>giveRanD(1.))
        {
#ifndef WARMUP
            const auto Dist=pos-right->pos;
            TEnergy-=Dist.norm();
            TWinding=TWinding+Dist;
#endif
            Rbead=this->right;
            ThereIsAWorm=true;
            Lbead=this;
            NOpen++;
            return true;
        }
        return false;
    }
}
bool Site::CloseWorm(double dU)
{



    dU+=mu*tao;

    if(this->TimeSliceOnBead!=Rbead->left->TimeSliceOnBead)
    {

        right->pos=position(true,this);
        double U=0;
        right->ChangeInU(false,dU,U);


        if(right->CloseWorm(dU))
        {
            right->active=true;
#ifndef WARMUP
            const auto Dist=right->pos-pos;
            TEnergy+=Dist.norm();
            TWinding=TWinding+Dist;
            TPotential+=U;
#endif
            return true;

        }
        return false;

    }
    else
    {
        if(propagator(Lbead->pos,Rbead->pos,NInactiveLinks(),dU)*position::getVolumen()/(eta*getNparti())>giveRanD(1.))
        {
#ifndef WARMUP
            const auto Dist=right->pos-pos;
            TEnergy+=Dist.norm();
            TWinding=TWinding+Dist;
#endif
            Rbead=nullptr;
            Lbead=nullptr;
            ThereIsAWorm=false;
            NClose++;
            return true;
        }
        return false;
    }
}
bool Site::Wiggle(double dU)
{

right->oldpos=right->pos;

    if(this->TimeSliceOnBead!=Rbead->left->TimeSliceOnBead)
    {
#ifndef WARMUP
        const auto Dist=oldpos-right->pos;
        TEnergyVar-=Dist.norm();
        TWindingVar=TWindingVar+Dist;
#endif
        double U=0;

        right->ChangeInU(true,dU,U);
#ifndef WARMUP
        TPotentialVar+=U;
#endif

        right->pos=position(true,this);


        right->ChangeInU(false,dU,U);

        if(right->Wiggle(dU))
        {
#ifndef WARMUP
            const auto Dist=right->pos-pos;
            TEnergy+=Dist.norm();
            TWinding=TWinding+Dist;
            TPotential+=U;
#endif

            return true;
        }
        right->pos=right->oldpos;
        return false;
    }
    else
    {
#ifndef WARMUP
        const auto Dist=oldpos-right->pos;
        TEnergyVar-=Dist.norm();
        TWindingVar=TWindingVar+Dist;
#endif

        if(exp(dU)>giveRanD(1.))
        {
            NWiggle++;
#ifndef WARMUP
            TPotential-=TPotentialVar;
            const auto Dist=right->pos-pos;
            TEnergy+=TEnergyVar+Dist.norm();
            TWinding=TWinding+TWindingVar+Dist;
#endif
            Lbead=nullptr;
            Rbead=nullptr;
#ifndef WARMUP
            TPotentialVar=0;
            TEnergyVar=0;
            TWindingVar=position(0);
#endif
            return true;
        }
#ifndef WARMUP
        TPotentialVar=0;
        TEnergyVar=0;
        TWindingVar=position(0);
#endif
        return false;
    }
}



bool Site::deleteToRight(const size_t step, double dU)
{

dU+=-mu*tao;
double U=0;
ChangeInU(true,dU,U);


        if(step>0)
        {
            if(right->deleteToRight(step-1,dU))
            {
                active=false;
#ifndef WARMUP
                const auto Dist=pos-right->pos;
                TEnergy-=Dist.norm();
                TWinding=TWinding+Dist;
                TPotential-=U;
#endif
                return true;
            }
            return false;

        }
        else
        {

            if(this==Lbead)dU+=mu*tao;


            if(exp(dU)>giveRanD(1.))
            {
#ifndef WARMUP
                if(this!=Lbead)
                {
                    const auto Dist=pos-right->pos;
                    TEnergy-=Dist.norm();
                    TWinding=TWinding+Dist;
                }
                TPotential-=U;
#endif
                active=false;
                Rbead=this->right;
                NMove++;
                return true;
            }
            return false;
        }


}
bool Site::deleteToLeft(const size_t step,double dU)
{

dU+=-mu*tao;
double U=0;
ChangeInU(true,dU,U);
            if(step>0)
            {
                if(left->deleteToLeft(step-1,dU))
                {
                    active=false;
#ifndef WARMUP
                    const auto Dis=left->pos-pos;
                    TEnergy-=Dis.norm();
                    TWinding=TWinding+Dis;
                    TPotential-=U;
#endif
                    return true;
                }
                return false;

            }
            else
            {

                if(exp(dU)>giveRanD(1.))
                {
                    active=false;
#ifndef WARMUP
                    TPotential-=U;
                    const auto Dis=left->pos-pos;
                    TEnergy-=Dis.norm();
                    TWinding=TWinding+Dis;
#endif
                    Lbead=this->left;
                    NMove++;
                    return true;
                }
                return false;
            }


}
bool Site::insertToRight(const size_t step,double dU)
{
    static bool aParticleisInserted=false;
    if (right->active)
    {
        if(!isGrandCanonical)
            return false;
        else {
            insertParticle();
            aParticleisInserted=true;
            right=theParticles->at(TimeSliceOnBead).back().right;
            auto var=right->left;
            right->left=this;
            Rbead->left=var;
            var->right=Rbead;
        }
    }

dU+=mu*tao;
right->pos=position(pos,variance);
double U=0;
right->ChangeInU(false,dU,U);


            if(step>0)
            {
                if(right->insertToRight(step-1,dU))
                {
                    right->active=true;
#ifndef WARMUP
                    const auto Dis=right->pos-pos;
                    TEnergy+=Dis.norm();
                    TWinding=TWinding+Dis;
                    TPotential+=U;
#endif
                    return true;
                }
                return false;

            }
            else
            {


                if(exp(dU)>giveRanD(1.))
                {
                    right->active=true;
#ifndef WARMUP
                    const auto Dis=right->pos-pos;
                    TEnergy+=Dis.norm();
                    TWinding=TWinding+Dis;
                    TPotential+=U;
#endif
                    Lbead=this->right;
                    if(aParticleisInserted)
                    {

                        aParticleisInserted=false;
                    }
                    NMove++;
                    return true;
                }

                if(aParticleisInserted)//but the move is not accepted
                {
                    auto var=&(theParticles->at(Rbead->TimeSliceOnBead).back());
                    auto varLeft=var->left;
                    var->left=Rbead->left;
                    var->left->right=var;
                    Rbead->left=varLeft;
                    varLeft->right=Rbead;
                    removeLastParticle();
                    aParticleisInserted=false;
                }
                return false;

            }

}
bool Site::insertToLeft(const size_t step,double dU)
{


    static bool aParticleisInserted=false;

    if (left->active)
    {
        if(!isGrandCanonical)
            return false;
        else {

            insertParticle();

            aParticleisInserted=true;
            left=theParticles->at(TimeSliceOnBead).back().left;
            auto var=left->right;
            left->right=this;
            Lbead->right=var;
            var->left=Lbead;
        }
    }

dU+=mu*tao;
    left->pos=position(pos,variance);
    double U=0;
    left->ChangeInU(false,dU,U);
            if(step>0)
            {

                if(left->insertToLeft(step-1,dU))
                {
                    left->active=true;
#ifndef WARMUP
                    const auto Dis=pos-left->pos;
                    TEnergy+=Dis.norm();
                    TWinding=TWinding+Dis;
                    TPotential+=U;
#endif
                    return true;
                }
                return false;

            }
            else
            {


                if(exp(dU)>giveRanD(1.))
                {
                    left->active=true;
#ifndef WARMUP
                    const auto Dis=pos-left->pos;
                    TEnergy+=Dis.norm();
                    TWinding=TWinding+Dis;
                    TPotential+=U;
#endif
                    Rbead=this->left;
                    if(aParticleisInserted)
                    {

                        aParticleisInserted=false;
                    }

                    NMove++;
                    return true;
                }
                if(aParticleisInserted)//but the move is not accepted
                {
                    auto var=&(theParticles->at(Lbead->TimeSliceOnBead).back());
                    auto varRight=var->right;
                    var->right=Lbead->right;
                    var->right->left=var;
                    Lbead->right=varRight;
                    varRight->left=Lbead;
                    removeLastParticle();
                    aParticleisInserted=false;
                }
                return false;
            }
}



Site* Site::chooseTheBead(double &SumI, const size_t& vae, const Site * const startBead)const
{
     multimap <double, Site*> prop;

    auto ptr=this->up;
    double var;
    while(ptr!=this)
    {


        if(ptr->active&&ptr!=Rbead&&ptr!=Lbead)
        {
            var=propagator(startBead->pos,ptr->pos,vae);
            if(var>0){
                SumI+=var;
                prop.insert(pair <double, Site*> (var, ptr));
            }

        }
        ptr=ptr->up;
    }
    if(prop.size()==0)
        return nullptr;
    auto ran=giveRanD(1);

    double pr=0;
    map<double, Site*>::reverse_iterator rit;
      for (rit=prop.rbegin(); rit!=prop.rend(); ++rit)
        {
            pr+=rit->first/SumI;
            if(pr>ran)
            {
                return rit->second;
            }

        }
}

bool Site::swap(Site* const zeta, const double& SumI, const double& SumZ, double dU, const bool &isRight)
{

if(isRight)
{
    static bool aParticleisInserted=false;
    static Site* ri;
    if(this->TimeSliceOnBead!=Rbead->left->TimeSliceOnBead)
    {
        //cout<<"pos="<<pos<<endl;
        if (right->active)
        {
           if(!isGrandCanonical)
           {
               theZeta=nullptr;
                return false;
           }
            else {
                insertParticle();
                aParticleisInserted=true;
                //cout<<"seinserto aparticle"<<endl;
                ri=right;
                right=theParticles->at(TimeSliceOnBead).back().right;
                auto var=right->left;
                right->left=this;
                ri->left=var;
                var->right=ri;
            }
        }

        right->pos=position(true,this);
        double Ualpha=0,Uzeta=0;
        theZeta=zeta->right;
        right->ChangeInU(false,dU,Ualpha);
        zeta->right->ChangeInU(true,dU,Uzeta);

        if(right->swap(zeta->right,SumI,SumZ,dU,true))
        {
            this->right->active=true;
            zeta->right->active=false;
#ifndef WARMUP
            const auto Dist1=right->pos-pos;
            const auto Dist2=zeta->pos-zeta->right->pos;
            TEnergy+=Dist1.norm()-Dist2.norm();
            TWinding=TWinding+Dist1+Dist2;
            TPotential+=Ualpha;
            TPotential-=Uzeta;
#endif
            return true;
        }

        return false;

    }
    else
    {

         if(exp(dU)*SumI/SumZ>giveRanD(1.))
         {
             Site* const prev=this->right;
             this->right=zeta->right;
#ifndef WARMUP
             const auto Dist1=right->pos-pos;
             const auto Dist2=zeta->pos-zeta->right->pos;
             TEnergy+=Dist1.norm()-Dist2.norm();
             TWinding=TWinding+Dist1+Dist2;
#endif
            zeta->right=prev;
            zeta->right->left=zeta;
            this->right->left=this;
            NSwap++;
            if(aParticleisInserted)
            {

                aParticleisInserted=false;

            }
            theZeta=nullptr;
            return true;
         }
         if(aParticleisInserted)//but the move is not accepted
         {

             const auto var=&(theParticles->at(ri->TimeSliceOnBead).back());
             auto varLeft=var->left;
             var->left=ri->left;
             var->left->right=var;
             ri->left=varLeft;
             varLeft->right=ri;
             removeLastParticle();
             aParticleisInserted=false;

         }
         theZeta=nullptr;
         return false;
    }

}
else {
    static bool aParticleisInserted=false;
    static Site* le;
    if(this->TimeSliceOnBead!=Lbead->right->TimeSliceOnBead)
    {
        //cout<<"pos="<<pos<<endl;
        if (left->active)
        {
            if(!isGrandCanonical)
            {
                theZeta=nullptr;
                return false;

            }
            else {
                insertParticle();
                aParticleisInserted=true;
                //cout<<"seinserto aparticle"<<endl;
                le=left;
                left=theParticles->at(TimeSliceOnBead).back().left;
                auto var=left->right;
                left->right=this;
                le->right=var;
                var->left=le;

            }
        }

        left->pos=position(false,this);
        double Ualpha=0,Uzeta=0;
        theZeta=zeta->left;
        left->ChangeInU(false,dU,Ualpha);
        zeta->left->ChangeInU(true,dU,Uzeta);

        if(left->swap(zeta->left,SumI,SumZ,dU,false))
        {
            this->left->active=true;
            zeta->left->active=false;
#ifndef WARMUP
            const auto Dist1=pos-left->pos;
            const auto Dist2=zeta->left->pos-zeta->pos;
            TEnergy+=Dist1.norm()-Dist2.norm();
            TWinding=TWinding+Dist1+Dist2;
            TPotential+=Ualpha;
            TPotential-=Uzeta;
#endif
            return true;
        }

        return false;
    }
    else
    {


         if(exp(dU)*SumI/SumZ>giveRanD(1.))
        {
             Site * const prev=this->left;
             this->left=zeta->left;
#ifndef WARMUP
             const auto Dist1=pos-left->pos;
             const auto Dist2=zeta->left->pos-zeta->pos;
             TEnergy+=Dist1.norm()-Dist2.norm();
             TWinding=TWinding+Dist1+Dist2;
#endif
            zeta->left=prev;
            zeta->left->right=zeta;
            this->left->right=this;
            NSwap++;
            if(aParticleisInserted)
            {

                aParticleisInserted=false;


            }
            theZeta=nullptr;
            return true;
        }
         if(aParticleisInserted)//but the move is not accepted
         {

             const auto var=&(theParticles->at(le->TimeSliceOnBead).back());
             auto varRight=var->right;
             var->right=le->right;
             var->right->left=var;
             le->right=varRight;
             varRight->left=le;
             removeLastParticle();
             aParticleisInserted=false;

         }
         theZeta=nullptr;
        return false;
    }

}



}


bool Site::shiftParticle(double dU, const position& shift)const
{
    right->oldpos=right->pos;

    double U=0;
    right->ChangeInU(true,dU,U);
#ifndef WARMUP
    TPotentialVar+=U;
#endif
    //cout<<"right="<<right->ParticleOnBead<<" "<<right->TimeSliceOnBead<<" "<<right->pos<<endl;
    right->pos=right->pos+shift;
    //cout<<"right="<<right->ParticleOnBead<<" "<<right->TimeSliceOnBead<<" "<<right->pos<<endl;
    right->ChangeInU(false,dU,U);

        if(this!=Lbead->left)
        {

            if(right->shiftParticle(dU,shift))
            {
#ifndef WARMUP
                TPotential+=U;
#endif
                //cout<<oldpos-right->oldpos<<" "<<pos-right->pos<<endl;
                return true;
            }
            right->pos=right->oldpos;
            return false;
        }
        else
        {


            if(exp(dU)>giveRanD(1.))
            {
                NShift++;
#ifndef WARMUP
                TPotential+=U;
                TPotential-=TPotentialVar;
                TPotentialVar=0.;
#endif
                //cout<<oldpos-right->oldpos<<" "<<pos-right->pos<<endl;
               return true;
            }
            right->pos=right->oldpos;
#ifndef WARMUP
            TPotentialVar=0.;
#endif
            return false;
        }


}

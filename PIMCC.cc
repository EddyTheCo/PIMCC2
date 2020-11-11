#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include<iostream>
#include"lattice.hh"



using namespace std;

void printinput(void);
int main()
{
auto start = chrono::high_resolution_clock::now();

if(ReadFromInput<string>(10)=="restart")    
Constants::readRandom();
else {
int va=system("cp input .input.start"); //Makes a copy of the input file
}
int va=system("cp input .input.ini"); //Makes a copy of the input file

   auto theLattice=lattice();

    theLattice.setup();


if(size_t Warmup=ReadFromInput<int>(21))
{
    cout<<"Starting The WarmingUp"<<endl;
    theLattice.Warm();
    cout<<"Finishing The WarmingUp"<<endl;
}
    theLattice.move();

auto finish = std::chrono::high_resolution_clock::now();
chrono::duration<double> elapsed = finish - start;
cout<<" Elapsed time: " << elapsed.count()<<"s"<<endl;

}



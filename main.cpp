//
//  main.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <ctime>
#include "Matrix.h"
#include "EulerEquations.h"
// Primitive
double rho(double x){
    if (x < 0.5) {
        return 1;
    }
    else {
        return  0.125;
    }
}
double p(double x){
    if (x < 0.5) {
        return 1;
    }
    else {
        return 0.1;
    }
}
double u(double x){
    return 0;
}

double rho123(double x){
    return 1;
}
double p123(double x){
    return 0.4;
}
double u123(double x) {
    if (x < 0.5) {
        return -2;
    }
    else {
        return 2;
    }
}

void WriteDensity(Profiles& u, int nCells, std::ofstream& myfile){

    for (int i = 1 ; i <= nCells; i++) {
        myfile << u.u1(i) << " ";
    }
}
void WriteVelocity(Profiles& u, int nCells, std::ofstream& myfile){
    
    for (int i = 1 ; i <= nCells; i++) {
        myfile << u.u2(i) << " ";
    }
}
void WritePressure(Profiles& u, int nCells, std::ofstream& myfile){
    for (int i = 1 ; i <= nCells; i++) {
        myfile << u.u3(i) << " ";
    }
}

void PrintDensity(Profiles& u, int nCells){
    for (int i = 1 ; i <= nCells; i++) {
        std::cout<< u.u1(i) << ", ";
    }
}
void PrintVelocity(Profiles& u, int nCells){
    for (int i = 1 ; i <= nCells; i++) {
        std::cout<< u.u2(i) << ", ";
    }
}
void PrintPressure(Profiles& u, int nCells){
    for (int i = 1 ; i <= nCells; i++) {
        std::cout<< u.u3(i) << ", ";
    }
}

int main(int argc, const char * argv[]) {
    EulerSolver sol(rho, p, u);
//    EulerSolver sol(rho123,p123,u123);
    int nCells = 1000;
    sol.SetCellNumber(nCells);
    sol.SetGamma(1.4);
    sol.SetRange(0, 1);
    sol.SetTime(0, 0.25);
    std::clock_t start;
    double duration;
    
    start = std::clock();
    
    Profiles res(nCells);
//    sol.FORCESolve(res);
    
    sol.LFSolve(res);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration "<< duration <<'\n';
    
    std::ofstream density;
    std::ofstream velocity;
    std::ofstream pressure;
    density.open("/Users/lixr/Documents/Codes/Research/EulerEquations/EulerEquations/FORCEdnew.dat");
    velocity.open("/Users/lixr/Documents/Codes/Research/EulerEquations/EulerEquations/FORCEvnew.dat");
    pressure.open("/Users/lixr/Documents/Codes/Research/EulerEquations/EulerEquations/FORCEpnew.dat");
    WriteDensity(res, nCells,density);
    WriteVelocity(res, nCells, velocity);
    WritePressure(res, nCells, pressure);
    return 0;
}

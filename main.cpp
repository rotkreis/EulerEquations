//
//  main.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include "Matrix.h"
#include "Vector.h"
#include "EulerEquations.h"
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
template<class T>
void WriteDensity(Profiles<T>& u, int nCells, std::ofstream& myfile){
    
    for (int i = 1 ; i <= nCells; i++) {
        myfile << u.u1(i) << " ";
    }
}

template<class T>
void WriteVelocity(Profiles<T>& u, int nCells, std::ofstream& myfile){
    
    for (int i = 1 ; i <= nCells; i++) {
        myfile << u.u2(i) << " ";
    }
}

template<class T>
void WritePressure(Profiles<T>& u, int nCells, std::ofstream& myfile){
    for (int i = 1 ; i <= nCells; i++) {
        myfile << u.u3(i) << " ";
    }
}


template<class T>
void PrintDensity(Profiles<T>& u, int nCells){
    for (int i = 1 ; i <= nCells; i++) {
        std::cout<< u.u1(i) << ", ";
    }
}
template<class T>
void PrintVelocity(Profiles<T>& u, int nCells){
    for (int i = 1 ; i <= nCells; i++) {
        std::cout<< u.u2(i) << ", ";
    }
}
template<class T>
void PrintPressure(Profiles<T>& u, int nCells){
    for (int i = 1 ; i <= nCells; i++) {
        std::cout<< u.u3(i) << ", ";
    }
}

int main(int argc, const char * argv[]) {
    EulerSolver<double> sol(rho, p, u);
    //    EulerSolver sol(rho123,p123,u123);
    int nCells = 200;
    sol.SetCellNumber(nCells);
    sol.SetGamma(1.4);
    sol.SetRange(0, 1);
    sol.SetTime(0, 0.25);
    std::clock_t start;
    double duration;
    Profiles<double> res(nCells);
    sol.LFSolve(res);
    
}

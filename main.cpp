//
//  main.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include <fstream>
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

void PrintDensity(Profiles& u, int nCells){

    for (int i = 1 ; i <= nCells; i++) {
        std::cout << u.u1(i) << ", ";
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
    int nCells = 1000;
    sol.SetCellNumber(nCells);
    sol.SetGamma(1.4);
    sol.SetRange(0, 1);
    sol.SetTime(0, 0.25);

    Profiles res(nCells);
    sol.HLLSolve(res);
    PrintVelocity(res, nCells);
    std::cout << std::endl;
    PrintPressure(res, nCells);

    return 0;
}

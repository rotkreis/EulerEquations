//
//  main.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
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



int main(int argc, const char * argv[]) {
    std::vector<mVector> nihao(3);
    nihao[1] = mVector(3);
    nihao[1][1] = 1.1;
    std::cout << nihao[1] << endl;
    
    EulerSolver sol(rho, p, u);
    sol.SetCellNumber(11);
    sol.SetGamma(1.4);
    sol.SetRange(0, 1);
    sol.SetTime(0, 0.25);
    sol.test();
    
    
    return 0;
}

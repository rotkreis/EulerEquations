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
int main(int argc, const char * argv[]) {
    // insert code here...
    Matrix<double> mat(3,3);
    mat(1,1)=1;
    std::cout << mat;
    mVector<double> vec(3);
    vec[1] = 1;
    std::cout << vec.NormInf();
    std::cout << "Hello, World!\n";
    return 0;
}

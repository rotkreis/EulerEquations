//
//  main.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include <iostream>
#include "Matrix.h"
int main(int argc, const char * argv[]) {
    // insert code here...
    Matrix<double> hei(3,3);
    hei(2,1) = (double)1;
    std::cout << hei;
    std::cout << "Hello, World!\n";
    return 0;
}

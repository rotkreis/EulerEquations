//
//  Vector.h
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#ifndef __EulerEquations__Vector__
#define __EulerEquations__Vector__

#include <stdio.h>
#include "Matrix.h"



class mVector : public Matrix {
public:
    using Matrix::operator=;
public:
    int dim(){
        return this -> mNumRows;
    }
    int dim() const{
        return this -> mNumRows;
    }
public:
    TYPE& operator[](int index){
        return data[0][index];
    }
    const TYPE& operator[](int index) const{
        return data[0][index];
    }
public:
    mVector() {}
    mVector(int dim){
        mNumRows = dim;
        mNumCols = 1;
        data = new TYPE*[1];
        data[0] = new TYPE[dim];
        for (int i = 0; i != dim; i++) {
            data[0][i] = 0;
        }
    }
    mVector(const mVector& rhs){
        mNumRows = rhs.dim();
        mNumCols = 1;
        data = new TYPE*[1];
        data[0] = new TYPE[rhs.dim()];
        for (int i = 0; i != rhs.dim(); i++) {
            data[0][i] = rhs[i];
        }
    }
    double NormInf();
    double NormInf(int& position);
    double Norm_2();
    double Norm_1();
    
};

double InnerProduct(mVector& v1, mVector& v2);

#endif /* defined(__EulerEquations__Vector__) */

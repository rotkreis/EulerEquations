//
//  Vector.cpp
//  EulerEquations
//
//  Created by Li Xinrui on 11/22/14.
//  Copyright (c) 2014 Li Xinrui. All rights reserved.
//

#include "Vector.h"

double mVector::NormInf(){
    double max = 0;
    for (int i = 0; i != this->mNumRows; i++) {
        if (std::abs((double) this->data[0][i]) > max ) {
            max = std::abs((double) this->data[0][i]);
        }
    }
    return max;
};

double mVector::NormInf(int& position){
    double max = 0;
    for (int i = 0; i != this->mNumRows; i++) {
        if (std::abs((double) this->data[0][i]) > max ) {
            max = std::abs((double) this->data[0][i]);
            position = i;
        }
    }
    return max;
};

double InnerProduct(mVector &v1, mVector &v2){
    double product = 0;
    assert(v1.dim() == v2.dim());
    for (int i = 0; i != v1.dim(); i++) {
        product += v1[i] * v2[i];
    }
    return product;
}

double mVector::Norm_2(){
    double sum = 0;
    for (int i = 0; i != this -> mNumRows; i++) {
        sum += this -> data[0][i] * data[0][i];
    }
    return sqrt(sum);
}

double mVector::Norm_1(){
    double sum = 0;
    for (int i = 0; i != this -> mNumRows; i++) {
        sum += std::abs(data[0][i]);
    }
    return sum;
}

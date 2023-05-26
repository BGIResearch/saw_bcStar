/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>

using namespace std;

class Sequence{
public:
    Sequence();
    Sequence(string seq);
    void print();
    int length();
    Sequence reverseComplement();

    Sequence operator~();

    static bool test();

public:
    string mStr;
};

#endif
/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H


#include "common.h"
#include "IncludeDefine.h"
#include <iostream>
#include <string.h>
class BloomFilter {
public:
    BloomFilter();
    bool push(uint64 key);
    bool get(uint64 key);
    bool push_mod(uint64 key);
    bool get_mod(uint64 key);

    bool push_xor(uint64 key);
    bool get_xor(uint64 key);

    bool push_Classification(uint64 key);
    bool get_Classification(uint64 key);

    bool push_wang(uint64 key);
    bool get_wang(uint64 key);

public:

    uint64* hashtable;
    uint64* hashtableClassification;
    uint64 size;


    const static uint32 HashTableMax = 1ll<<26;



    const uint64 Bloom_MOD = 73939133;

};

#endif //BLOOMFILTER_H

/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef CHIPMASKHDF5_H
#define CHIPMASKHDF5_H

#include <iostream>
#include <hdf5.h>
#include <unordered_map>
#include <folly/container/F14Map.h>
#include "slideMask.h"
#include "util.h"
#include "omp.h"
#include <libdeflate.h>
#include <unistd.h>
//#include "robin_hood.h"
#include "bloomFilter.h"

#define RANK 3
#define CDIM0 1537
#define CDIM1 1537
#define DATASETNAME "bpMatrix_"
#define ATTRIBUTEDIM 6   //rows, cols, rowShift, colShift, barcodeLength, slide pitch
#define ATTRIBUTEDIM1 2
#define ATTRIBUTENAME "dnbInfo"
#define ATTRIBUTENAME1 "chipInfo"

using namespace std;

class ChipMaskHDF5{
public:
    ChipMaskHDF5(std::string FileName);
    ~ChipMaskHDF5();

    void creatFile();
    herr_t writeDataSet(std::string chipID, slideRange& sliderange, unordered_map<uint64_t, Position1>& bpMap, int BarcodeLen, int segment, int slidePitch, uint compressionLevel = 6, int index = 1);
    void openFile();
    void readDataSet(folly::F14ValueMap<uint64_t, Position1>& bpMap, int index = 1);
    int readDataSetParallelize(folly::F14ValueMap<uint64_t, Position1>& bpMap,BloomFilter *&bloomFilter, int index = 1);
//    int readDataSetParallelize(folly::F14ValueMap<uint64_t, uint64_t>& bpMap,BloomFilter *&bloomFilter, int index = 1);
    void readDataSetParallelize(unordered_map<uint64_t, Position1>& bpMap,BloomFilter *&bloomFilter, int index = 1);

public:
    std::string fileName;
    hid_t fileID;
    uint64_t*** bpMatrix;
};

#endif 
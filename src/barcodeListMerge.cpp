/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "barcodeListMerge.h"
//#include "heatMap.h"

BarcodeListMerge::BarcodeListMerge(Options* opt){
    mOptions = opt;
    split(opt->in, dnbMapFiles, ",");
    mergedDnbMapFile = opt->out;
}

BarcodeListMerge::~BarcodeListMerge(){
    mergedDnbMap.clear();
    unordered_map<uint64_t, int>().swap(mergedDnbMap);
}

void BarcodeListMerge::mergeBarcodeLists(){
    unordered_map<uint64_t, int> dnbMap;
    ifstream inDnb;
    for (auto dnbFile = dnbMapFiles.begin(); dnbFile!=dnbMapFiles.end(); dnbFile++){
        dnbMap.clear();
        if (ends_with(*dnbFile, ".bin")){
            inDnb.open(dnbFile->c_str(), ios::binary);
            boost::archive::binary_iarchive ia(inDnb);
            ia >> dnbMap;
            inDnb.close();
            addBarcodeList(dnbMap);
        }else{
            inDnb.open(dnbFile->c_str());
            int x;
            int y;
            int count;
            uint64_t encodePos;
            unordered_map<uint64_t, int>::iterator mergedDnbIter;
            while(inDnb >> x >> y >> count){
                encodePos = ((uint64_t)x << 32) | (uint64_t)y;
                mergedDnbIter = mergedDnbMap.find(encodePos);
                if (mergedDnbIter != mergedDnbMap.end()){
                    mergedDnbIter->second += count;
                }else{
                    mergedDnbMap[encodePos] = count;
                }
            }
            inDnb.close();
        }
    }
    unordered_map<uint64_t, int>().swap(dnbMap);
    dumpMergedBarcodeList(mergedDnbMapFile);
}

void BarcodeListMerge::addBarcodeList(unordered_map<uint64_t, int>& dnbMap){
    unordered_map<uint64_t, int>::iterator mergedDnbIter;
    for (auto dnbIter = dnbMap.begin(); dnbIter != dnbMap.end(); dnbIter++){
        mergedDnbIter = mergedDnbMap.find(dnbIter->first);
        if (mergedDnbIter != mergedDnbMap.end()){
            mergedDnbIter->second += dnbIter->second;
        }else{
            mergedDnbMap[dnbIter->first] = dnbIter->second;
        }
    }
    dnbMap.clear();
}

void BarcodeListMerge::dumpMergedBarcodeList(string& outfile){
    ofstream outDnb;
    if (ends_with(outfile, ".bin")){
        outDnb.open(outfile, ios::binary);
        boost::archive::binary_oarchive oa(outDnb);
        oa << mergedDnbMap;
    }else{
        outDnb.open(outfile);
        for (auto dnbIter = mergedDnbMap.begin(); dnbIter != mergedDnbMap.end(); dnbIter++){
            uint32_t x = dnbIter->first >> 32;
			uint32_t y = dnbIter->first & 0x00000000FFFFFFFF;
            outDnb << x << "\t" << y << "\t" << dnbIter->second << endl;
        }
    }
    outDnb.close();
}
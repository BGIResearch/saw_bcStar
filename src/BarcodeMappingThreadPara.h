/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
//
// Created by berry on 2021/3/14.
//

#ifndef SOURCE_BARCODEMAPPINGTHREADPARA_H
#define SOURCE_BARCODEMAPPINGTHREADPARA_H
#include "mgzip/mgzf.h"

class BarcodeMappingThreadPara{
public:
	BarcodeMappingThreadPara(Result* stBcPara,Options opt,BarcodeToPositionMulti* bpm){
	  this->stBcPara=stBcPara;
	  this->opt=opt;
	  this->bpm=bpm;
	  this->mgzParser=nullptr;
	}
    BarcodeMappingThreadPara(Result* stBcPara,Options opt,BarcodeToPositionMulti* bpm,mgzipParser* mgzParser){
        this->stBcPara=stBcPara;
        this->opt=opt;
        this->bpm=bpm;
        this->mgzParser=mgzParser;
    }
    Result* stBcPara;
    Options opt;
    BarcodeToPositionMulti* bpm;
    mgzipParser* mgzParser;
    static uint64_t currentTotalReadsNum;
    static uint64_t currentMappedReadsNum;
    static bool checkBcMapRate;
};


#endif //SOURCE_BARCODEMAPPINGTHREADPARA_H

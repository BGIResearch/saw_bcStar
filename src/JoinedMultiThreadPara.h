/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
//
// Created by berry on 2021/3/12.
//

#ifndef SOURCE_JOINEDMULTITHREADPARA_H
#define SOURCE_JOINEDMULTITHREADPARA_H

#include "ReadAlignChunk.h"
#include "BarcodeMappingThreadPara.h"


class JoinedMultiThreadPara{
public:
    JoinedMultiThreadPara(ReadAlignChunk* starPara,BarcodeMappingThreadPara* bmPara){
        this->starPara=starPara;
        this->bmPara=bmPara;
    }
    ReadAlignChunk* starPara;
    BarcodeMappingThreadPara* bmPara;

};


#endif //SOURCE_JOINEDMULTITHREADPARA_H

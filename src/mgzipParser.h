/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
//
// Created by berry on 2021/12/13.
//

#ifndef STAR__MGZIPPARSER_H_
#define STAR__MGZIPPARSER_H_
#include "mgzip/mgzf.h"
#include "common.h"
#include "read.h"
#include "util.h"
#include <string.h>

#define FTEXT 1<<0
#define FHCRC 1<<1
#define FEXTRA 1<<2
#define FNAME 1<<3
#define FCOMMENT 1<<4

class mgzipParser {
 public:
//  string fileName1;
//  string fileName2;
  mgzipParser(string f1,string f2,int threadIdx,int totalThreadsNum);
  ~mgzipParser();
//  ReadPair* getOneReadSE();
  int getReadsFromBlock(ReadPair** readsBuf,int readsNumber);
  int inflateFromBlock(MgzBlock* fp, uint8_t* compressedData,int block_length);
  int moveToNextBlock();
  int readHeaderOfBlock(MGZF* fp);
  MgzBlock* readBlock(MGZF* fp);
  MgzBlock* r1Block;
  MgzBlock* r2Block;
  MGZF* r1;
  MGZF* r2;
 private:
  int threadIdx;
 
  
  uint64_t start1,start2;
  uint64_t end1,end2;
  uint64_t startBlockIdx,endBlockIdx;
  uint64_t curBlockIdx;
  int bufIdx1,bufIdx2;
};

#endif //STAR__MGZIPPARSER_H_

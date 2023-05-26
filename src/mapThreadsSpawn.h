/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef CODE_mapThreadsSpawn
#define CODE_mapThreadsSpawn
#include "Parameters.h"
#include "ReadAlignChunk.h"
#include "options.h"
#include "mgzipParser.h"
#include <set>
void mapThreadsSpawn (Parameters &P, ReadAlignChunk** RAchunk);
Options getStBCpara(Parameters &P);
#endif
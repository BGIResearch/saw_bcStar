#ifndef CODE_BAMbinSortByCoordinate
#define CODE_BAMbinSortByCoordinate
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Genome.h"
#include "zstd/zstd.h"
#include "zstd/common.h"

#include SAMTOOLS_BGZF_H

void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &mapGen);
#if defined(TESTCOMPRESS) || defined(TESTZSTD)
void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &mapGen,uint64_t* arr);
#endif
#endif
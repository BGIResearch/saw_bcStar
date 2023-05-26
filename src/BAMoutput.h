#ifndef CODE_BAMoutput
#define CODE_BAMoutput

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include "Parameters.h"
#include <ext/stdio_filebuf.h>
#include <type_traits>
#include "zstd/zstd.h"
#include "zstd/common.h"

class BAMoutput {//
public:
    //sorted output
    BAMoutput (int iChunk, string tmpDir, Parameters &Pin);
    void coordOneAlign (char *bamIn, uint bamSize, uint iRead);
    void coordBins ();
    void coordFlush ();
    //unsorted output
    BAMoutput (BGZF *bgzfBAMin, Parameters &Pin);
    void unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2);
    void unsortedFlush ();
    void coordUnmappedPrepareBySJout();
    void sortedFlush();

    uint32 nBins; //number of bins to split genome into
    uint* binTotalN; //total number of aligns in each bin
    uint* binTotalBytes;//total size of aligns in each bin
    uint* binTotalUnCompressBytes;//total uncompressed size of aligns in each bin
private:
    uint64 bamArraySize; //this size will be allocated
    char* bamArray; //large array to store the bam alignments, pre-sorted
    uint64 binSize, binSize1;//storage size of each bin
    uint64 binGlen;//bin genomic length
    char **binStart; //pointers to starts of the bins
    uint64 *binBytes, binBytes1;//number of bytes currently written to each bin

    

#ifdef TESTCOMPRESS
    gzFile *gzStream;//output streams for each bin
#else
#ifdef TESTZSTD
	int* zstdFd;
  /* Create the context. */
  ZSTD_CCtx** cctx;
  
#else
  ofstream **binStream;//output streams for each bin
#endif
#endif
    BGZF *bgzfBAM;
    Parameters &P;
    string bamDir;
};

#endif

#include "BAMoutput.h"
#include <sys/stat.h>
#include "GlobalVariables.h"
#include <pthread.h>
#include "serviceFuns.cpp"
#include "ThreadControl.h"
#include "streamFuns.h"

BAMoutput::BAMoutput(int iChunk, string tmpDir, Parameters &Pin) : P(Pin)
{ // allocate bam array

  nBins = P.outBAMcoordNbins;
  binSize = P.chunkOutBAMsizeBytes / nBins;
  bamArraySize = binSize * nBins;
  bamArray = new char[bamArraySize];

  bamDir = tmpDir + to_string((uint)iChunk); // local directory for this thread (ibinStreamChunk)

  mkdir(bamDir.c_str(), P.runDirPerm);
  binStart = new char *[nBins];
  binBytes = new uint64[nBins];

#ifdef TESTCOMPRESS
  gzStream = new gzFile[nBins];
#else
#ifdef TESTZSTD
  cctx = new ZSTD_CCtx *[nBins];
  //    zstdStream=new FILE*[nBins];
  zstdFd = new int[nBins];
#else
  binStream = new ofstream *[nBins];
#endif
#endif
  binTotalN = new uint[nBins];
  binTotalBytes = new uint[nBins];
  binTotalUnCompressBytes = new uint[nBins];
  for (uint ii = 0; ii < nBins; ii++)
  {
    binStart[ii] = bamArray + bamArraySize / nBins * ii;
    binBytes[ii] = 0;

#ifdef TESTCOMPRESS
    gzStream[ii] = gzopen((bamDir + "/" + to_string(ii)).c_str(), "w");
#else
#ifdef TESTZSTD
    //	  zstdStream[ii]=fopen((bamDir +"/"+to_string(ii)).c_str(),"wb");
    zstdFd[ii] = open((bamDir + "/" + to_string(ii)).c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU);
    if (zstdFd[ii] == -1)
    {
      string errMsg = strerror(errno);
      errMsg += "," + bamDir + "/" + to_string(ii);
      sawErrCode(err_fileOpen_failed, errMsg);
      cerr << "Error, " << strerror(errno) << "," << bamDir << "/" << ii << endl;
      exit(errno);
    }
    cctx[ii] = ZSTD_createCCtx();
    if (cctx[ii] == NULL)
    {
      cerr << "Error, ZSTD_createCCtx() failed!" << endl;
      string errMsg = "ZSTD_createCCtx() failed";
      sawErrCode(err_otherAPI_failed, errMsg);
      exit(1);
    }
    //	  CHECK(cctx != NULL, "ZSTD_createCCtx() failed!");

    /* Set any parameters you want.
     * Here we set the compression level, and enable the checksum.
     */
//	  ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, 3);
//	  ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
//	  ZSTD_CCtx_setParameter(cctx, 1, 1);
#else
    binStream[ii] = &ofstrOpen((bamDir + "/" + to_string(ii)).c_str(), ERROR_OUT, P); // open temporary files
#endif
#endif
    binTotalN[ii] = 0;
    binTotalBytes[ii] = 0;
    binTotalUnCompressBytes[ii] = 0;
  };

  binSize1 = binStart[nBins - 1] - binStart[0];
  nBins = 1; // start with one bin to estimate genomic bin sizes
};

BAMoutput::BAMoutput(BGZF *bgzfBAMin, Parameters &Pin) : P(Pin)
{ // allocate BAM array with one bin, streamed directly into bgzf file

  bamArraySize = P.chunkOutBAMsizeBytes;
  bamArray = new char[bamArraySize];
  binBytes1 = 0;
  bgzfBAM = bgzfBAMin;
  // not used
  binSize = 0;

#ifdef TESTCOMPRESS
  gzStream = NULL;
#else
#ifdef TESTZSTD
  cctx = NULL;
  //  	zstdStream=NULL;
  zstdFd = NULL;
#else
  binStream = NULL;
#endif
#endif
  binStart = NULL;
  binBytes = NULL;
  binTotalBytes = NULL;
  binTotalUnCompressBytes = NULL;
  binTotalN = NULL;
  nBins = 0;
};

void BAMoutput::unsortedOneAlign(char *bamIn, uint bamSize, uint bamSize2)
{ // record one alignment to the buffer, write buffer if needed

  if (bamSize == 0)
    return; // no output, could happen if one of the mates is not mapped

  if (binBytes1 + bamSize2 > bamArraySize)
  { // write out this buffer

    if (g_threadChunks.threadBool)
      pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
    if (bgzf_write(bgzfBAM, bamArray, binBytes1) < 0)
    {
      string errMsg = "bgzf_write error";
      sawErrCode(err_fileIO_failed, errMsg);
      cerr << "Error, bgzf_write error" << endl;
      exit(1);
    }
    if (g_threadChunks.threadBool)
      pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);

    binBytes1 = 0; // rewind the buffer
  };

  memcpy(bamArray + binBytes1, bamIn, bamSize);
  binBytes1 += bamSize;
};

void BAMoutput::unsortedFlush()
{ // flush all alignments
  if (g_threadChunks.threadBool)
    pthread_mutex_lock(&g_threadChunks.mutexOutSAM);
  if (bgzf_write(bgzfBAM, bamArray, binBytes1) < 0)
  {
    string errMsg = "bgzf_write error";
    sawErrCode(err_fileIO_failed, errMsg);
    cerr << "Error, bgzf_write error" << endl;
    exit(1);
  }
  if (g_threadChunks.threadBool)
    pthread_mutex_unlock(&g_threadChunks.mutexOutSAM);
  binBytes1 = 0; // rewind the buffer
};

void BAMoutput::coordOneAlign(char *bamIn, uint bamSize, uint iRead)
{

  uint32 *bamIn32;
  uint alignG;
  uint32 iBin = 0;

  if (bamSize == 0)
  {
    return; // no output, could happen if one of the mates is not mapped
  }
  else
  {
    // determine which bin this alignment belongs to
    bamIn32 = (uint32 *)bamIn;
    alignG = (((uint)bamIn32[1]) << 32) | ((uint)bamIn32[2]);
    if (bamIn32[1] == ((uint32)-1))
    { // unmapped
      iBin = P.outBAMcoordNbins - 1;
    }
    else if (nBins > 1)
    { // bin starts have already been determined
      iBin = binarySearch1a<uint64>(alignG, P.outBAMsortingBinStart, (int32)(nBins - 1));
    };
  };

  //     if ( alignG == (uint32) -1 ) {//unmapped alignment, last bin
  //         iBin=nBins-1;
  //     } else {
  //         iBin=(alignG + chrStart)/binGlen;
  //     };

  // write buffer is filled
  if (binBytes[iBin] + bamSize + sizeof(uint) > ((iBin > 0 || nBins > 1) ? binSize : binSize1))
  { // write out this buffer
    if (nBins > 1 || iBin == (P.outBAMcoordNbins - 1))
    { // normal writing, bins have already been determined

#ifdef TESTCOMPRESS
      if (gzwrite(gzStream[iBin], binStart[iBin], binBytes[iBin]) == 0)
      {
        cerr << "Error, gzwrite error" << endl;
        exit(1);
      }
#else
#ifdef TESTZSTD
      //		  ZSTD_inBuffer input = { binStart[iBin], binBytes[iBin], 0 };
      //		  size_t const buffOutSize = binBytes[iBin];
      //		  void *const buffOut = malloc(buffOutSize);
      //            do {
      //			  ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
      //			  size_t const remaining = ZSTD_compressStream2(cctx[iBin], &output, &input, ZSTD_e_end);
      ////			  size_t const wSize = fwrite(buffOut, 1, output.pos, zstdStream[iBin]);
      ////			  if(wSize!=output.pos) {
      ////				cerr << "Error, fwrite error" << endl;
      ////				exit(1);
      ////			  }
      //			  ssize_t wSize=write(zstdFd[iBin],buffOut,output.pos);
      //			  if(wSize<0 || wSize!=output.pos){
      //				cerr<<"Error code: "<<sawErrCode(162)<<endl;
      //				cerr<<"Error, "<<strerror(errno)<<","<<bamDir <<"/"<<iBin<<endl;
      //				exit(errno);
      //			  }
      //			  binTotalBytes[iBin]+=wSize;
      //			}while(input.pos != input.size);
      //		  free(buffOut);

      size_t const cBuffSize = ZSTD_compressBound(binBytes[iBin]);
      char *const cBuff = new char[cBuffSize];
      size_t const cSize = ZSTD_compress(cBuff, cBuffSize, binStart[iBin], binBytes[iBin], 1);
      ssize_t wSize = write(zstdFd[iBin], cBuff, cSize);
      if (wSize < 0 || (size_t)wSize != cSize)
      {
        string errMsg = strerror(errno);
        errMsg += "," + bamDir + "/" + to_string(iBin);
        sawErrCode(err_fileOpen_failed, errMsg);
        cerr << "Error, " << strerror(errno) << "," << bamDir << "/" << iBin << endl;
        exit(errno);
      }
      binTotalBytes[iBin] += wSize;
      binTotalUnCompressBytes[iBin] += binBytes[iBin];
      delete[] cBuff;
#else
      binStream[iBin]->write(binStart[iBin], binBytes[iBin]);
#endif
#endif
      binBytes[iBin] = 0; // rewind the buffer
    }
    else
    { // the first chunk of reads was written in one bin, need to determine bin sizes, and re-distribute reads into bins
      coordBins();
      coordOneAlign(bamIn, bamSize, iRead); // record the current align into the new bins
      return;
    };
  };

  // record this alignment in its bin
  memcpy(binStart[iBin] + binBytes[iBin], bamIn, bamSize);
  binBytes[iBin] += bamSize;
  memcpy(binStart[iBin] + binBytes[iBin], &iRead, sizeof(uint));
  binBytes[iBin] += sizeof(uint);
  //    binTotalBytes[iBin] += bamSize+sizeof(uint);
  binTotalN[iBin] += 1;
  return;
};

void BAMoutput::coordBins()
{                             // define genomic starts for bins
  nBins = P.outBAMcoordNbins; // this is the true number of bins

  // mutex here
  if (P.runThreadN > 1)
    pthread_mutex_lock(&g_threadChunks.mutexBAMsortBins);
  if (P.outBAMsortingBinStart[0] != 0)
  { // it's set to 0 only after the bin sizes are determined
    // extract coordinates and sort
    uint *startPos = new uint[binTotalN[0] + 1]; // array of aligns start positions
    for (uint ib = 0, ia = 0; ia < binTotalN[0]; ia++)
    {
      uint32 *bamIn32 = (uint32 *)(binStart[0] + ib);
      startPos[ia] = (((uint)bamIn32[1]) << 32) | ((uint)bamIn32[2]);
      ib += bamIn32[0] + sizeof(uint32) + sizeof(uint); // note that size of the BAM record does not include the size record itself
    };
    qsort((void *)startPos, binTotalN[0], sizeof(uint), funCompareUint1);

    // determine genomic starts of the bins
    P.inOut->logMain << "BAM sorting: " << binTotalN[0] << " mapped reads\n";
    P.inOut->logMain << "BAM sorting bins genomic start loci:\n";

    P.outBAMsortingBinStart[0] = 0;
    for (uint32 ib = 1; ib < (nBins - 1); ib++)
    {
      P.outBAMsortingBinStart[ib] = startPos[binTotalN[0] / (nBins - 1) * ib];
      P.inOut->logMain << ib << "\t" << (P.outBAMsortingBinStart[ib] >> 32) << "\t" << ((P.outBAMsortingBinStart[ib] << 32) >> 32) << endl;
      // how to deal with equal boundaries???
    };
    delete[] startPos;
  };
  // mutex here
  if (P.runThreadN > 1)
    pthread_mutex_unlock(&g_threadChunks.mutexBAMsortBins);

  // re-allocate binStart
  uint binTotalNold = binTotalN[0];
  char *binStartOld = new char[binSize1];
  memcpy(binStartOld, binStart[0], binBytes[0]);

  binBytes[0] = 0;
  binTotalN[0] = 0;
  binTotalBytes[0] = 0;
  binTotalUnCompressBytes[0] = 0;
  // re-bin all aligns
  for (uint ib = 0, ia = 0; ia < binTotalNold; ia++)
  {
    uint32 *bamIn32 = (uint32 *)(binStartOld + ib);
    uint ib1 = ib + bamIn32[0] + sizeof(uint32); // note that size of the BAM record does not include the size record itself
    coordOneAlign(binStartOld + ib, (uint)(bamIn32[0] + sizeof(uint32)), *((uint *)(binStartOld + ib1)));
    ib = ib1 + sizeof(uint); // iRead at the end of the BAM record
  };
  delete[] binStartOld;
  return;
};

void BAMoutput::coordFlush()
{ // flush all alignments
  if (nBins == 1)
  {
    coordBins();
  };
  //  cout<<"iBin:\t";
  for (uint32 iBin = 0; iBin < P.outBAMcoordNbins; iBin++)
  {

#ifdef TESTCOMPRESS
    gzwrite(gzStream[iBin], binStart[iBin], binBytthread indxes[iBin]);
//	  gzflush(gzStream[iBin]);
#else
#ifdef TESTZSTD
    //	  ZSTD_inBuffer input = { binStart[iBin], binBytes[iBin], 0 };
    //	  size_t const buffOutSize = binBytes[iBin];
    //	  void *const buffOut = malloc(buffOutSize);
    ////	  cout<<iBin<<"\t";
    //	  int finished;
    ////	  ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
    //	  do {
    //		ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
    //		size_t const remaining = ZSTD_compressStream2(cctx[iBin], &output, &input, ZSTD_e_end);
    ////			  size_t const wSize = fwrite(buffOut, 1, output.pos, zstdStream[iBin]);
    ////			  if(wSize!=output.pos) {
    ////				cerr << "Error, fwrite error" << endl;
    ////				exit(1);
    ////			  }
    //		ssize_t wSize=write(zstdFd[iBin],buffOut,output.pos);
    //		if(wSize<0 || wSize!=output.pos){
    //		  cerr<<"Error code: "<<sawErrCode(162)<<endl;
    //		  cerr<<"Error, "<<strerror(errno)<<","<<bamDir <<"/"<<iBin<<endl;
    //		  exit(errno);
    //		}
    //		binTotalBytes[iBin]+=wSize;
    //	  }while(input.pos != input.size);
    ////	  binTotalBytes[iBin]+=wSize;
    ////	  ZSTD_flushStream(cctx[iBin],&output);
    //	  free(buffOut);
    //	  binTotalUnCompressBytes[iBin]+=binBytes[iBin];

    size_t const cBuffSize = ZSTD_compressBound(binBytes[iBin]);
    char *const cBuff = new char[cBuffSize];
    size_t const cSize = ZSTD_compress(cBuff, cBuffSize, binStart[iBin], binBytes[iBin], 1);
    ssize_t wSize = write(zstdFd[iBin], cBuff, cSize);
    if (wSize < 0 || (size_t)wSize != cSize)
    {
      string errMsg = strerror(errno);
      errMsg += "," + bamDir + "/" + to_string(iBin);
      sawErrCode(err_fileOpen_failed, errMsg);
      cerr << "Error, " << strerror(errno) << "," << bamDir << "/" << iBin << endl;
      exit(errno);
    }
    binTotalBytes[iBin] += wSize;
    binTotalUnCompressBytes[iBin] += binBytes[iBin];
    delete[] cBuff;
#else
    binStream[iBin]->write(binStart[iBin], binBytes[iBin]);
    binStream[iBin]->flush();
#endif
#endif
    binBytes[iBin] = 0; // rewind the buffer
  };
  //    cout<<endl;
};

void BAMoutput::coordUnmappedPrepareBySJout()
{ // flush all alignments
}
void BAMoutput::sortedFlush()
{
  using FilebufType = __gnu_cxx::stdio_filebuf<std::ifstream::char_type>;
  static_assert(std::is_base_of<ifstream::__filebuf_type, FilebufType>::value &&
                    (sizeof(FilebufType) == sizeof(ifstream::__filebuf_type)),
                "The filebuf type appears to have extra data members, the cast might be unsafe");
  for (uint32 iBin = 0; iBin < P.outBAMcoordNbins; iBin++)
  {
#ifdef TESTCOMPRESS
    gzclose(gzStream[iBin]);
#else
#ifdef TESTZSTD
    /*
      cout<<iBin<<"\t";
      char* in=new char[1];
      in[0]='a';
      char* buffOut=new char[1024];
      ZSTD_inBuffer input = { in, 0, 0 };
      ZSTD_outBuffer output = {buffOut, 1024, 0};
      size_t const remaining = ZSTD_compressStream2(cctx[iBin], &output, &input, ZSTD_e_end);
      size_t const wSize = fwrite(buffOut, 1, output.pos, zstdStream[iBin]);
      if(wSize!=output.pos) {
        cerr << "Error, fwrite error" << endl;
        exit(1);
      }
      free(buffOut);
      delete[] in;
     */
    ZSTD_freeCCtx(cctx[iBin]);
    if (fsync(zstdFd[iBin]) < 0)
    {
      string errMsg = strerror(errno);
      errMsg += "," + bamDir + "/" + to_string(iBin);
      sawErrCode(err_fileOpen_failed, errMsg);
      cerr << "Error, " << strerror(errno) << "," << bamDir << "/" << iBin << endl;
      exit(errno);
    }
    if (close(zstdFd[iBin]) < 0)
    {
      string errMsg = strerror(errno);
      errMsg += "," + bamDir + "/" + to_string(iBin);
      sawErrCode(err_fileOpen_failed, errMsg);
      cerr << "Error, " << strerror(errno) << "," << bamDir << "/" << iBin << endl;
      exit(errno);
    }
    // check file size on disk
    // binTotalBytes[iBin] record the expected file size
    struct stat fileStat;
    string tmpBinFile = bamDir + "/" + to_string(iBin);
    if (stat(tmpBinFile.c_str(), &fileStat) != 0)
    { // check if file exists
      string errMsg = strerror(errno);
      errMsg += "," + tmpBinFile;
      sawErrCode(err_fileOpen_failed, errMsg);
      cerr << "Error, " << strerror(errno) << "," << tmpBinFile << endl;
      exit(errno);
    };
    if ((uint)fileStat.st_size != (uint)binTotalBytes[iBin])
    {
      string errMsg = "expected size is differ from the file size" + tmpBinFile + "," + to_string(binTotalBytes[iBin]) + "," + to_string(fileStat.st_size);
      sawErrCode(err_sw_exception, errMsg);
      cerr << "Error, expected size is differ from the file size," << tmpBinFile << "," << binTotalBytes[iBin] << "," << fileStat.st_size << endl;
      exit(1);
    }
#else
    binStream[iBin]->rdbuf()->pubsync();
    const int fd = static_cast<FilebufType *>(binStream[iBin]->rdbuf())->fd();
    assert(fd >= 0);
    if (fsync(fd) != 0)
    {
      cerr << "Error, fsync error" << endl;
      exit(1);
    }
    binStream[iBin]->close();

#endif
#endif
  };
};

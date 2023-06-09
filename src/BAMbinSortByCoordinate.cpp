#include "BAMbinSortByCoordinate.h"
#include "ErrorWarning.h"
#include "serviceFuns.cpp"
#include "BAMfunctions.h"

void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &mapGen)
{

  if (binS == 0)
    return; // nothing to do for empty bins
  // allocate arrays
  char *bamIn = new char[binS + 1];
  uint *startPos = new uint[binN * 3];

  uint bamInBytes = 0;
  // load all aligns
  for (uint it = 0; it < nThreads; it++)
  {
    string bamInFile = dirBAMsort + to_string(it) + "/" + to_string((uint)iBin);
    ifstream bamInStream;
    bamInStream.open(bamInFile.c_str(), std::ios::binary | std::ios::ate); // open at the end to get file size
    int64 s1 = bamInStream.tellg();
    if (s1 > 0)
    {
      bamInStream.seekg(std::ios::beg);
      bamInStream.read(bamIn + bamInBytes, s1); // read the whole file
    }
    else if (s1 < 0)
    {
      ostringstream errOut;
      string errMsg="EXITING because of FATAL ERROR: failed reading from temporary file: " + dirBAMsort + to_string(it) + "/" + to_string((uint)iBin);
      sawErrCode(err_fileIO_failed,errMsg);
      errOut << "EXITING because of FATAL ERROR: failed reading from temporary file: " << dirBAMsort + to_string(it) + "/" + to_string((uint)iBin);
      exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
    };
    bamInBytes += bamInStream.gcount();
    bamInStream.close();
    remove(bamInFile.c_str());
  };
  if (bamInBytes != binS)
  {
    //      cout<<"bin number "<<iBin<<endl;
    ostringstream errOut;
    
    string errMsg = "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: ";
    errMsg += "Expected bin size=" + to_string(binS) + " ; size on disk=" + to_string(bamInBytes) + " ; bin number=" + to_string(iBin) + "\n";
    sawErrCode(err_sw_exception,errMsg);
    errOut << "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: ";
    errOut << "Expected bin size=" << binS << " ; size on disk=" << bamInBytes << " ; bin number=" << iBin << "\n";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
  };

  // extract coordinates

  for (uint ib = 0, ia = 0; ia < binN; ia++)
  {
    uint32 *bamIn32 = (uint32 *)(bamIn + ib);
    startPos[ia * 3] = (((uint)bamIn32[1]) << 32) | ((uint)bamIn32[2]);
    startPos[ia * 3 + 2] = ib;
    ib += bamIn32[0] + sizeof(uint32);              // note that size of the BAM record does not include the size record itself
    startPos[ia * 3 + 1] = *((uint *)(bamIn + ib)); // read order
    ib += sizeof(uint);
  };

  // sort
  qsort((void *)startPos, binN, sizeof(uint) * 3, funCompareArrays<uint, 3>);

  BGZF *bgzfBin;
  bgzfBin = bgzf_open((dirBAMsort + "/b" + to_string((uint)iBin)).c_str(), ("w" + to_string((long long)P.outBAMcompression)).c_str());
  if (bgzfBin == NULL)
  {
    ostringstream errOut;
    string errMsg="EXITING because of fatal ERROR: could not open temporary bam file: " + dirBAMsort + "/b" + to_string((uint)iBin) + "\n";
    errMsg+="SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running BCSTAR";
    sawErrCode(err_fileOpen_failed,errMsg);
    errOut << "EXITING because of fatal ERROR: could not open temporary bam file: " << dirBAMsort + "/b" + to_string((uint)iBin) << "\n";
    errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running BCSTAR";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
  };

  outBAMwriteHeader(bgzfBin, P.samHeaderSortedCoord, mapGen.chrNameAll, mapGen.chrLengthAll);
  // send ordered aligns to bgzf one-by-one
  for (uint ia = 0; ia < binN; ia++)
  {
    char *ib = bamIn + startPos[ia * 3 + 2];
    if (bgzf_write(bgzfBin, ib, *((uint32 *)ib) + sizeof(uint32)) < 0)
    {
      string errMsg="bgzf_write error";
      sawErrCode(err_fileIO_failed,errMsg);
      cerr << "Error, bgzf_write error" << endl;
      exit(1);
    }
  };

  if (bgzf_flush(bgzfBin) < 0)
  {
    string errMsg="bgzf_flush error";
    sawErrCode(err_fileIO_failed,errMsg);
    cerr << "Error, bgzf_flush error" << endl;
    exit(1);
  }
  bgzf_close(bgzfBin);
  // release memory
  delete[] bamIn;
  delete[] startPos;
};
#if defined(TESTCOMPRESS) || defined(TESTZSTD)
void BAMbinSortByCoordinate(uint32 iBin, uint binN, uint binS, uint nThreads, string dirBAMsort, Parameters &P, Genome &mapGen, uint64_t *arr)
{

  if (binS == 0)
    return; // nothing to do for empty bins
  // allocate arrays
  char *bamIn = new char[binS + 1];
  uint *startPos = new uint[binN * 3];

  uint bamInBytes = 0;
  // load all aligns
  for (uint it = 0; it < nThreads; it++)
  {
    string bamInFile = dirBAMsort + to_string(it) + "/" + to_string((uint)iBin);
#ifdef TESTCOMPRESS
    gzFile gzIn = gzopen(bamInFile.c_str(), "r");
    uint64_t uncompressedSize = arr[it];
    if (gzread(gzIn, bamIn + bamInBytes, uncompressedSize) < uncompressedSize)
    {
      cerr << "Error, gzread error" << endl;
      exit(1);
    }
    bamInBytes += uncompressedSize;
    gzclose(gzIn);
#else
#ifdef TESTZSTD
    size_t cSize;
    ifstream bamInStream;
    bamInStream.open(bamInFile.c_str(), std::ios::binary | std::ios::ate); // open at the end to get file size
    cSize = bamInStream.tellg();
    char *cBuff = new char[cSize];
    bamInStream.seekg(std::ios::beg);
    bamInStream.read(cBuff, cSize);
    bamInStream.close();

    size_t const dSize = ZSTD_decompress(bamIn + bamInBytes, arr[it], cBuff, cSize);
    //	CHECK_ZSTD(dSize);
    /* When zstd knows the content size, it will error if it doesn't match. */
    //	CHECK(dSize == rSize, "Impossible because zstd will check this condition!");
    if (dSize != arr[it])
    {
      string errMsg="uncompressed size are different between file and recorded in memory!";
      sawErrCode(err_sw_exception,errMsg);
      cerr << bamInFile << "\t" << dSize << "\t" << arr[it] << "\t" << cSize << endl;
      cerr << cBuff[0] << cBuff[1] << cBuff[2] << cBuff[cSize - 3] << cBuff[cSize - 2] << cBuff[cSize - 1] << endl;
      cerr << "Error, uncompressed size are different between file and recorded in memory!" << endl;
      exit(1);
    }
    //	free(rBuff);
    bamInBytes += arr[it];
    delete[] cBuff;
#else
    ifstream bamInStream;
    bamInStream.open(bamInFile.c_str(), std::ios::binary | std::ios::ate); // open at the end to get file size
    int64 s1 = bamInStream.tellg();
    if (s1 > 0)
    {
      bamInStream.seekg(std::ios::beg);
      bamInStream.read(bamIn + bamInBytes, s1); // read the whole file
    }
    else if (s1 < 0)
    {
      ostringstream errOut;
      errOut << "EXITING because of FATAL ERROR: failed reading from temporary file: " << dirBAMsort + to_string(it) + "/" + to_string((uint)iBin);
      sawErrCode(err_fileIO_failed,errOut.str());
      exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
    };
    bamInBytes += bamInStream.gcount();
    bamInStream.close();
#endif
#endif
    remove(bamInFile.c_str());
  };
  if (bamInBytes != binS)
  {
    //      cout<<"bin number "<<iBin<<endl;
    ostringstream errOut;
    string errMsg="EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk.";
    sawErrCode(err_sw_exception,errMsg);
    errOut << "EXITING because of FATAL ERROR: number of bytes expected from the BAM bin does not agree with the actual size on disk: ";
    errOut << "Expected bin size=" << binS << " ; size on disk=" << bamInBytes << " ; bin number=" << iBin << "\n";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
  };

  // extract coordinates

  for (uint ib = 0, ia = 0; ia < binN; ia++)
  {
    uint32 *bamIn32 = (uint32 *)(bamIn + ib);
    startPos[ia * 3] = (((uint)bamIn32[1]) << 32) | ((uint)bamIn32[2]);
    startPos[ia * 3 + 2] = ib;
    ib += bamIn32[0] + sizeof(uint32);              // note that size of the BAM record does not include the size record itself
    startPos[ia * 3 + 1] = *((uint *)(bamIn + ib)); // read order
    ib += sizeof(uint);
  };

  // sort
  qsort((void *)startPos, binN, sizeof(uint) * 3, funCompareArrays<uint, 3>);

  BGZF *bgzfBin;
  bgzfBin = bgzf_open((dirBAMsort + "/b" + to_string((uint)iBin)).c_str(), ("w" + to_string((long long)P.outBAMcompression)).c_str());
  if (bgzfBin == NULL)
  {
    ostringstream errOut;
    string errMsg="EXITING because of fatal ERROR: could not open temporary bam file: " + dirBAMsort + "/b" + to_string((uint)iBin) + "\n";
    errMsg+="SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running BCSTAR";
    sawErrCode(err_sw_exception,errMsg);
    errOut << "EXITING because of fatal ERROR: could not open temporary bam file: " << dirBAMsort + "/b" + to_string((uint)iBin) << "\n";
    errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running BCSTAR";
    exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
  };

  outBAMwriteHeader(bgzfBin, P.samHeaderSortedCoord, mapGen.chrNameAll, mapGen.chrLengthAll);
  // send ordered aligns to bgzf one-by-one
  for (uint ia = 0; ia < binN; ia++)
  {
    char *ib = bamIn + startPos[ia * 3 + 2];
    if (bgzf_write(bgzfBin, ib, *((uint32 *)ib) + sizeof(uint32)) < 0)
    {
      string errMsg="bgzf_write error";
      sawErrCode(err_fileIO_failed,errMsg);
      cerr << "Error, bgzf_write error" << endl;
      exit(1);
    }
  };

  if (bgzf_flush(bgzfBin) < 0)
  {
    string errMsg="bgzf_flush error";
    sawErrCode(err_fileIO_failed,errMsg);
    cerr << "Error, bgzf_flush error" << endl;
    exit(1);
  }
  bgzf_close(bgzfBin);
  // release memory
  delete[] bamIn;
  delete[] startPos;
};
#endif
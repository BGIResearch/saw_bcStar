/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include "ReadAlignChunk.h"
#include "GlobalVariables.h"
#include "ThreadControl.h"
#include "ErrorWarning.h"
#include "SequenceFuns.h"

void ReadAlignChunk::processChunks(BarcodeMappingThreadPara *bmPara)
{ // read-map-write chunks
  // do stBarcodeMapping
  // todo
#ifdef TIME_STAT
  bcTime = 0;
  mapTime = 0;
#endif
  BarcodeToPositionMulti *barcodeToPosMulti = bmPara->bpm;
  //
  noReadsLeft = false;  // true if there no more reads left in the file
  bool newFile = false; // new file marker in the input stream
  //    P.chunkInSizeBytes=1000;
  int readsNumInBlock = P.chunkInSizeBytes / 400;
  //    ofstream tmpOut("bm.out"+to_string(iThread));
  //    ofstream tmpOut2("debug.out"+to_string(iThread));
  int threadReadsInNum = 0;
  string polyAstr = "";
  for (int i = 0; i < bmPara->opt.transBarcodeToPos.umiLen; i++)
  {
    polyAstr += "A";
  }
  int chunkNum = 0;
  //	string mRNApolyA="";
  //	for(int i=0;i<bmPara->opt.polyAnum;i++){
  //	  mRNApolyA+="A";
  //	}
  while (!noReadsLeft)
  { // continue until the input EOF
    //////////////read a chunk from input files and store in memory
    string outstr = "";
#ifdef TIME_STAT
    auto tStart = chrono::high_resolution_clock::now();
#endif
    if (P.outFilterBySJoutStage < 2)
    { // read chunks from input file
      if (!bmPara->opt.isMgzInput)
      {
        if (P.runThreadN > 1)
          pthread_mutex_lock(&g_threadChunks.mutexInRead);
      }

      uint chunkInSizeBytesTotal[2] = {0, 0};
      uint sizeInRead1 = 0;

      string unmappedOut;
      bool hasPosition;
      int num = 0;
      int readNum = 0;
      ReadPair **readsBuf = new ReadPair *[readsNumInBlock];
      if (!bmPara->opt.isMgzInput)
      {
        //            while (num < readsNumInBlock&& !P.reader->mLeft->eof() && !P.reader->mRight->eof()){
        //            while (sizeInRead1 < (float)P.chunkInSizeBytes*0.8){
        while (readNum < readsNumInBlock)
        {
          // read a readpair
          ReadPair *readPair = nullptr;
          if (P.pe)
          {
            readPair = P.reader->read();
          }
          else
          {
            readPair = P.reader->readR1();
          }
          if (readPair == nullptr)
          {
            // check read done all files or not
            if (P.gzfileIdx < P.inputFiles.size())
            {
              string firstFile = P.inputFiles[P.gzfileIdx];
              //					P.inputFiles.pop_front();
              //					if(mgzf_is_mgzf(firstFile.c_str())){
              //					  P.bcOpt.isMgzInput=true;
              //					}else{
              //					  P.bcOpt.isMgzInput=false;
              //					}
              if (!P.bcOpt.isMgzInput)
              {
                if (P.pe)
                {
                  string secondFile = P.inputFiles2[P.gzfileIdx];
                  //						P.inputFiles2.pop_front();
                  P.reader = new FastqReaderPair(firstFile, secondFile, true);
                }
                else
                {
                  P.reader = new FastqReaderPair(firstFile, true);
                  //                            noReadsLeft=false;
                }
              }
              //					P.inOut->logMain  <<timeMonthDayTime()<< " ....."<<P.gzfileIdx<<"\t"<<firstFile<<"\t"<<P.inputFiles2[P.gzfileIdx]<<endl;
              P.inOut->logMain << timeMonthDayTime() << " ....." << P.gzfileIdx << "\t" << firstFile << endl;
              P.gzfileIdx++;
            }
            else
            {
              noReadsLeft = true;
            }
            break;
          }
          readsBuf[readNum] = readPair;
          readNum++;
          threadReadsInNum++;
          bmPara->stBcPara->mTotalRead++;
        }
      }
      else
      {
        Read *l = nullptr;
        Read *r = nullptr;
        for (uint64_t i = 0; i < (uint)readsNumInBlock; i++)
        {
          readsBuf[i] = new ReadPair(l, r);
        }
        int actualReadsNum = bmPara->mgzParser->getReadsFromBlock(readsBuf, readsNumInBlock);
        if (actualReadsNum == 0)
        {
          if (bmPara->mgzParser->r1 != nullptr)
          {
            block_destroy(bmPara->mgzParser->r1Block);
            if (bmPara->mgzParser->r2 != nullptr)
            {
              block_destroy(bmPara->mgzParser->r2Block);
            }
          }
          delete bmPara->mgzParser;
          // check read done all files or not
          if (fileIdx < P.inputFiles.size())
          {
            string firstFile = P.inputFiles[fileIdx];

            string secondFile;
            //				  P.inputFiles.pop_front();
            if (P.pe)
            {
              secondFile = P.inputFiles2[fileIdx];
              //					P.inputFiles2.pop_front();
            }
            fileIdx++;
            //				  delete bmPara->mgzParser;
            if (mgzf_is_mgzf(firstFile.c_str()))
            {
              //					P.bcOpt.isMgzInput=true;
            }
            else
            {
              cerr << "Warning: input file is not in mgz format" << endl;
              P.bcOpt.isMgzInput = false;
            }
            if (P.bcOpt.isMgzInput)
            {
              bmPara->mgzParser = new mgzipParser(firstFile, secondFile, iThread, P.runThreadN);
            }
          }
          else
          {
            noReadsLeft = true;
            //				  delete bmPara->mgzParser;
          }
        }
        else
        {
          bmPara->stBcPara->mTotalRead += actualReadsNum;

          //				for(int i=0;i<readsNumInBlock;i++){
          //				  delete readsBuf[i];
          //				}
          //				delete[] readsBuf;
          //				continue;
        }
        readNum = actualReadsNum;
      }
      //            cout<<iThread<<"\t"<<bmPara->stBcPara->mTotalRead<<endl;
      // put the IO lock at front
      if (!bmPara->opt.isMgzInput)
      {
        if (P.runThreadN > 1)
        {
          pthread_mutex_unlock(&g_threadChunks.mutexInRead);
        }
      }
      // process reads parallel
      int bclen = bmPara->opt.barcodeLen % 2 == 0 ? bmPara->opt.barcodeLen / 2 : bmPara->opt.barcodeLen / 2 + 1;
      int umilen = bmPara->opt.transBarcodeToPos.umiLen % 2 == 0 ? bmPara->opt.transBarcodeToPos.umiLen / 2 : bmPara->opt.transBarcodeToPos.umiLen / 2 + 1;
      int qualSys = 33;
      uint64_t bcMappedNum = 0;
      for (int i = 0; i < readNum; i++)
      {
        ReadPair *readPair = readsBuf[i];
        if (readPair == nullptr)
        {
          string errMsg="code error";
          sawErrCode(err_sw_exception,errMsg);
          cerr << "Error, code error" << endl;
          exit(1);
        }
        Read *or1 = readPair->mLeft;
        Read *or2 = readPair->mRight;
        string umi, umiQual, bcStr;
        if (P.pe)
        {
          // trans PE to writeFq
          // filter low quality umi, add barcode and umi to the read2 name
          string writeFqEncodedBarcode, writeFqEncodedUmi;
          //                    getBarcodeAndUmi(or1,writeFqEncodedBarcode,writeFqEncodedUmi);

          string bcQual = or1->mQuality.substr(bmPara->opt.barcodeStart, bmPara->opt.barcodeLen);
          bcStr = or1->mSeq.mStr.substr(bmPara->opt.barcodeStart, bmPara->opt.barcodeLen);
          if (bmPara->opt.transBarcodeToPos.umiLen == 0)
          {
            umi = "TTTTTTTTTT";
            //					or1->mName+=umi;
          }
          else
          {
            umi = or1->mSeq.mStr.substr(bmPara->opt.transBarcodeToPos.umiStart, bmPara->opt.transBarcodeToPos.umiLen);
          }
          umiQual = or1->mQuality.substr(bmPara->opt.transBarcodeToPos.umiStart, bmPara->opt.transBarcodeToPos.umiLen);
          bmPara->stBcPara->mBarcodeProcessor->barcodeStatAndFilter(bcQual);

          for (uint64_t i = 0; i < umiQual.size(); i++)
          {
            if (umiQual[i] - qualSys >= 30)
            {
              bmPara->stBcPara->mBarcodeProcessor->umiQ30++;
              bmPara->stBcPara->mBarcodeProcessor->umiQ20++;
              bmPara->stBcPara->mBarcodeProcessor->umiQ10++;
            }
            else if (umiQual[i] - qualSys >= 20)
            {
              bmPara->stBcPara->mBarcodeProcessor->umiQ20++;
              bmPara->stBcPara->mBarcodeProcessor->umiQ10++;
            }
            else if (umiQual[i] - qualSys >= 10)
            {
              bmPara->stBcPara->mBarcodeProcessor->umiQ10++;
            }
          }
          uint64_t bcInt = 0;
          if (bmPara->opt.encodeRule == "ACTG")
          {
            bcInt = seqEncode(or1->mSeq.mStr.c_str(), bmPara->opt.barcodeStart, bmPara->opt.barcodeLen);
          }
          else if (bmPara->opt.encodeRule == "ACGT")
          {
            bcInt = seqEncode2(or1->mSeq.mStr.c_str(), bmPara->opt.barcodeStart, bmPara->opt.barcodeLen);
          }
          else
          {
            string errMsg="encodeRule only can be set to ACGT or ACTG";
            sawErrCode(err_param_invalid_bcmap,errMsg);
            cerr << "Error, encodeRule only can be set to ACGT or ACTG" << endl;
            exit(1);
          }
          //                    cout<<or1->mSeq.mStr<<"\t"<<bmPara->opt.transBarcodeToPos.umiStart<<"\t"<<bmPara->opt.transBarcodeToPos.umiLen<<endl;
          uint64_t umiInt = 0;
          if (bmPara->opt.transBarcodeToPos.umiLen == 0)
          {
            umiInt = seqEncode2(umi.c_str(), 0, 10);
          }
          else
          {
            umiInt = seqEncode2(or1->mSeq.mStr.c_str(), bmPara->opt.transBarcodeToPos.umiStart, bmPara->opt.transBarcodeToPos.umiLen);
          }
          //                    cout<<umiInt<<endl;
          char *bcBuf = new char[bclen];
          char *umiBuf = new char[umilen];
          sprintf(bcBuf, "%llX", (long long unsigned int)bcInt);
          sprintf(umiBuf, "%llX", (long long unsigned int)umiInt);
          //                    cout<<umiBuf<<endl;
          // remove "/2" in read2 name
          if (or2->mName.rfind("/2") != string::npos)
          {
            or2->mName.erase(or2->mName.rfind("/2"));
          }
          else if (or2->mName.rfind("/1") != string::npos)
          {
            or2->mName.erase(or2->mName.rfind("/1"));
          }
          or2->mName += " ";
          or2->mName += bcBuf;
          or2->mName += " ";
          or2->mName += umiBuf;
          //                    cout<<or2->mName<<endl;
          delete[] bcBuf;
          delete[] umiBuf;
          or1 = or2;
        }
        if (umi.empty())
        {
          bmPara->stBcPara->mBarcodeProcessor->getUmi(or1->mName, umi);
        }
        //				bc=or1->mSeq.mStr.substr(bmPara->opt.barcodeStart,bmPara->opt.barcodeLen);

        // stat Quality of sequence
        if (bmPara->stBcPara->readLen == 0)
        {
          bmPara->stBcPara->readLen = or1->mQuality.size();
        }
        else
        {
          //#ifndef MRNAPOLYA
          if (bmPara->stBcPara->readLen != or1->mQuality.size())
          {
            cerr << "read length are not the same in the fastq file, " << or1->toString() << endl;
            exit(1);
          }
          //#endif
        }
        for (uint64_t i = 0; i < or1->mQuality.size(); i++)
        {
          if ((or1->mQuality)[i] - qualSys >= 30)
          {
            bmPara->stBcPara->seqQ30++;
            bmPara->stBcPara->seqQ20++;
            bmPara->stBcPara->seqQ10++;
          }
          else if ((or1->mQuality)[i] - qualSys >= 20)
          {
            bmPara->stBcPara->seqQ20++;
            bmPara->stBcPara->seqQ10++;
          }
          else if ((or1->mQuality)[i] - qualSys >= 10)
          {
            bmPara->stBcPara->seqQ10++;
          }
        }
        bmPara->stBcPara->totalBases += or1->mQuality.size();

        // polyA detect
        if (bmPara->opt.polyAnum > 0)
        {
          int polyApos = findPolyA(or1->mSeq.mStr, bmPara->opt.polyAnum, bmPara->opt.mismatchInPolyA);
          if (polyApos >= 0)
          {
            bmPara->stBcPara->containPolyA++;
            if (polyApos < 30)
            {
              bmPara->stBcPara->filterByPolyA++;
              continue;
            }
            else
            {
              or1->mSeq.mStr = or1->mSeq.mStr.substr(0, polyApos);
              or1->mQuality = or1->mQuality.substr(0, polyApos);
            }
          }
        }
        int Nindex = -1;
        int Nnum = count(bcStr.begin(), bcStr.end(), 'N');
        if (Nnum > 1)
        {
          bmPara->stBcPara->tooManyNinBarcode++;
          continue;
        }
        else if (Nnum == 1)
        {
          Nindex = bcStr.find('N');
        }
        hasPosition = bmPara->stBcPara->mBarcodeProcessor->process(or1, P.outSAMattributes.at(0), Nindex);
        //#ifdef TESTLRMOVE1
        ////			  for(int i=0;i<readNum;i++){
        ////                delete readsBuf[i];
        ////            }
        ////            delete[] readsBuf;
        //                continue;
        //#endif
        //			  hasPosition=chunkNum%3==0?true:false;
        //                hasPosition=bmPara->stBcPara->mBarcodeProcessor->process(or1,or2);
        if (hasPosition)
        {
          bcMappedNum++;
          //                    if(P.pe){

          int qualSys = 33;
          if (umi.find("N") != string::npos)
          {
            bmPara->stBcPara->mBarcodeProcessor->umiNFilterReads++;
            continue;
          }
          else
          {
            if (umi == polyAstr)
            {
              bmPara->stBcPara->mBarcodeProcessor->umiPloyAFilterReads++;
              continue;
            }
            if (P.pe)
            {
              int lowQualNum = 0;
              bool filterFlag = false;
              for (uint64_t i = 0; i < umiQual.size(); i++)
              {
                if (umiQual[i] - qualSys < 11)
                {
                  lowQualNum++;
                }
                if (lowQualNum >= 2)
                {
                  filterFlag = true;
                  break;
                }
              }
              if (filterFlag)
              {
                bmPara->stBcPara->mBarcodeProcessor->umiQ10FilterReads++;
                continue;
              }
            }
          }
          //                    }else{
          //					  if(umi==polyAstr){
          //						bmPara->stBcPara->mBarcodeProcessor->umiPloyAFilterReads++;
          //						continue;
          //					  }
          //                    }

          num++;
          outstr += or1->toString();
          sizeInRead1 += outstr.size();
        }
        else if (barcodeToPosMulti->mUnmappedWriter)
        {
          //                        unmappedOut+=or1->toString();
          continue;
        }
        //                delete readPair;
        //                delete or1;
        //                delete or2;
      }

      // check barcode mapping rate

      if (BarcodeMappingThreadPara::checkBcMapRate)
      {
        pthread_mutex_lock(&g_threadChunks.mutexInRead);
        BarcodeMappingThreadPara::currentTotalReadsNum += readNum;
        BarcodeMappingThreadPara::currentMappedReadsNum += bcMappedNum;
        uint64_t checkPoint = 100 * 1000 * 1000;
        //			  cout<<BarcodeMappingThreadPara::currentTotalReadsNum<<endl;
        if (BarcodeMappingThreadPara::currentTotalReadsNum >= checkPoint)
        {
          double bcMapRate = double(BarcodeMappingThreadPara::currentMappedReadsNum) / double(BarcodeMappingThreadPara::currentTotalReadsNum);
          //				cout<<bcMapRate<<"\t"<<BarcodeMappingThreadPara::currentTotalReadsNum<<"\t"<<BarcodeMappingThreadPara::currentMappedReadsNum<<endl;
          if (bcMapRate < bmPara->bpm->mOptions->bcMappingRateCutoff)
          {
            cerr << fixed << setprecision(2);
            string errMsg="barcode mapping rate is too low:"+to_string(bcMapRate);
            sawErrCode(err_cidRate_low,errMsg);
            cerr << "Error, barcode mapping rate is too low," << bcMapRate << endl;
            BarcodeMappingThreadPara::checkBcMapRate = false;
            exit(1);
          }
          else
          {
            BarcodeMappingThreadPara::checkBcMapRate = false;
          }
        }
        pthread_mutex_unlock(&g_threadChunks.mutexInRead);
      }
      //#ifndef TESTLRMOVE1
      if (bmPara->opt.isMgzInput)
      {
        for (int i = 0; i < readsNumInBlock; i++)
        {
          delete readsBuf[i];
        }
      }
      else
      {
        for (int i = 0; i < readNum; i++)
        {
          delete readsBuf[i];
        }
      }
      //#endif
      delete[] readsBuf;
      chunkNum++;
      //		  cout<<get_local_time()<<"...thread "<<iThread<<" : processed chunk number "<<chunkNum<<endl;
      //			continue;
      //            cout<<get_local_time()<<"...thread "<<iThread<<". Number of reads have been processed by ST_barcodeMapping: "<<threadReadsInNum<<"\t"<<P.iReadAll<<endl;
      //            tmpOut<<outstr<<endl;
      //            cout<<outstr<<endl;
      istringstream read2stream(outstr);
      bool startStatus = read2stream.peek();
      P.readNmatesIn = 1; // SE data
      P.readNmates = 1;
      while (startStatus && !read2stream.eof())
      {
        P.iReadAll++;
        if (P.outFilterBySJoutStage != 2)
        { // not the 2nd stage of the 2-stage mapping, read ID from the 1st read

          string readID;
          read2stream >> readID;
          if (P.outSAMreadIDnumber)
          {
            readID = "@" + to_string(P.iReadAll);
          };
          // read the second field of the read name line
          char passFilterIllumina = 'N';
          if (read2stream.peek() != '\n')
          { // 2nd field exists
            string field2;
            read2stream >> field2;
            if (field2.length() >= 3 && field2.at(1) == ':' && field2.at(2) == 'Y' && field2.at(3) == ':')
              passFilterIllumina = 'Y';
          };
          readID += ' ' + to_string(P.iReadAll) + ' ' + passFilterIllumina + ' ' + to_string(P.readFilesIndex);

          // ignore the rest of the read name for both mates
          for (uint imate = 0; imate < P.readNmatesIn; imate++)
            read2stream.ignore(DEF_readNameSeqLengthMax, '\n');

          if (P.pSolo.type == 1)
          { // record barcode sequence
            string seq1;
            getline(P.inOut->readIn[1], seq1);
            if (seq1.size() != P.pSolo.bL && P.pSolo.bL > 0)
            {
              ostringstream errOut;
              errOut << "EXITING because of FATAL ERROR in input read file: the total length of barcode sequence is " << seq1.size() << " not equal to expected " << P.pSolo.bL << "\n";
              errOut << "Read ID=" << readID << "   Sequence=" << seq1 << "\n";
              errOut << "SOLUTION: make sure that the barcode read is the second in --readFilesIn and check that is has the correct formatting\n";
              errOut << "          If UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength\n";
              sawErrCode(err_sw_exception,errOut.str());
              exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };
            readID += ' ' + seq1;
            P.inOut->readIn[1].ignore(DEF_readNameSeqLengthMax, '\n'); // skip to the end of 3rd ("+") line
            getline(P.inOut->readIn[1], seq1);                         // read qualities
            readID += ' ' + seq1;
          };

          // copy the same readID to both mates
          for (uint imate = 0; imate < P.readNmates; imate++)
          {
            chunkInSizeBytesTotal[imate] += 1 + readID.copy(chunkIn[imate] + chunkInSizeBytesTotal[imate], readID.size(), 0);
            chunkIn[imate][chunkInSizeBytesTotal[imate] - 1] = '\n';
          };
        };
        // copy 3 (4 for stage 2) lines: sequence, dummy, quality
        for (uint imate = 0; imate < P.readNmates; imate++)
        {
          for (uint iline = (P.outFilterBySJoutStage == 2 ? 0 : 1); iline < 4; iline++)
          {
            read2stream.getline(chunkIn[imate] + chunkInSizeBytesTotal[imate], DEF_readNameSeqLengthMax + 1);
            chunkInSizeBytesTotal[imate] += read2stream.gcount();
            chunkIn[imate][chunkInSizeBytesTotal[imate] - 1] = '\n';
          };
        };

        if (newFile)
        {
          P.inOut->readIn[0] >> P.readFilesIndex;
          pthread_mutex_lock(&g_threadChunks.mutexLogMain);
          P.inOut->logMain << "Starting to map file # " << P.readFilesIndex << "\n";
          for (uint imate = 0; imate < P.readFilesNames.size(); imate++)
          {
            P.inOut->logMain << "mate " << imate + 1 << ":   " << P.readFilesNames.at(imate).at(P.readFilesIndex) << "\n";
            P.inOut->readIn[imate].ignore(numeric_limits<streamsize>::max(), '\n');
          };
          P.inOut->logMain << flush;
          pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
          newFile = false;
        };
      };
      read2stream.clear();
      // TODO: check here that both mates are zero or non-zero
      //            if (chunkInSizeBytesTotal[0]==0) {
      //                noReadsLeft=true; //true if there no more reads left in the file
      iChunkIn = g_threadChunks.chunkInN; // to keep things consistent
      g_threadChunks.chunkInN++;
      //                if(!P.inputFiles.empty()){
      //                    string firstFile=P.inputFiles.front();
      //                    P.inputFiles.pop_front();
      //                    if(P.pe){
      //                        string secondFile=P.inputFiles2.front();
      //                        P.inputFiles2.pop_front();
      //                        P.reader=new FastqReaderPair(firstFile,secondFile,true);
      //                    }else{
      //                        P.reader=new FastqReaderPair(firstFile,true);
      //                        noReadsLeft=false;
      //                    }
      //                }
      //            } else {
      //                noReadsLeft=false;
      //                iChunkIn=g_threadChunks.chunkInN;
      //                g_threadChunks.chunkInN++;
      //            };

      for (uint imate = 0; imate < P.readNmates; imate++)
      {
        chunkIn[imate][chunkInSizeBytesTotal[imate]] = '\n'; // extra empty line at the end of the chunks
        //              cout<<P.chunkInSizeBytesArray<<"\t"<<chunkInSizeBytesTotal[imate]<<endl;
      }
    }
    else
    { // read from one file per thread
      noReadsLeft = true;
      for (uint imate = 0; imate < P.readNmates; imate++)
      {
        RA->chunkOutFilterBySJoutFiles[imate].flush();
        RA->chunkOutFilterBySJoutFiles[imate].seekg(0, ios::beg);
        //                RA->readInStream[imate]=& RA->chunkOutFilterBySJoutFiles[imate];
      };
    };
//        RA->readInStream[0]->str(outstr);
//		for(uint i=0;i<P.readNmates;i++) {
//		  RA->readInStream[i]->clear();
////		  RA->readInStream[i]->str(chunkIn[i]);
//		  RA->readInStream[0]->rdbuf()->pubsetbuf(chunkIn[i],P.chunkInSizeBytesArray);
//		}
//        RA->readInStream[0]->rdbuf()->pubsetbuf(chunkIn[0],P.chunkInSizeBytesArray);
// readInStream[ii]->rdbuf()->pubsetbuf(chunkIn[ii],P.chunkInSizeBytesArray);
//       RA->readInStream[ii]=readInStream[ii];
//        tmpOut<<chunkIn[0]<<endl;
//        exit(1);
//        cout<<get_local_time()<<"...thread "<<iThread<<". before mapChunk "<<endl;
//		continue;
#ifdef TIME_STAT
    auto tEnd = chrono::high_resolution_clock::now();
    bcTime += chrono::duration<double, nano>(tEnd - tStart).count();
#endif
#ifdef TIME_STAT
    tStart = chrono::high_resolution_clock::now();
#endif
    mapChunk();
#ifdef TIME_STAT
    tEnd = chrono::high_resolution_clock::now();
    mapTime += chrono::duration<double, nano>(tEnd - tStart).count();
#endif
    //        cout<<get_local_time()<<"...thread "<<iThread<<". after mapChunk "<<endl;
    if (iThread == 0 && P.runThreadN > 1 && P.outSAMorder == "PairedKeepInputOrder")
    { // concatenate Aligned.* files
      chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
    };
    //        cout<<get_local_time()<<"...thread "<<iThread<<". cycle done"<<endl;
    //        cout<<get_local_time()<<"...thread "<<iThread<<". whether noReadsLeft:"<<noReadsLeft<<endl;

  }; // cycle over input chunks

  //    cout<<get_local_time()<<"...thread "<<iThread<<" before P.outFilterBySJoutStage"<<endl;
  if (P.outFilterBySJoutStage != 1 && RA->iRead > 0)
  { // not the first stage of the 2-stage mapping
    if (P.outBAMunsorted)
      chunkOutBAMunsorted->unsortedFlush();
    if (P.outBAMcoord)
      chunkOutBAMcoord->coordFlush();
    if (chunkOutBAMquant != NULL)
      chunkOutBAMquant->unsortedFlush();

    // the thread is finished mapping reads, concatenate the temp files into output files
    if (P.pCh.segmentMin > 0)
    {
      chunkFstreamCat(RA->chunkOutChimSAM, P.inOut->outChimSAM, P.runThreadN > 1, g_threadChunks.mutexOutChimSAM);
      chunkFstreamCat(*RA->chunkOutChimJunction, P.inOut->outChimJunction, P.runThreadN > 1, g_threadChunks.mutexOutChimJunction);
    };
    if (P.outReadsUnmapped == "Fastx")
    {
      if (P.runThreadN > 1)
        pthread_mutex_lock(&g_threadChunks.mutexOutUnmappedFastx);

      for (uint ii = 0; ii < P.readNmatesIn; ii++)
      {
        chunkFstreamCat(RA->chunkOutUnmappedReadsStream[ii], P.inOut->outUnmappedReadsStream[ii], false, g_threadChunks.mutexOutUnmappedFastx);
      };

      if (P.runThreadN > 1)
        pthread_mutex_unlock(&g_threadChunks.mutexOutUnmappedFastx);
    };
  };
  //    cout<<get_local_time()<<"...thread "<<iThread<<" processChunks almost done"<<endl;
  //    tmpOut.close();
  //    tmpOut2.close();
  //    cout<<get_local_time()<<"...thread "<<iThread<<" threadReadsInNum: "<<threadReadsInNum<<endl;
  if (P.runThreadN > 1)
    pthread_mutex_lock(&g_threadChunks.mutexLogMain);
  //    if(iThread==0){
  //        bmPara->stBcPara->mBarcodeProcessor->bpmap->clear();
  //    }
  //  cout <<"query time of bc: "<<bmPara->stBcPara->mBarcodeProcessor->queryTime<<endl;
  // cout<<"thread indx:\t"<<iThread<<endl;
  chunkOutBAMcoord->sortedFlush();
  P.inOut->logMain << "Completed: thread #" << iThread << endl;
  if (P.runThreadN > 1)
    pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
  //    cout<<get_local_time()<<"...thread "<<iThread<<" processChunks done"<<endl;
};

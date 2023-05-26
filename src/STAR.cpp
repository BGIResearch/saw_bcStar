/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "zlib.h"
#include <signal.h>
#include <execinfo.h>
#include "IncludeDefine.h"
#include "Parameters.h"
#include "SequenceFuns.h"
#include "Genome.h"
#include "Chain.h"
#include "ReadAlignChunk.h"
#include "ReadAlign.h"
#include "Stats.h"
#include "genomeGenerate.h"
#include "outputSJ.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "sysRemoveDir.h"
#include "BAMfunctions.h"
#include "Transcriptome.h"
#include "BAMbinSortByCoordinate.h"
#include "BAMbinSortUnmapped.h"
#include "signalFromBAM.h"
#include "mapThreadsSpawn.h"
#include "ErrorWarning.h"
#include "SjdbClass.h"
#include "sjdbInsertJunctions.h"
#include "Variation.h"
#include "Solo.h"
#include "mgzip/mgzf.h"

#include "bam_cat.h"
#include "htslib/sam.h"
#include "parametersDefault.xxd"
#define BTLOG "_bt_log"
#define MAJOR_VERSION "1"
#define MINOR_VERSION "0.11"
typedef void (*sighandler_t)(int);
sighandler_t signal(int signum, sighandler_t handler);
void printBackTrace()
{
  int size = 16;
  void *array[16];
  int stack_num = backtrace(array, size);
  char **stacktrace = backtrace_symbols(array, stack_num);
  ofstream f(BTLOG, ios::app);
  if (!f)
  {
    string errMsg = BTLOG;
    errMsg = "cannot open such file," + errMsg;
    sawErrCode(err_fileOpen_failed, errMsg);
    cerr << "Error, cannot open such file," << BTLOG << endl;
    for (int i = 0; i < stack_num; ++i)
    {
      cerr << stacktrace[i] << endl;
    }
  }
  else
  {
    for (int i = 0; i < stack_num; ++i)
    {
      f << stacktrace[i] << endl;
    }
    f.close();
  }
  free(stacktrace);
}
void processSignal(int sig)
{
  string errMsg = "get signal:\t" + to_string(sig);
  sawErrCode(err_sw_exception, errMsg);
  cerr << "get signal:\t" << sig << endl;
  switch (sig)
  {
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 10:
  case 11:
  case 12:
  {
    printBackTrace();
    break;
  }
  default:
  {
    break;
  }
  }
  if (sig == 1)
  {
    printBackTrace();
  }
  signal(sig, SIG_DFL);
  kill(getpid(), sig);
}
void myExit()
{
  ofstream f(BTLOG);
  if (!f)
  {

    string errMsg = BTLOG;
    errMsg = "cannot open such file," + errMsg;
    sawErrCode(err_fileOpen_failed, errMsg);
    cerr << "Error, cannot open such file," << BTLOG << endl;
    cerr << "program exit" << endl;
    printBackTrace();
  }
  else
  {
    f << "program exit" << endl;
    f.close();
    printBackTrace();
  }
}
void printVersionAndAuthor()
{
  cout << "Program: bcStar\n";
  cout << "Version: " << MAJOR_VERSION << "." << MINOR_VERSION << endl;
  cout << "Contact: GongChun<gongchun@genomics.cn>" << endl;
}
void printBcParaUsage()
{
  cmdline::parser cmd;
  // input/output
  cmd.add<string>("in", 'i', "the first sequencing fastq file path of read or bin file of  first sequencing", "");
  cmd.add<string>("in1", 'I', "the second sequencing fastq file path of read1 or bin file of second sequencing", false, "");
  cmd.add<string>("in2", 0, "the second sequencing fastq file or bin file path of read1", false, "");
  cmd.add<string>("encodeRule", 0, "encode Rule of 4 bases, default:ACGT which means A=>0,C=>1,G=>2,T=>3, and ACTG means T=>2,G=>3", false, "ACGT");
  cmd.add<string>("barcodeReadsCount", 0, "the mapped barcode list file with reads count per barcode.", false, "");
  cmd.add<string>("maskFile", 'M', "mask file path", false, "");
  cmd.add<string>("platform", 0, "sequence platform [SEQ500, T1, T5, T10]", false, "T1");
  cmd.add<string>("out", 'O', "output file prefix or fastq output file of read1", true, "");
  cmd.add<string>("report", 0, "logging file path.", false, "");
  cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
  cmd.add<string>("unmappedOut", 0, "output file path for barcode unmapped reads, if this path isn't given, discard the reads.", false, "");
  cmd.add<string>("q10Tiff", 0, "generate q10 heatmap tiff when get barcode to position map.", false, "");
  cmd.add<string>("fov", 0, "fov range in format minCol-maxCol_minRow-maxRow.", false, "0-17_0-96");
  cmd.add<string>("chipID", 0, "chip id", false, "");
  cmd.add<int>("barcodeLen", 'l', "barcode length, default  is 25", false, 25);
  cmd.add<int>("barcodeStart", 0, "barcode start position", false, 0);
  cmd.add<int>("umiRead", 0, "read1 or read2 contains the umi sequence.", false, 1);
  cmd.add<int>("umiStart", 0, "umi start position. if the start postion is negative number, no umi sequence will be found", false, 40);
  cmd.add<int>("umiLen", 0, "umi length.", false, 10);
  cmd.add<string>("fixedSequence", 0, "fixed sequence in read1 that will be filtered.", false, "");
  cmd.add<int>("fixedStart", 0, "fixed sequence start position can by specied.", false, -1);
  cmd.add<string>("fixedSequenceFile", 0, "file contianing the fixed sequences and the start position, one sequence per line in the format: TGCCTCTCAG\t-1. when position less than 0, means wouldn't specified", false, "");
  cmd.add<int>("turnFovDegree", 0, "degree of every fov will be turned.", false, 0);
  cmd.add<long>("mapSize", 0, "bucket size of the new unordered_map.", false, 0);
  cmd.add<int>("mismatch", 0, "max mismatch is allowed for barcode overlap find.", false, 0);
  cmd.add<int>("segment", 0, "barcode segment number of first sequencing read.", false, 1);
  cmd.add<string>("rc", 0, "true means get the reverse complement barcode to barcode map. all means get both forward and reverse complement barcode to barcode map", false, "false");
  cmd.add<string>("readidSep", 0, "number of characters will be trimed from readid to get read position in slide. If not given this will be automatically get from fastq file.", false, "/");
  cmd.add<int>("action", 0, "chose one action you want to run [barcode_stat = 1, barcode_overlap = 2, get_barcode_position_map = 3, map_barcode_to_slide = 4, merge_barcode_list = 5].", false, 1);
  cmd.add<int>("thread", 'w', "number of thread that will be used to run.", false, 2);
  cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");
  //  cmd.add("useF14",0,"use F14 to replace unordered_map");
  cmd.add<uint64_t>("bcNum", 0, "barcode number in maskFile.", false, 0);
  cmd.add<float>("bcMappingRateCutoff", 0, "lower bound of barcode mapping rate.", false, 0);
  cmd.add<int>("polyAnum", 0, "cutoff of continuous A bases number which would be considered as a polyA.", false, 0);
  cmd.add<int>("mismatchInPolyA", 0, "mismatch in polyA.", false, 0);
  //  cmd.add("isMgzInput",0,"the format of input fastq files is mGzip");
  cout << "Parameters of barcode mapping. (these should be written in a file as the value of parameter --bcPara)" << endl;
  cmd.parse_check("");
}
void usage(int usageType)
{
  printVersionAndAuthor();
  cout << "Usage: bcStar  [options]... --genomeDir REFERENCE   --readFilesIn R1.fq R2.fq --bcPara <barcodeMapping parameter file>\n";
  cout << "Parameters of barcode mapping should be written in a file as the value of parameter --bcPara" << endl;
  cout << "e.g. echo --in=input.gz >bcp.txt, and run bcStar [other options] --bcPara bcp.txt" << endl;
  //    cout <<"Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2019\n\n";
  //    cout <<"For more details see:\n";
  //    cout <<"<https://github.com/alexdobin/BCSTAR>\n";
  //    cout <<"<https://github.com/alexdobin/BCSTAR/blob/master/doc/STARmanual.pdf>\n";
  cout << "To list parameters of barcode mapping, run bcStar --helpBcPara" << endl;
  if (usageType == 0)
  { // brief
    cout << "To list all parameters, run bcStar --help\n";
  }
  else
  { // full
    cout << '\n'
         << parametersDefault;
    printBcParaUsage();
  };
  exit(0);
}
uint64_t BarcodeMappingThreadPara::currentTotalReadsNum = 0;
uint64_t BarcodeMappingThreadPara::currentMappedReadsNum = 0;
bool BarcodeMappingThreadPara::checkBcMapRate = true;
int main(int argInN, char *argIn[])
{
  atexit(myExit);
  for (int i = 1; i < 31; i++)
  {
    if (i != 9 && i != 17 && i != 1 && i != 16 && i != 19 && i != 20 && i != 23 && i != 28 && i != 29)
    {
      signal(i, processSignal);
    }
  }
  //	uint64_t testBcNum=atoi(argIn[1]);
  ////  testBcNum=testBcNum*1024*1024;
  //	folly::F14ValueMap<uint64_t,Position1> testF14;
  //	testF14.reserve(testBcNum);
  //	cout<<testBcNum<<"\t"<<testF14.getAllocatedMemorySize()<<endl;
  //	exit(1);
  // If no argument is given, or the first argument is either '-h' or '--help', run usage()
  if (argInN == 1)
  {
    usage(0);
  }
  else if (argInN == 2)
  {
    string firstArt = argIn[1];
    if (strcmp("-h", argIn[1]) == 0 || strcmp("--help", argIn[1]) == 0)
    {
      usage(1);
    }
    else if (firstArt == "--helpBcPara")
    {
      cout << "bcStar version: " << MAJOR_VERSION << "." << MINOR_VERSION << endl;
      printBcParaUsage();
    }
    else
    {
      usage(1);
    }
  };

  //  	MPI_Init(&argInN, &argIn);
  time(&g_statsAll.timeStart);
  string stBCparaFile = "";
  char **paraStr = new char *[argInN + 1];
  bool st = false;
  int j = 0;
  for (int i = 0; i < argInN; i++)
  {
    string tmpStr = argIn[i];
    if (tmpStr == "--stParaFile")
    {
      st = true;
    }
    else
    {
      if (st)
      {
        stBCparaFile = tmpStr;
        st = false;
      }
      else
      {
        paraStr[j] = argIn[i];
        j++;
      }
    }
  }
  Parameters P; // all parameters
  P.stBCparaFile = stBCparaFile;
  ifstream paraIn(stBCparaFile);
  string tmpLine;
  string r2file;
  while (getline(paraIn, tmpLine))
  {
    if (tmpLine.find("in1") == 0)
    {
      vector<string> eles;
      split(tmpLine, eles, "=");
      r2file = eles[1];
      break;
    }
  }
  //    paraStr[j++]="--readFilesIn";
  //    paraStr[j++]=(char*)r2file.c_str();
  //    paraStr[j++]="--readFilesCommand gzip -dc";
  P.inputParameters(j, paraStr);

  //    P.inputParameters(argInN, argIn);
  *(P.inOut->logStdOut) << get_local_time() << " ..... started BCSTAR run\n"
                        << flush;

  // generate genome
  if (P.runMode == "alignReads")
  { // continue
  }
  else if (P.runMode == "genomeGenerate")
  {
    Genome mainGenome(P);
    mainGenome.genomeGenerate();
    (void)sysRemoveDir(P.outFileTmp);
    P.inOut->logMain << "DONE: Genome generation, EXITING\n"
                     << flush;
    exit(0);
  }
  else if (P.runMode == "liftOver")
  {
    for (uint ii = 0; ii < P.pGe.gChainFiles.size(); ii++)
    {
      Chain chain(P, P.pGe.gChainFiles.at(ii));
      chain.liftOverGTF(P.pGe.sjdbGTFfile, P.outFileNamePrefix + "GTFliftOver_" + to_string(ii + 1) + ".gtf");
      P.inOut->logMain << "DONE: lift-over of GTF file, EXITING\n"
                       << flush;
      exit(0);
    };
  }
  else
  {
    P.inOut->logMain << "EXITING because of INPUT ERROR: unknown value of input parameter runMode=" << P.runMode << endl
                     << flush;
    exit(1);
  };
#ifdef TIME_STAT
  P.loadGenomeTime = P.loadMaskTime = P.sortBamTime = 0;
  auto tStart = chrono::high_resolution_clock::now();
#endif
  Genome mainGenome(P);
  mainGenome.genomeLoad();
#ifdef TIME_STAT
  auto tEnd = chrono::high_resolution_clock::now();
  P.loadGenomeTime = chrono::duration<double, nano>(tEnd - tStart).count();
  cout << "time cost of loading genome(s):\t" << P.loadGenomeTime / 1000000000 << endl;
#endif
  if (P.pGe.gLoad == "LoadAndExit" || P.pGe.gLoad == "Remove")
  {
    return 0;
  };

  // calculate genome-related parameters
  Transcriptome *mainTranscriptome = NULL;

  mainGenome.Var = new Variation(P, mainGenome.chrStart, mainGenome.chrNameIndex);

  if (P.pGe.gFastaFiles.at(0) != "-")
  { // insert sequences in the genome
  };

  SjdbClass sjdbLoci;
  if (P.sjdbInsert.pass1)
  {
    Genome mainGenome1 = mainGenome; // not sure if I need to create the copy - mainGenome1 below should not be changed
    sjdbInsertJunctions(P, mainGenome, mainGenome1, sjdbLoci);
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////START
  if (P.runThreadN > 1)
  {
    g_threadChunks.threadArray = new pthread_t[P.runThreadN];
    pthread_mutex_init(&g_threadChunks.mutexInRead, NULL);
    pthread_mutex_init(&g_threadChunks.mutexOutSAM, NULL);
    pthread_mutex_init(&g_threadChunks.mutexOutBAM1, NULL);
    pthread_mutex_init(&g_threadChunks.mutexOutUnmappedFastx, NULL);
    pthread_mutex_init(&g_threadChunks.mutexOutFilterBySJout, NULL);
    pthread_mutex_init(&g_threadChunks.mutexStats, NULL);
    pthread_mutex_init(&g_threadChunks.mutexBAMsortBins, NULL);
    pthread_mutex_init(&g_threadChunks.mutexError, NULL);
  };
  g_statsAll.progressReportHeader(P.inOut->logProgress);
  P.bcOpt = getStBCpara(P);
  //    P.inputFileNumber=0;
  Options opt = P.bcOpt;
  // check input is a list or not
  gzFile input = gzopen(opt.transBarcodeToPos.in1.c_str(), "r");
  if (!input)
  {
    string errMsg = "Error, cannot open the file which be expected in gz or ascii format:" + opt.transBarcodeToPos.in1;
    sawErrCode(err_fileOpen_failed, errMsg);
    cerr << "Error, cannot open the file which be expected in gz or ascii format:" << opt.transBarcodeToPos.in1 << endl;
    exit(1);
  }
  char *tmpBuf = new char[1000];
  if (gzgets(input, tmpBuf, 1000) != nullptr)
  {
    string str = tmpBuf;
    if(str.size()==0){
      string errMsg="Error, failed to parse input file,"+opt.transBarcodeToPos.in1;
      sawErrCode(err_fileParse_failed,errMsg);
      exit(1);
    }
    str.erase(str.size() - 1, 1);
    ifstream tmpIf(str.c_str());
    if (!tmpIf)
    {
      // not list
      P.inputFiles.push_back(opt.transBarcodeToPos.in1);
      //		  inputFileNumber++;
    }
    else
    {
      // list
      P.inputFiles.push_back(str);
      while ((gzgets(input, tmpBuf, 1000) != nullptr))
      {
        str = tmpBuf;
        str.erase(str.size() - 1, 1);
        P.inputFiles.push_back(str);
        //			  inputFileNumber++;
      }
    }
    tmpIf.close();
  }
  else
  {
    string errMsg = "empty file or failed to get content";
    sawErrCode(err_fileIO_failed, errMsg);
    exit(1);
  }
  gzclose(input);
  delete[] tmpBuf;

  // if pe reads in, not writeFq format
  if (!opt.transBarcodeToPos.in2.empty())
  {
    P.pe = true;
    input = gzopen(opt.transBarcodeToPos.in2.c_str(), "r");
    if (!input)
    {
      string errMsg = "failed to open the file which be expected in gz or ascii format:" + opt.transBarcodeToPos.in2;
      sawErrCode(err_fileOpen_failed, errMsg);
      cerr << "Error, cannot open the file which be expected in gz or ascii format:" << opt.transBarcodeToPos.in2 << endl;
      exit(1);
    }
    tmpBuf = new char[1000];
    if (gzgets(input, tmpBuf, 1000) != nullptr)
    {
      string str = tmpBuf;

      str.erase(str.size() - 1, 1);
      ifstream tmpIf(str.c_str());
      if (!tmpIf)
      {
        // not list
        P.inputFiles2.push_back(opt.transBarcodeToPos.in2);
      }
      else
      {
        // list
        P.inputFiles2.push_back(str);
        while ((gzgets(input, tmpBuf, 1000) != nullptr))
        {
          str = tmpBuf;
          str.erase(str.size() - 1, 1);
          P.inputFiles2.push_back(str);
        }
      }
      tmpIf.close();
    }
    gzclose(input);
    delete[] tmpBuf;
    if (P.inputFiles.size() != P.inputFiles2.size())
    {
      cerr << "Error, file number are different between in1 and in2" << endl;
      exit(1);
    }
  }
  P.bcOpt.isMgzInput = true;
  for (auto &iter : P.inputFiles)
  {
    if (!mgzf_is_mgzf(iter.c_str()))
    {
      P.bcOpt.isMgzInput = false;
      break;
    }
  }
  string firstFile = P.inputFiles[0];
  //    P.inputFiles.pop_front();
  //	if(mgzf_is_mgzf(firstFile.c_str())){
  //	  P.bcOpt.isMgzInput=true;
  //	}else{
  //	  if(P.bcOpt.isMgzInput==true){
  //	    cerr<<"Warning: input file is not in mgz format"<<endl;
  //		P.bcOpt.isMgzInput=false;
  //	  }
  //	}
  if (P.bcOpt.isMgzInput)
  {
    P.inOut->logMain << timeMonthDayTime() << " ..... mgz input\n"
                     << flush;
  }
  else
  {
    P.inOut->logMain << timeMonthDayTime() << " ..... gz input\n"
                     << flush;
  }
  if (P.pe)
  {
    string secondFile = P.inputFiles2[0];
    //        P.inputFiles2.pop_front();
    if (P.inputFiles.size() != P.inputFiles2.size())
    {
      string errMsg = "number of Read1 file is not equal to that of Read2 file";
      sawErrCode(err_fileParse_failed, errMsg);
      cerr << "Error, number of Read1 file is not equal to that of Read2 file" << endl;
      exit(1);
    }
    P.reader = new FastqReaderPair(firstFile, secondFile);
  }
  else
  {
    P.reader = new FastqReaderPair(firstFile, true);
  }
  if (P.twoPass.yes)
  { // 2-pass
    // re-define P for the pass1

    Genome mainGenome1 = mainGenome;

    Parameters P1 = P;
    // turn off unnecessary calculations
    P1.outSAMtype[0] = "None";
    P1.outSAMbool = false;
    P1.outBAMunsorted = false;
    P1.outBAMcoord = false;

    P1.pCh.segmentMin = 0;

    P1.quant.yes = false;
    P1.quant.trSAM.yes = false;
    P1.quant.geCount.yes = false;
    P1.outSAMunmapped.within = false;
    P1.quant.trSAM.bamYes = false;

    P1.outFilterBySJoutStage = 0;

    P1.outReadsUnmapped = "None";

    P1.outFileNamePrefix = P.twoPass.dir;

    P1.readMapNumber = min(P.twoPass.pass1readsN, P.readMapNumber);
    //         P1.inOut->logMain.open((P1.outFileNamePrefix + "Log.out").c_str());

    P1.wasp.outputMode = "None"; // no WASP filtering on the 1st pass
    P1.pSolo.type = 0;           // no solo in the first pass

    g_statsAll.resetN();
    time(&g_statsAll.timeStartMap);
    P.inOut->logProgress << timeMonthDayTime(g_statsAll.timeStartMap) << "\tStarted 1st pass mapping\n"
                         << flush;
    *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started 1st pass mapping\n"
                        << flush;

    // run mapping for Pass1
    ReadAlignChunk *RAchunk1[P.runThreadN];
    for (int ii = 0; ii < P1.runThreadN; ii++)
    {
      RAchunk1[ii] = new ReadAlignChunk(P1, mainGenome, mainTranscriptome, ii);
    };
    mapThreadsSpawn(P1, RAchunk1);
    outputSJ(RAchunk1, P1); // collapse and output junctions
    //         for (int ii=0;ii<P1.runThreadN;ii++) {
    //             delete [] RAchunk[ii];
    //         };

    time_t rawtime;
    time(&rawtime);
    P.inOut->logProgress << timeMonthDayTime(rawtime) << "\tFinished 1st pass mapping\n";
    *P.inOut->logStdOut << timeMonthDayTime(rawtime) << " ..... finished 1st pass mapping\n"
                        << flush;
    ofstream logFinal1((P.twoPass.dir + "/Log.final.out").c_str());
    g_statsAll.reportFinal(logFinal1);

    P.twoPass.pass2 = true; // starting the 2nd pass
    P.twoPass.pass1sjFile = P.twoPass.dir + "/SJ.out.tab";

    sjdbInsertJunctions(P, mainGenome, mainGenome1, sjdbLoci);

    // reopen reads files
    P.closeReadsFiles();
    P.openReadsFiles();
    delete P.reader;
    P.reader = new FastqReaderPair(opt.transBarcodeToPos.in1, opt.transBarcodeToPos.in2, true);
  }
  else
  { // not 2-pass
    // nothing for now
  };

  if (P.quant.yes)
  { // load transcriptome
    mainTranscriptome = new Transcriptome(P);
  };
  P.gzfileIdx = 1;
  // initialize Stats
  g_statsAll.resetN();
  time(&g_statsAll.timeStartMap);
  *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeStartMap) << " ..... started mapping\n"
                      << flush;

  g_statsAll.timeLastReport = g_statsAll.timeStartMap;

  // open SAM/BAM files for output
  if (P.outSAMmode != "None")
  { // open SAM file and write header
    ostringstream samHeaderStream;

    for (uint ii = 0; ii < mainGenome.nChrReal; ii++)
    {
      samHeaderStream << "@SQ\tSN:" << mainGenome.chrName.at(ii) << "\tLN:" << mainGenome.chrLength[ii] << "\n";
    };

    mainGenome.chrNameAll = mainGenome.chrName;
    mainGenome.chrLengthAll = mainGenome.chrLength;
    { // add exra references
      ifstream extrastream(P.pGe.gDir + "/extraReferences.txt");
      while (extrastream.good())
      {
        string line1;
        getline(extrastream, line1);
        istringstream stream1(line1);
        string field1;
        stream1 >> field1; // should check for @SQ

        if (field1 != "")
        { // skip blank lines
          samHeaderStream << line1 << "\n";

          stream1 >> field1;
          mainGenome.chrNameAll.push_back(field1.substr(3));
          stream1 >> field1;
          mainGenome.chrLengthAll.push_back((uint)stoll(field1.substr(3)));
        };
      };
      extrastream.close();
    };

    if (P.outSAMheaderPG.at(0) != "-")
    {
      samHeaderStream << P.outSAMheaderPG.at(0);
      for (uint ii = 1; ii < P.outSAMheaderPG.size(); ii++)
      {
        samHeaderStream << "\t" << P.outSAMheaderPG.at(ii);
      };
      samHeaderStream << "\n";
    };

    samHeaderStream << "@PG\tID:BCSTAR\tPN:BCSTAR\tVN:" << STAR_VERSION << "\tCL:" << P.commandLineFull << "\n";

    if (P.outSAMheaderCommentFile != "-")
    {
      ifstream comstream(P.outSAMheaderCommentFile);
      while (comstream.good())
      {
        string line1;
        getline(comstream, line1);
        if (line1.find_first_not_of(" \t\n\v\f\r") != std::string::npos)
        { // skip blank lines
          samHeaderStream << line1 << "\n";
        };
      };
      comstream.close();
    };

    for (uint32 ii = 0; ii < P.outSAMattrRGlineSplit.size(); ii++)
    { //@RG lines
      samHeaderStream << "@RG\t" << P.outSAMattrRGlineSplit.at(ii) << "\n";
    };

    samHeaderStream << "@CO\t"
                    << "user command line: " << P.commandLine << "\n";

    samHeaderStream << P.samHeaderExtra;

    if (P.outSAMheaderHD.at(0) != "-")
    {
      P.samHeaderHD = P.outSAMheaderHD.at(0);
      for (uint ii = 1; ii < P.outSAMheaderHD.size(); ii++)
      {
        P.samHeaderHD += "\t" + P.outSAMheaderHD.at(ii);
      };
    }
    else
    {
      P.samHeaderHD = "@HD\tVN:1.4";
    };

    P.samHeader = P.samHeaderHD + "\n" + samHeaderStream.str();
    // for the sorted BAM, need to add SO:cooridnate to the header line
    P.samHeaderSortedCoord = P.samHeaderHD + (P.outSAMheaderHD.size() == 0 ? "" : "\tSO:coordinate") + "\n" + samHeaderStream.str();
    // add ACGT encode rule
    if (P.outSAMattributes.at(0) == "spatial")
    {
      P.samHeaderSortedCoord += "@CO\tUR is the umi(MID) sequence which was encoded to a hex number, encode rule:A->0 C->1 G->2 T->3\n";
    }
    if (P.outSAMbool)
    { //
      *P.inOut->outSAM << P.samHeader;
    };
    if (P.outBAMunsorted)
    {
      outBAMwriteHeader(P.inOut->outBAMfileUnsorted, P.samHeader, mainGenome.chrNameAll, mainGenome.chrLengthAll);
    };
    //             if (P.outBAMcoord){
    //                 outBAMwriteHeader(P.inOut->outBAMfileCoord,P.samHeader,mainGenome.chrName,mainGenome.chrLength);
    //             };

    if (P.quant.trSAM.bamYes)
    {
      samHeaderStream.str("");
      vector<uint> trlength;
      for (uint32 ii = 0; ii < mainTranscriptome->trID.size(); ii++)
      {
        uint32 iex1 = mainTranscriptome->trExI[ii] + mainTranscriptome->trExN[ii] - 1; // last exon of the transcript
        trlength.push_back(mainTranscriptome->exLenCum[iex1] + mainTranscriptome->exSE[2 * iex1 + 1] - mainTranscriptome->exSE[2 * iex1] + 1);
        samHeaderStream << "@SQ\tSN:" << mainTranscriptome->trID.at(ii) << "\tLN:" << trlength.back() << "\n";
      };
      for (uint32 ii = 0; ii < P.outSAMattrRGlineSplit.size(); ii++)
      { //@RG lines
        samHeaderStream << "@RG\t" << P.outSAMattrRGlineSplit.at(ii) << "\n";
      };
      outBAMwriteHeader(P.inOut->outQuantBAMfile, samHeaderStream.str(), mainTranscriptome->trID, trlength);
    };
  };

  // initialize chimeric parameters here

  // P.inOut->logMain << "mlock value="<<mlockall(MCL_CURRENT|MCL_FUTURE) <<"\n"<<flush;

  // prepare chunks and spawn mapping threads
  P.inOut->logMain << timeMonthDayTime() << " ..... started barcode mapping and star mapping\n"
                   << flush;
  ReadAlignChunk *RAchunk[P.runThreadN];
  for (int ii = 0; ii < P.runThreadN; ii++)
  {
    RAchunk[ii] = new ReadAlignChunk(P, mainGenome, mainTranscriptome, ii);
  };

  mapThreadsSpawn(P, RAchunk);
  P.inOut->logMain << timeMonthDayTime() << " ..... finished barcode mapping and star mapping\n"
                   << flush;
  if (P.outFilterBySJoutStage == 1)
  { // completed stage 1, go to stage 2
    P.inOut->logMain << "Completed stage 1 mapping of outFilterBySJout mapping\n"
                     << flush;
    outputSJ(RAchunk, P); // collapse novel junctions
    P.readFilesIndex = -1;

    P.outFilterBySJoutStage = 2;
    if (P.outBAMcoord)
    {
      for (int it = 0; it < P.runThreadN; it++)
      { // prepare the unmapped bin
        RAchunk[it]->chunkOutBAMcoord->coordUnmappedPrepareBySJout();
      };
    };

    mapThreadsSpawn(P, RAchunk);
  };

  // close some BAM files
  if (P.inOut->outBAMfileUnsorted != NULL)
  {
    if (bgzf_flush(P.inOut->outBAMfileUnsorted) == -1)
    {

      string errMsg = "bgzf_flush error";
      sawErrCode(err_fileIO_failed, errMsg);
      cerr << "Error, bgzf_flush error" << endl;
      exit(1);
    }
    bgzf_close(P.inOut->outBAMfileUnsorted);
  };
  if (P.inOut->outQuantBAMfile != NULL)
  {
    if (bgzf_flush(P.inOut->outQuantBAMfile) == -1)
    {
      string errMsg = "bgzf_flush error";
      sawErrCode(err_fileIO_failed, errMsg);
      cerr << "Error, bgzf_flush error" << endl;
      exit(1);
    }
    bgzf_close(P.inOut->outQuantBAMfile);
  };

  if (P.outBAMcoord && P.limitBAMsortRAM == 0)
  { // make it equal ot the genome size
    P.limitBAMsortRAM = mainGenome.nGenome + mainGenome.SA.lengthByte + mainGenome.SAi.lengthByte;
  };

  time(&g_statsAll.timeFinishMap);
  *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinishMap) << " ..... finished mapping\n"
                      << flush;

  // no need for genome anymore, free the memory
#ifdef D_PRINT
  cout << "before freeGenome:\t" << (float)get_proc_virtualMem(getpid()) / (1024 * 1024) << "\t" << (float)get_proc_VmRSS(getpid()) / (1024 * 1024) << endl;
#endif
  mainGenome.freeMemory();
#ifdef D_PRINT
  cout << "after freeGenome:\t" << (float)get_proc_virtualMem(getpid()) / (1024 * 1024) << "\t" << (float)get_proc_VmRSS(getpid()) / (1024 * 1024) << endl;
#endif
  // aggregate output junctions
  // collapse splice junctions from different threads/chunks, and output them
  outputSJ(RAchunk, P);

  // solo genes
  Solo soloMain(RAchunk, P, *RAchunk[0]->chunkTr); // solo for genes
  soloMain.processAndOutput();

  if (P.quant.geCount.yes)
  { // output gene quantifications
    for (int ichunk = 1; ichunk < P.runThreadN; ichunk++)
    { // sum counts from all chunks into 0th chunk
      RAchunk[0]->chunkTr->quants->addQuants(*(RAchunk[ichunk]->chunkTr->quants));
    };
    RAchunk[0]->chunkTr->quantsOutput();
  };

  if (P.runThreadN > 1 && P.outSAMorder == "PairedKeepInputOrder")
  { // concatenate Aligned.* files
    RAchunk[0]->chunkFilesCat(P.inOut->outSAM, P.outFileTmp + "/Aligned.out.sam.chunk", g_threadChunks.chunkOutN);
  };
#ifdef TIME_STAT
  tStart = chrono::high_resolution_clock::now();
#endif
  if (P.outBAMcoord)
  { // sort BAM if needed
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... started sorting BAM\n"
                        << flush;
    P.inOut->logMain << timeMonthDayTime() << " ..... started sorting BAM\n"
                     << flush;
    uint32 nBins = P.outBAMcoordNbins;

    // check max size needed for sorting
    uint maxMem = 0;
    for (uint32 ibin = 0; ibin < nBins - 1; ibin++)
    { // check all bins
      uint binS = 0;
      for (int it = 0; it < P.runThreadN; it++)
      { // collect sizes from threads
        binS += RAchunk[it]->chunkOutBAMcoord->binTotalUnCompressBytes[ibin] + 24 * RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
      };
      if (binS > maxMem)
        maxMem = binS;
    };
    P.inOut->logMain << "Max memory needed for sorting = " << maxMem << endl;
    uint availableMem = 0;
    uint vmPeak = get_proc_virtualMemPeak(getpid()) * 1024;
    uint vmSize = get_proc_virtualMem(getpid()) * 1024;
#ifdef D_PRINT
    uint vmRss = get_proc_VmRSS(getpid()) * 1024;
#endif
    availableMem = vmPeak - vmSize;
#ifdef D_PRINT
    cout << "availableMem\t" << (float)availableMem / (1024 * 1024 * 1024) << "\t" << (float)vmPeak / (1024 * 1024 * 1024) << "\t" << (float)vmSize / (1024 * 1024 * 1024) << "\t" << (float)vmRss / (1024 * 1024 * 1024) << endl;
#endif
    if (maxMem < availableMem)
    {
      P.limitBAMsortRAM = availableMem;
      //	  	}else {
      //		  if(maxMem > P.limitBAMsortRAM) {
      //			P.limitBAMsortRAM = availableMem;
      //		  }
    }
    if (maxMem > P.limitBAMsortRAM)
    {
      ostringstream errOut;
      errOut << "EXITING because of fatal ERROR: not enough memory for BAM sorting: \n";
      errOut << "SOLUTION: re-run BCSTAR with at least --limitBAMsortRAM " << maxMem + 1000000000;
      sawErrCode(err_resource_oom, errOut.str());
      exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }
    else if (maxMem == 0)
    {
      P.inOut->logMain << "WARNING: nothing to sort - no output alignments" << endl;
      BGZF *bgzfOut;
      bgzfOut = bgzf_open(P.outBAMfileCoordName.c_str(), ("w" + to_string((long long)P.outBAMcompression)).c_str());
      if (bgzfOut == NULL)
      {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open output bam file: " << P.outBAMfileCoordName << "\n";
        errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running BCSTAR";
        sawErrCode(err_fileOpen_failed, errOut.str());
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
      };
      outBAMwriteHeader(bgzfOut, P.samHeaderSortedCoord, mainGenome.chrNameAll, mainGenome.chrLengthAll);
      bgzf_close(bgzfOut);
    }
    else
    { // sort
      uint totalMem = 0;
      //            #pragma omp parallel num_threads(P.outBAMsortingThreadNactual)
#pragma omp parallel num_threads(P.runThreadN)
#pragma omp for schedule(static)
      //            #pragma omp for schedule (dynamic,1)
      for (uint32 ibin1 = 0; ibin1 < nBins; ibin1++)
      {
        uint32 ibin = nBins - 1 - ibin1; // reverse order to start with the last bin - unmapped reads

        uint binN = 0, binS = 0;
#if defined(TESTCOMPRESS) || defined(TESTZSTD)
        uint64_t *uncompressSizeArr = new uint64_t[P.runThreadN];
#endif
#pragma omp critical
        {
#ifdef D_PRINT
          P.inOut->logProgress << omp_get_thread_num() << "\t" << ibin << "\t";
#endif
          //					  }

          //                cout<<ibin1<<"\t:\t";

          for (int it = 0; it < P.runThreadN; it++)
          { // collect sizes from threads
            uncompressSizeArr[it] = RAchunk[it]->chunkOutBAMcoord->binTotalUnCompressBytes[ibin];
            binN += RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin];
            binS += RAchunk[it]->chunkOutBAMcoord->binTotalUnCompressBytes[ibin];
#ifdef D_PRINT
            P.inOut->logProgress << ibin << "\t" << RAchunk[it]->chunkOutBAMcoord->binTotalBytes[ibin] << "\t";
#endif
            //					  }
          };
#ifdef D_PRINT
          P.inOut->logProgress << binS << endl;
#endif
        }
        //				cout<<endl;
        if (binS == 0)
          continue; // empty bin

        if (ibin == nBins - 1)
        { // last bin for unmapped reads
          BAMbinSortUnmapped(ibin, P.runThreadN, P.outBAMsortTmpDir, P, mainGenome);
        }
        else
        {
          uint newMem = binS + binN * 24;
          bool boolWait = true;
          while (boolWait)
          {
#pragma omp critical
            {
              if (totalMem + newMem < P.limitBAMsortRAM)
              {
                boolWait = false;
                totalMem += newMem;
              }
            }
            usleep(100000);
          };
#if defined(TESTCOMPRESS) || defined(TESTZSTD)
          BAMbinSortByCoordinate(ibin, binN, binS, P.runThreadN, P.outBAMsortTmpDir, P, mainGenome, uncompressSizeArr);
          delete[] uncompressSizeArr;
#else
          BAMbinSortByCoordinate(ibin, binN, binS, P.runThreadN, P.outBAMsortTmpDir, P, mainGenome);
#endif
#pragma omp critical
          totalMem -= newMem; //"release" RAM
        };
      };

      // concatenate all BAM files, using bam_cat
      char **bamBinNames = new char *[nBins];
      vector<string> bamBinNamesV;
      for (uint32 ibin = 0; ibin < nBins; ibin++)
      {
        bamBinNamesV.push_back(P.outBAMsortTmpDir + "/b" + std::to_string((uint)ibin));
        struct stat buffer;
        if (stat(bamBinNamesV.back().c_str(), &buffer) != 0)
        { // check if file exists
          bamBinNamesV.pop_back();
        };
      };
      for (uint32 ibin = 0; ibin < bamBinNamesV.size(); ibin++)
      {
        bamBinNames[ibin] = (char *)bamBinNamesV.at(ibin).c_str();
      };
      bam_cat(bamBinNamesV.size(), bamBinNames, 0, P.outBAMfileCoordName.c_str());
    };
  };
#ifdef TIME_STAT
  tEnd = chrono::high_resolution_clock::now();
  P.sortBamTime = chrono::duration<double, nano>(tEnd - tStart).count();
  cout << "time cost of sorting BAM(s):\t" << P.sortBamTime / 1000000000 << endl;
#endif
  // wiggle output
  if (P.outWigFlags.yes)
  {
    *(P.inOut->logStdOut) << timeMonthDayTime() << " ..... started wiggle output\n"
                          << flush;
    P.inOut->logMain << timeMonthDayTime() << " ..... started wiggle output\n"
                     << flush;
    string wigOutFileNamePrefix = P.outFileNamePrefix + "Signal";
    signalFromBAM(P.outBAMfileCoordName, wigOutFileNamePrefix, P);
  };

  g_statsAll.writeLines(P.inOut->outChimJunction, P.pCh.outJunctionFormat, "#", STAR_VERSION + string("   ") + P.commandLine);

  g_statsAll.progressReport(P.inOut->logProgress);
  P.inOut->logProgress << "ALL DONE!\n"
                       << flush;
  P.inOut->logFinal.open((P.outFileNamePrefix + "Log.final.out").c_str());
  g_statsAll.reportFinal(P.inOut->logFinal);

  P.closeReadsFiles(); // this will kill the readFilesCommand processes if necessary
  delete P.reader;
  // mainGenome.~Genome(); //need explicit call because of the 'delete P.inOut' below, which will destroy P.inOut->logStdOut
  if (mainGenome.sharedMemory != NULL)
  { // need explicit call because this destructor will write to files which are deleted by 'delete P.inOut' below
    delete mainGenome.sharedMemory;
    mainGenome.sharedMemory = NULL;
  };
  if (P.outTmpKeep == "None")
  {
    P.inOut->logMain << "remove tmp dir:\t" << P.outFileTmp << endl;
    sysRemoveDir(P.outFileTmp);
    FILE *tmpDirFile = fopen(P.outFileTmp.c_str(), "r");
    if (tmpDirFile != nullptr)
    {
      P.inOut->logMain << "Warning, cannot remove tmp directory. Try again..." << endl;
      fclose(tmpDirFile);
      string cmd = "rm -rf " + P.outFileTmp;
      if (system(cmd.c_str()) == -1)
      {
        string errMsg = "fail to remove tmp directory";
        sawErrCode(err_fileRemove_failed, errMsg);
        cerr << "Error, fail to remove tmp directory." << endl;
        exit(1);
      }
    }
    else
    {
      P.inOut->logMain << "Warning, cannot open such file," << P.outFileTmp << endl;
    }
    string cmd = "rm -rf " + P.outFileTmp;
    system(cmd.c_str());
  }
  else
  {
    P.inOut->logMain << "don't remove tmp dir\n"
                     << flush;
  }
  // create .bai for bam
  P.inOut->logMain << timeMonthDayTime() << " ..... started building index for bam!\n"
                   << flush;
  if (P.outBAMcoord)
  {
    if (bam_index_build(P.outBAMfileCoordName.c_str(), 16) != 0)
    {
      string errMsg = "bam_index_build failed";
      sawErrCode(err_indexBuildeForBAM_failed, errMsg);
      cerr << "bam_index_build failed: " << P.outBAMfileCoordName << endl;
      exit(1);
    }
  }
  P.inOut->logMain << timeMonthDayTime() << " ..... finished building index!\n"
                   << flush;
  *P.inOut->logStdOut << timeMonthDayTime(g_statsAll.timeFinish) << " ..... finished successfully\n"
                      << flush;
  P.inOut->logMain << "ALL DONE!\n"
                   << flush;
  delete P.inOut; // to close files
  return 0;
};

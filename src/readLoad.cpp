#include "readLoad.h"
#include "ErrorWarning.h"
#include "util.h"

int readLoad(istream& readInStream, Parameters& P, uint iMate, uint& Lread, uint& LreadOriginal, char* readName, char* Seq, char* SeqNum, char* Qual, char* QualNum, uint &clip3pNtotal, uint &clip5pNtotal, uint &clip3pAdapterN, uint &iReadAll, uint &readFilesIndex, char &readFilter, string &readNameExtra,uint32_t &xCor,uint32_t &yCor,string &umiStr){
    //load one read from a stream
    int readFileType=0;

//     readInStream.getline(readName,DEF_readNameLengthMax); //extract name

    if (readInStream.peek()!='@' && readInStream.peek()!='>') return -1; //end of the stream

    readName[0]=0;//clear char array
  if(P.outSAMattributes.at(0)=="spatial") {
	string rawReadName;
	readInStream >> rawReadName; //TODO check that it does not overflow the array
 
	vector<string> eles;
 
	split(rawReadName, eles, "|||");
	if(eles.size()!=3) {
        string errMsg=readName;
        errMsg="read name format error:" + errMsg;
	  sawErrCode(err_fileParse_failed,errMsg);
	  cerr << "Error,read name format error:" << readName << endl;
	  exit(1);
	}
	istringstream tmps(eles[0]);
	tmps >> readName;
	eles[1] = eles[1].erase(0, 5);
	eles[2] = eles[2].erase(0, 5);
	vector<string> eles2;
	split(eles[1], eles2, "_");
	xCor = atoi(eles2[0].c_str());
	yCor = atoi(eles2[1].c_str());
	umiStr = eles[2];
  }else{
	readInStream>>readName;
  }
    if (strlen(readName)>=DEF_readNameLengthMax-1) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in reads input: read name is too long:" << readInStream.gcount()<<"\n";
        errOut << "Read Name="<<readName<<"\n";
        errOut << "DEF_readNameLengthMax="<<DEF_readNameLengthMax<<"\n";
        errOut << "SOLUTION: increase DEF_readNameLengthMax in IncludeDefine.h and re-compile BCSTAR\n";
        sawErrCode(err_sw_exception,errOut.str());
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    readInStream >> iReadAll >> readFilter >> readFilesIndex; //extract read number

    getline(readInStream, readNameExtra);
    if (!readNameExtra.empty()) {
        size_t n1=readNameExtra.find_first_not_of(" \t");
        if (n1!=std::string::npos) {
            readNameExtra=readNameExtra.substr(n1);
        } else {
            readNameExtra="";
        };
    };

//     readInStream.ignore(DEF_readNameSeqLengthMax,'\n');//ignore the resit of the line - just in case

    readInStream.getline(Seq,DEF_readSeqLengthMax+1); //extract sequence

    Lread=0;
    for (int ii=0; ii<readInStream.gcount()-1; ii++) {
        if (int(Seq[ii])>=32) {
            Seq[Lread]=Seq[ii];
            ++Lread;
        };
    };

    if (Lread<1) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in reads input: short read sequence line: " << Lread <<"\n";
        errOut << "Read Name="<<readName<<"\n";
        errOut << "Read Sequence="<<Seq<<"===\n";
        errOut << "DEF_readNameLengthMax="<<DEF_readNameLengthMax<<"\n";
        errOut << "DEF_readSeqLengthMax="<<DEF_readSeqLengthMax<<"\n";
        sawErrCode(err_sw_exception,errOut.str());
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    if (Lread>DEF_readSeqLengthMax) {
        ostringstream errOut;
        errOut << "EXITING because of FATAL ERROR in reads input: Lread>=" << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<"\n";
        errOut << "Read Name="<<readName<<"\n";
        errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile BCSTAR\n";
        sawErrCode(err_sw_exception,errOut.str());
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

//     //was trying to read multi-line
//     char nextChar='A';
//     Lread=0;
//     while (nextChar!='@' && nextChar!='>' && nextChar!='+' && nextChar!=' ' && nextChar!='\n' && !readInStream.eof()) {//read multi-line fasta
//         readInStream.getline(Seq+Lread,DEF_readSeqLengthMax+1); //extract sequence
//         Lread+=(uint) readInStream.gcount() - 1;    //count chars in the sequence line, but do not read yet
//         nextChar=readInStream.peek();
//     };
//     if (Lread>DEF_readSeqLengthMax) {
//         ostringstream errOut;
//         errOut << "EXITING because of FATAL ERROR in reads input: Lread>=" << Lread << "   while DEF_readSeqLengthMax=" << DEF_readSeqLengthMax <<"\n";
//         errOut << "Read Name="<<readName<<"\n";
//         errOut << "SOLUTION: increase DEF_readSeqLengthMax in IncludeDefine.h and re-compile BCSTAR\n";
//         exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
//     };
//     LreadOriginal=Lread;

    LreadOriginal=Lread;
    if ( Lread>(P.clip5pNbases[iMate]+P.clip3pNbases[iMate]) ) {
        Lread=Lread-(P.clip5pNbases[iMate]+P.clip3pNbases[iMate]);
    } else {
        Lread=0;
    };
    convertNucleotidesToNumbersRemoveControls(Seq+P.clip5pNbases[iMate],SeqNum,Lread);

    //clip the adapter
    if (P.clip3pAdapterSeq.at(iMate).length()>0) {
        clip3pAdapterN = Lread-localSearch(SeqNum,Lread,P.clip3pAdapterSeqNum[iMate],P.clip3pAdapterSeq.at(iMate).length(),P.clip3pAdapterMMp[iMate]);
        Lread = Lread>clip3pAdapterN ? Lread-clip3pAdapterN : 0;
    } else {
        clip3pAdapterN = 0;
    };

    //final read length, trim 3p after the adapter was clipped
    if (Lread>P.clip3pAfterAdapterNbases[iMate]) {
        Lread =Lread - P.clip3pAfterAdapterNbases[iMate];
    } else {
        Lread=0;
    };

    clip3pNtotal=P.clip3pNbases[iMate] + clip3pAdapterN + P.clip3pAfterAdapterNbases[iMate];
    clip5pNtotal=P.clip5pNbases[iMate];

    if (readName[0]=='@') {//fastq format, read qualities
        readFileType=2;
        readInStream.ignore(DEF_readNameLengthMax,'\n'); //extract header line
        readInStream.getline(Qual,DEF_readSeqLengthMax);//read qualities
        if ((uint) readInStream.gcount() != LreadOriginal+1) {//inconsistent read sequence and quality
            ostringstream errOut;
            errOut << "EXITING because of FATAL ERROR in reads input: quality string length is not equal to sequence length\n";
            errOut << readName<<"\n";
            errOut << Seq <<"\n";
            errOut << Qual <<"\n";
            errOut << "SOLUTION: fix your fastq file\n";
            sawErrCode(err_fileParse_failed,errOut.str());
            exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
        };
        if (P.outQSconversionAdd!=0) {
            for (uint ii=0;ii<LreadOriginal;ii++) {
                int qs=int(Qual[ii])+P.outQSconversionAdd;
                if (qs<33) {
                    qs=33;
                } else if (qs>126) {
                    qs=126;
                };
                Qual[ii]=qs;
            };
        };

    } else if (readName[0]=='>') {//fasta format, assign Qtop to all qualities
        readFileType=1;
        for (uint ii=0;ii<LreadOriginal;ii++) Qual[ii]='A';
        Qual[LreadOriginal]=0;
    } else {//header
        ostringstream errOut;
        errOut <<"Unknown reads file format: header line does not start with @ or > : "<< readName<<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    for (uint ii=0;ii<Lread;ii++) {//for now: qualities are all 1
        if (SeqNum[ii]<4) {
            QualNum[ii]=1;
        } else {
            QualNum[ii]=0;
        };
    };

//     for (uint ii=0;ii<Lread;ii++) {//simply cut too high Qs
//         QualNum[ii]=(Qual[ii+P.clip5pNbases[iMate]] > P.QasciiSubtract) ? (Qual[ii+P.clip5pNbases[iMate]] - P.QasciiSubtract) : 0; //substract QasciiSubtract
//         QualNum[ii]=P.QSconv[(int) QualNum[ii]];
//         QualNum[ii]=min(QualNum[ii], P.Qtop);//cut QSs at the Qtop
// //         if (QualNum[ii]==2) QualNum[ii]=P.Qtop;
//         if (SeqNum[ii]>3) QualNum[ii]=0; //QS=0 for Ns
//         Qual1[1][Lread-ii-1]=QualNum[ii]; //reverse
//     };


    //trim read name
    for (uint ii=0; ii<P.readNameSeparatorChar.size(); ii++)
    {
        char* pSlash=strchr(readName,P.readNameSeparatorChar.at(ii)); //trim everything after ' '
        if (pSlash!=NULL) *pSlash=0;
    };
    return readFileType;
};

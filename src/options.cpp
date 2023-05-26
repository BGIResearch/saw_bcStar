/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "options.h"

Options::Options()
{
}

void Options::init()
{
	
	isSeq500 = getIsSeq500(platform);
	if (!drawHeatMap.q10TiffFile.empty()){
		drawHeatMap.getQ10Tiff = true;
	}
	setFovRange(drawHeatMap.fovRange);
}

bool Options::validate()
{
	if (in.empty() || out.empty()) {
		string errMsg="please give input file and output file";
	  sawErrCode(err_param_unset_bcmap,errMsg);
		cerr << "please give input file and output file" << endl;
		exit(-1);
		//return false;
	}
	if (actionInt == 4) {
//		if (transBarcodeToPos.in1.empty() || transBarcodeToPos.in2.empty()) {
		if (transBarcodeToPos.in1.empty()) {
			string errMsg="please give fastq files";
		  sawErrCode(err_param_unset_bcmap,errMsg);
			cerr << "please give fastq files";
			exit(-1);
			//return false;
		}
		check_file_valid(transBarcodeToPos.in1);
		if(!transBarcodeToPos.in2.empty()){
		  gzFile gzin1=gzopen(transBarcodeToPos.in1.c_str(),"r");
		  char* lineBuf=new char[500];
		  for(int i=0;i<4;i++) {
			if(gzgets(gzin1, lineBuf, 500)!=NULL) {
			  string line=lineBuf;
			  line=line.erase(line.size()-1,1);
			  if(i==1){
//			    cout<<line<<endl;
//			    exit(1);
				if(uint(barcodeStart+barcodeLen)>line.size()){
					string errMsg="please check the barcode position and length";
				  sawErrCode(err_param_invalid_bcmap,errMsg);
				  cerr<<"Error, please check the barcode position and length"<<endl;
				  exit(1);
				}
				if(uint(transBarcodeToPos.umiStart+transBarcodeToPos.umiLen)>line.size()){
					string errMsg="please check the umi position and length";
				  sawErrCode(err_param_invalid_bcmap,errMsg);
				  cerr<<"Error, please check the umi position and length"<<endl;
				  exit(1);
				}
			  }
			}
		  }
		  
		}
	 
//		check_file_valid(transBarcodeToPos.in2);
	}
	if (! maskFile.empty()){
		check_file_valid(maskFile);
	}

	if (barcodeSegment<=0){
		cerr << "barcodeSegment should >0, but get: " << barcodeSegment << ". set to be the default value 1"<<endl;
		barcodeSegment = 1;
		barcodeStat.segment = 1;
	}
	
	return true;
}

int Options::transRC(string isRC){
	int rc = 0;
	transform(isRC.begin(), isRC.end(), isRC.begin(), towlower);
	if (isRC.compare("false")==0){
		rc = 0;
	}else if (isRC.compare("true")==0){
		rc = 1;
	}else if (isRC.compare("all")==0){
		rc = 2;
	}else{
//		loginfo("RC option should one of [true, false, all], but got: " + isRC);
string errMsg="RC option should one of [true, false, all], but got: "+isRC;
	  sawErrCode(err_param_invalid_bcmap,errMsg);
		cerr<<"Error, RC option should one of [true, false, all], but got: "<<isRC<<endl;
		exit(1);
	}
	return rc;
}

bool Options::getIsSeq500(string& platform){
	if (platform.compare("SEQ500") == 0 || platform.compare("seq500") == 0 || platform.compare("Seq500") == 0 || platform.compare("SEq500") == 0 || platform.compare("SeQ500") == 0){
		return true;
	}else{
		return false;
	}
}

void Options::setFovRange(string fovRange){
	vector<string> ranges;
	split(fovRange, ranges, "_");
	if (ranges.size() < 2){
//		loginfo("fov format is wrong");
string errMsg="fov format is wrong";
	  sawErrCode(err_param_invalid_bcmap,errMsg);
		cerr<<"Error, fov format is wrong"<<endl;
		exit(1);
	}
	vector<string> colRange;
	vector<string> rowRange;
	split(ranges.at(0), colRange, "-");
	split(ranges.at(1), rowRange, "-");
	if (colRange.size() < 2 || rowRange.size() < 2){
//		loginfo("fov format is wrong");
	  string errMsg="fov format is wrong";
	  sawErrCode(err_param_invalid_bcmap,errMsg);
	  cerr<<"Error, fov format is wrong"<<endl;
		exit(1);
	}
	drawHeatMap.minCol = stoi(colRange.at(0));
	drawHeatMap.maxCol = stoi(colRange.at(1));
	drawHeatMap.minRow = stoi(rowRange.at(0));
	drawHeatMap.maxRow = stoi(rowRange.at(1));
}

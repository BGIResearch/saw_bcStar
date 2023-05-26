/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cctype>
#include <algorithm>
#include "util.h"
#include "cmdline.h"

using namespace std;

class DrawHeatMapOptions {
public:
	DrawHeatMapOptions() {
	}

public:
	//hashmap file of the overlaped barcode list
	string map;
	//mask file path
	string maskFile;
	string fovRange;
	//min fov colimn
	int minCol;
	// max fov column
	int maxCol;
	//min fov row
	int minRow;
	//max fov row
	int maxRow;	
	//wether generate q10 heatmap tiff
	bool getQ10Tiff = false;
	string q10TiffFile;
};

class BarcodeOverlapOptions {
public:
	BarcodeOverlapOptions() {

	}
public:
	// map file contianning the barcode hash map of second sequencing
	string in2;
	int mismatch;
	
	// number of base trimed on the front of sequence
	//int frontTrim = 0;
	// number of base trimed on the tail of sequence
	//int tailTrim = 0;
	// map bucket size
};

class BarcodeStatOptions {
public:
	BarcodeStatOptions(){
	}
public:
	int segment;
	// wheather the barcodes of two sequencing are reverse complement
	string rcString;
	int rc;
	string readidSep;
};

class TransBarcodeToPosOptions {
public:
	TransBarcodeToPosOptions() {

	}
public:
	//first sequencing fastq file or barcode map file
	string in;
	//second sequencing fastq file or barcode map file of read1
	string in1;
	//second sequencing fastq file of read2
	string in2;
	// second sequencing output fastq file of read1
	//string out1;
	//second sequencing output fastq file of read2
	//string out2;
	//allowed max mismatch
	int mismatch;
	//barcode to position map dump file path
	//string bpMapOutFile;
	//file path for reads with unmapped barcode
	string unmappedOutFile;
	//which reads contains the umi
	int umiRead;
	//umi start position
	int umiStart;
	//umi length
	int umiLen;
	//mapped dnb list file
	string mappedDNBOutFile;
	//fixed sequence that will be detected in the read1.
	string fixedSequence;
	//fixed sequence file contianing fixed sequence
	string fixedSequenceFile;
	//fixed sequence start position
	int fixedStart;
};

class Options {
public:
	Options();
	void init();
	bool validate();
	int transRC(string isRC);
	bool getIsSeq500(string& platform);
	void setFovRange(string fovRange);
public:
	enum actions { barcode_stat = 1, barcode_overlap = 2, get_barcode_position_map = 3, map_barcode_to_slide = 4, merge_barcode_list = 5 } action;
	int actionInt = 1;
	// file name of first sequencing read
	string in;
	// output file name
	string out;
	// encode rule
	string encodeRule;
	// report output file path
	string report;
	//compression level for gzip output
	int compression;
	int barcodeLen;
	int barcodeStart;
	int barcodeSegment;
	int turnFovDegree;
	string platform;
	string chipID;
	bool isSeq500;
	string maskFile;
	long mapSize = 100000000;
	bool verbose;
	bool useF14=true;
	bool useBf;
	int thread;
	uint64_t bcNum;
	float bcMappingRateCutoff=0.0;
	int mismatchInPolyA=0;
	int polyAnum=0;
	string rcString;
	int rc;
	bool splitBarcode=false;
	DrawHeatMapOptions drawHeatMap;
	BarcodeOverlapOptions barcodeOverlap;
	BarcodeStatOptions barcodeStat;
	TransBarcodeToPosOptions transBarcodeToPos;
  	cmdline::parser cmd;
	//add mgz support
	bool isMgzInput=false;
};

#endif

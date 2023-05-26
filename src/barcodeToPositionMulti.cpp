/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "barcodeToPositionMulti.h"

BarcodeToPositionMulti::BarcodeToPositionMulti(Options* opt)
{
	mOptions = opt;
	mProduceFinished = false;
	mFinishedThreads = 0;
	mOutStream = NULL;
	mZipFile = NULL;
	mWriter = NULL;
	mUnmappedWriter = NULL;
//	bool isSeq500 = opt->isSeq500;
	mbpmap = new BarcodePositionMap(opt);
	
}

BarcodeToPositionMulti::~BarcodeToPositionMulti()
{
}


void BarcodeToPositionMulti::initOutput() {
	mWriter = new WriterThread(mOptions->out, mOptions->compression);
	if (!mOptions->transBarcodeToPos.unmappedOutFile.empty()) {
		mUnmappedWriter = new WriterThread(mOptions->transBarcodeToPos.unmappedOutFile, mOptions->compression);
	}
}

void BarcodeToPositionMulti::closeOutput()
{
	if (mWriter) {
		delete mWriter;
		mWriter = NULL;
	}
	if (mUnmappedWriter) {
		delete mUnmappedWriter;
		mUnmappedWriter = NULL;
	}
}

void BarcodeToPositionMulti::initPackRepositoey()
{
	mRepo.packBuffer = new ReadPairPack * [PACK_NUM_LIMIT];
	memset(mRepo.packBuffer, 0, sizeof(ReadPairPack*) * PACK_NUM_LIMIT);
	mRepo.writePos = 0;
	mRepo.readPos = 0;
}

void BarcodeToPositionMulti::destroyPackRepository() {
	delete mRepo.packBuffer;
	mRepo.packBuffer = NULL;
}

void BarcodeToPositionMulti::producePack(ReadPairPack* pack) {
	mRepo.packBuffer[mRepo.writePos] = pack;
	mRepo.writePos++;
}


void BarcodeToPositionMulti::producerTask() {
	if (mOptions->verbose){
//		loginfo("start to load data");
    }
	long lastReported = 0;
	int slept = 0;
	long readNum = 0;
	ReadPair** data = new ReadPair * [PACK_SIZE];
	memset(data, 0, sizeof(ReadPair*) * PACK_SIZE);
	FastqReaderPair reader(mOptions->transBarcodeToPos.in1, mOptions->transBarcodeToPos.in2, true);
	int count = 0;
	bool needToBreak = false;
	while (true) {
		ReadPair* read = reader.read();
		if (!read || needToBreak) {
			ReadPairPack* pack = new ReadPairPack;
			pack->data = data;
			pack->count = count;
			producePack(pack);
			data = NULL;
			if (read) {
				delete read;
				read = NULL;
			}
			break;
		}
		data[count] = read;
		count++;
		if (mOptions->verbose && count + readNum >= lastReported + 1000000) {
			lastReported = count + readNum;
			string msg = "loaded " + to_string((lastReported / 1000000)) + "M read pairs";
//			loginfo(msg);
		}
		if (count == PACK_SIZE || needToBreak) {
			ReadPairPack* pack = new ReadPairPack;
			pack->data = data;
			pack->count = count;
			producePack(pack);
			//re-initialize data for next pack
			data = new ReadPair * [PACK_SIZE];
			memset(data, 0, sizeof(ReadPair*) * PACK_SIZE);
			// if the consumer is far behind this producer, sleep and wait to limit memory usage
			while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
				slept++;
				usleep(100);
			}
			readNum += count;
			// if the writer threads are far behind this producer, sleep and wait
			// check this only when necessary
			if (readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mWriter) {
				while ((mWriter && mWriter->bufferLength() > PACK_IN_MEM_LIMIT)) {
					slept++;
					usleep(1000);
				}
			}
			// reset count to 0
			count = 0;
		}
	}
	mProduceFinished = true;
	if (mOptions->verbose) {
//		loginfo("all reads loaded, start to monitor thread status");
	}
	if (data != NULL)
		delete[] data;
}


void BarcodeToPositionMulti::writeTask(WriterThread* config) {
	while (true) {
		//loginfo("writeTask running: " + config->getFilename());
		if (config->isCompleted()) {
			config->output();
			break;
		}
		config->output();
	}

	if (mOptions->verbose) {
		string msg = config->getFilename() + " writer finished";
//		loginfo(msg);
	}
}

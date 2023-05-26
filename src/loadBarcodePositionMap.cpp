

/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "loadBarcodePositionMap.h"

LoadBarcodePositionMap::LoadBarcodePositionMap(Options* opt)
{
	mOptions = opt;
	mProduceFinished = false;
	mFinishedThreads = 0;
	mZipFile = NULL;
	barcodeStart = opt->barcodeStart;
	barcodeLen = opt->barcodeLen;
	rc = opt->rc;
	if (opt->platform.compare(SEQ500) == 0) {
		isSeq500 = true;
	}
	else {
		isSeq500 = false;
	}
	sm = new SlideMask(opt->maskFile, isSeq500, opt->turnFovDegree);
	barcodePositionMap = new BarcodePositionMap(opt);
}

LoadBarcodePositionMap::~LoadBarcodePositionMap()
{
	if (sm) {
		delete sm;
	}
	if (barcodePositionMap) {
		delete barcodePositionMap;
	}
}

BarcodePositionMap* LoadBarcodePositionMap::load()
{
	initPackRepository();
	std::thread  producer(std::bind(&LoadBarcodePositionMap::producerTask, this));

	BarcodePositionConfig** configs = new BarcodePositionConfig* [mOptions->thread];
	for (int t = 0; t < mOptions->thread; t++) {
		configs[t] = new BarcodePositionConfig(mOptions, t);
	}

	std::thread** threads = new thread * [mOptions->thread];
	for (int t = 0; t < mOptions->thread; t++) {
		threads[t] = new std::thread(std::bind(&LoadBarcodePositionMap::consumerTask, this, configs[t]));
	}

	producer.join();
	for (int t = 0; t < mOptions->thread; t++) {
		threads[t]->join();
	}

	BarcodePositionConfig* finalConfig = BarcodePositionConfig::merge(configs);

	if (mOptions->verbose)
//		loginfo("start to generate reports\n");
	
	finalConfig->print();


	//clean up
	for (int t = 0; t < mOptions->thread; t++) {
		delete threads[t];
		threads[t] = NULL;
		delete configs[t];
		configs[t] = NULL;
	}

	barcodePositionMap->bpmap = finalConfig->bpmap;
	delete sm;
	return barcodePositionMap;
}

bool LoadBarcodePositionMap::processSingleEnd(ReadPack* pack, BarcodePositionConfig* config)
{
	pair<int, int> slidePos;
	for (int p = 0; p < pack->count; p++) {
		Read* read = pack->data[p];
		config->totalReads++;
		string barcode = read->mSeq.mStr.substr(barcodeStart, barcodeLen);
		if (barcode.find("N") != barcode.npos) {
			config->lowQReads++;
			delete read;
			continue;
		}
		uint64_t barcodeInt = seqEncode(barcode.c_str(), 0, barcodeLen, rc);

		read->getDNBidx(isSeq500);
		Position1 position ;
		slidePos = sm->getIdx(read->dnbIdx[0], read->dnbIdx[1], read->dnbIdx[2]);
		if (slidePos.first < 0 || slidePos.second < 0) {
			config->withoutPositionRead++;
			delete read;
			continue;
		}
		//position.fov_c = read->dnbIdx[1];
		//position.fov_r = read->dnbIdx[2];
		position.x = slidePos.first;
		position.y = slidePos.second;
		config->addBarcode(barcodeInt, position);
		delete read;
	}
	delete pack->data;
	delete pack;
	return true;
}

void LoadBarcodePositionMap::initPackRepository()
{
	mRepo.packBuffer = new ReadPack * [PACK_NUM_LIMIT];
	memset(mRepo.packBuffer, 0, sizeof(ReadPack*) * PACK_NUM_LIMIT);
	mRepo.writePos = 0;
	mRepo.readPos = 0;
}

void LoadBarcodePositionMap::destoryPackRepository()
{
	delete mRepo.packBuffer;
	mRepo.packBuffer = NULL;
}

void LoadBarcodePositionMap::producePack(ReadPack* pack)
{
	mRepo.packBuffer[mRepo.writePos] = pack;
	mRepo.writePos++;
}

void LoadBarcodePositionMap::consumePack(BarcodePositionConfig* config)
{
	ReadPack* data;
	mInputMtx.lock();
	while (mRepo.writePos <= mRepo.readPos) {
		usleep(1000);
		if (mProduceFinished) {
			mInputMtx.unlock();
			return;
		}
	}
	data = mRepo.packBuffer[mRepo.readPos];
	mRepo.readPos++;
	mInputMtx.unlock();
	processSingleEnd(data, config);
}

void LoadBarcodePositionMap::producerTask()
{
	if (mOptions->verbose){
//		loginfo("start to load single end data");
    }
	long lastReported = 0;
	int slept = 0;
	long readNum = 0;
	Read** data = new Read * [PACK_SIZE];
	memset(data, 0, sizeof(Read*) * PACK_SIZE);
	FastqReader reader(mOptions->transBarcodeToPos.in);
	int count = 0;
	bool needToBreak = false;
	while (true) {
		Read* read = reader.read();
		// TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
		if (!read || needToBreak) {
			//the last pack
			ReadPack* pack = new ReadPack;
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
			string msg = "loaded " + to_string(lastReported / 1000000) + "M reads";
//			loginfo(msg);
		}
		//a full pack
		if (count == PACK_SIZE || needToBreak) {
			ReadPack* pack = new ReadPack;
			pack->data = data;
			pack->count = count;
			producePack(pack);
			//re-initialize data for next pack
			data = new Read * [PACK_SIZE];
			memset(data, 0, sizeof(Read*) * PACK_SIZE);
			while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
				slept++;
				usleep(100);
			}
			readNum += count;
			count = 0;
		}
	}
	mProduceFinished = true;
	if (mOptions->verbose)
//		loginfo("all single end reads loaded, start to monitor thread status");
	if (data != NULL)
		delete[] data;
}

void LoadBarcodePositionMap::consumerTask(BarcodePositionConfig* config)
{
	while (true) {
		while (mRepo.writePos <= mRepo.readPos) {
			if (mProduceFinished)
				break;
			usleep(1000);
		}
		if (mProduceFinished && mRepo.writePos == mRepo.readPos) {
			mFinishedThreads++;
			if (mOptions->verbose) {
				string msg = "finished " + to_string(mFinishedThreads) + " threads. Data processing completed.";
//				loginfo(msg);
			}
			break;
		}
		if (mProduceFinished) {
			if (mOptions->verbose) {
				string msg = "thread is processing the " + to_string(mRepo.readPos) + "/" + to_string(mRepo.writePos) + " pack";
//				loginfo(msg);
			}
		}
		consumePack(config);	
	}
	if (mOptions->verbose) {
		string msg = "finished one thread";
//		loginfo(msg);
	}
}

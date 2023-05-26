/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#ifndef BARCODEPOSITIONMAPGET_H
#define BARCODEPOSITIONMAPGET_H

#include <string>
#include "barcodePositionMap.h"
#include "options.h"
#include "htmlreporter.h"

using namespace std;

class BarcodePositionMapGet {
public:
	BarcodePositionMapGet(Options* opt);
	~BarcodePositionMapGet();
	void process();
private:
	Options* mOptions;
};

#endif

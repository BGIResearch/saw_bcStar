/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "slideMask.h"

SlideMask::SlideMask(string maskFile, bool IsSeq500, int turnFovDegree, int FovGap)
{
	isSeq500 = IsSeq500;
	fovGap = FovGap;
	makeIdx(maskFile, turnFovDegree);
}
SlideMask::~SlideMask() {

}

pair<int32_t, int32_t> SlideMask::getIdx(int idx,  int fovCol, int fovRow)
{
	int32_t idxCol;
	int32_t idxRow;
	if (isSeq500 && idx > intOffset) {
		return pair<int32_t, int32_t>(-1, -1);
	}
	if (idx <= totaldnbs) {
		dnbPos = dnbIdx[idx];
		//cout << "getIdxPos: col " << dnbPos.first << " row " << dnbPos.second << endl;
	}
	else {
		cerr << "the dnb index: "<<  idx << " is out of the mask file, please check the data" << endl;
		return pair<int32_t, int32_t>(-1, -1);
	}
	idxCol = dnbPos.first + (fovCol - 1) * maxdnbCol + (fovCol - 1) * fovGap;
	idxRow = dnbPos.second + (fovRow - 1) * maxdnbRow + (fovRow - 1) * fovGap;	
	//cout << "dnb col: " << dnbPos.first << "\tdnb row: " << dnbPos.second << endl;
	return pair<int32_t, int32_t>(idxCol, idxRow);
}

bool SlideMask::isInIntermediate(int idx)
{
	if (idx <= intOffset)
		return true;
	else
		return false;
}

pair<int32_t, int32_t> SlideMask::getShift(int block_c, int block_r)
{
	int32_t rowShift = 0;
	int32_t colShift = 0;
	for (int row = 0; row < block_r; row++) {
		rowShift += maskOrder[row][block_c].second;
		rowShift += FOV_GAP;
	}
	for (int col = 0; col < block_c; col++) {
		colShift += maskOrder[block_r][col].first;
		colShift += FOV_GAP;
	}
	//cout << "block index: " << block_c << "\t" << block_r << "\tshift: " << colShift << "\t" << rowShift << endl;
	return pair<int32_t, int32_t>(colShift, rowShift);
}

void SlideMask::makeIdx(string maskFile, int turnFovDegree)
{

    //nc::NdArray<int> mask;
    //cout << "mask file: " << maskFile << endl;
    //mask = nc::fromfile<int>(maskFile, "\t");
    //mask.reshape(rows, cols);
    //nc::NdArray<int> line;
    //cout << line << endl;

    cout<<"mask file: "<<maskFile<<endl;
    int blockIndex;
    int blockCol;
    int blockRow;
    int dnbCol;
    int dnbRow;
    int rowShift;
    int colShift;
    string line;
    ifstream maskInStream;
    maskInStream.open(maskFile);
//    int i=0;
    getline(maskInStream,line);
//cout << line << endl;
    while(maskInStream>>blockIndex>>blockCol>>blockRow>>dnbCol>>dnbRow){
        totaldnbs+=dnbCol*dnbRow;
        //cout << "dnbCol: " << dnbCol << "  dnbRow: " << dnbRow << " total: " << dnbCol*dnbRow << endl;
        if(blockCol!=1&&blockCol!=10&&blockRow!=1&&blockRow!=10){
            intOffset+=dnbCol*dnbRow;
            //cout << "inner block: " << blockIndex << endl;
        }
        if(rowOrder.count(blockRow)<1){
            rowOrder[blockRow]=dnbRow;
            maxdnbRow+=dnbRow;
        }
        if(colOrder.count(blockCol)<1){
            colOrder[blockCol]=dnbCol;
            maxdnbCol+=dnbCol;
        }
    }
    dnbIdx=new pair<int, int>[totaldnbs];
    maskInStream.close();
    maxdnbRow+=(rowOrder.size()-1)*TRACK_WIDTH;
    maxdnbCol+=(colOrder.size()-1)*TRACK_WIDTH;
    maskInStream.open(maskFile);
    getline(maskInStream,line);
    int iter=0;
    while(maskInStream>>blockIndex>>blockCol>>blockRow>>dnbCol>>dnbRow){
        rowShift=0;
        colShift=0;
        iter++;
        for(int r=1;r<blockRow;r++){
            rowShift+=rowOrder.at(r);
            rowShift+=TRACK_WIDTH;
        }
        for(int c=1;c<blockCol;c++){
            colShift+=colOrder.at(c);
            colShift+=TRACK_WIDTH;
        }
        for(int dr=0;dr<dnbRow;dr++){
            for(int dc=0;dc<dnbCol;dc++){
                int xPos=0;
                int yPos=0;
                if(turnFovDegree==90){
                    xPos=rowShift+dr;
                    yPos=maxdnbCol-(colShift+dc)-1;
                }else
                    if(turnFovDegree==0){
                        xPos=colShift+dc;
                        yPos=rowShift+dr;
                    }else
                        if(turnFovDegree==180){
                            xPos=maxdnbCol-(colShift+dc)-1;
                            yPos=maxdnbRow-(rowShift+dr)-1;
                        }else{
                            string errMsg="turnFovDegree should be one of [0, 90, 180], but got: " + to_string(turnFovDegree);
						  sawErrCode(err_param_invalid_bcmap,errMsg);
                            cerr << "turnFovDegree should be one of [0, 90, 180], but got: " << turnFovDegree << endl;
//                            assert(false);
                            exit(1);
                        }
                pair<int,int> DNBPos(xPos,yPos);
                dnbIdx[idx]=DNBPos;
                idx++;

                //cout << "idx: " << idx << " x: " << xPos << " y: " << yPos << endl;
            }
        }
    }
    maskInStream.close();
    if(turnFovDegree==90){
        int temp=maxdnbRow;
        maxdnbRow=maxdnbCol;
        maxdnbCol=temp;
    }

    cout<<"makeIdx finished, dnbnumber: "<<totaldnbs<<"\t"<<"\tintOffset: "<<intOffset<<"\tturn fov degree: "<<turnFovDegree
            <<endl;
    cout<<"maxDNBCol: "<<maxdnbCol<<"  maxDNBRow: "<<maxdnbRow<<endl;
}

slideRange SlideMask::getSlideRange(int minCol, int maxCol, int minRow, int maxRow){
	slideRange sliderange;
	if (minCol == 0){
		minCol = 1;
	}
	if (minRow == 0){
		minRow = 1;
	}
	
	sliderange.colStart = (maxdnbCol + FOV_GAP)*(minCol-1);
	sliderange.colEnd = (maxdnbCol + FOV_GAP)*(maxCol)+1;
	sliderange.rowStart = (maxdnbRow + FOV_GAP)*(minRow-1);
	sliderange.rowEnd = (maxdnbRow + FOV_GAP)*(maxRow)+1;
	
	return sliderange;
}

uint16 SlideMask::getSlidePitch(std::string platform){
	return slidePitchMap.find(platform)->second;
}



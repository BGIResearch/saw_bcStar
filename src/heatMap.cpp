#include "heatMap.h"

HeatMap::HeatMap(){
    colStart = 0;
    colEnd = 0;
    rowStart = 0;
    rowEnd = 0;
}

HeatMap::HeatMap(long ColStart, long ColEnd, long RowStart, long RowEnd){
    init(ColStart, ColEnd, RowStart, RowEnd);
}

HeatMap::~HeatMap(){
    
}

void HeatMap::init(long ColStart, long ColEnd, long RowStart, long RowEnd){
    colStart = ColStart;
    colEnd = ColEnd;
    rowStart = RowStart;
    rowEnd = RowEnd;
    cout << "colStart: " << colStart << " colEnd: " << colEnd << " rowStart: " << rowStart << " rowEnd: " << rowEnd << endl;
    cout << "matrix: " << rowEnd-rowStart+1 << " : " << colEnd-colStart+1 << endl;
    heatMap = cv::Mat::zeros(int(rowEnd-rowStart+1), int(colEnd-colStart+1), CV_8UC1);
    
}

void HeatMap::setPointValue(long col, long row, int value){
    long shiftCol = col - colStart;
    long shiftRow = row - rowStart;
    //cout << "col: " << col << " row: " << row << " shiftCol: " << shiftCol << " shiftRow: " << shiftRow << endl;
    heatMap.at<uchar>(shiftRow, shiftCol) = value;
}

void HeatMap::addPointValue(long col, long row, int value){
    long shiftCol = col - colStart;
    long shiftRow = row - rowStart;
    //cout << "row: " << row << "\tcol: " << col << "\tshiftRow: " << shiftRow << "\tshiftCol: " << shiftCol << "\tvalue: " <<  value << endl;
    heatMap.at<uchar>(shiftRow, shiftCol) += value;
}

void HeatMap::saveHeatMap(string heatMapFile){
    cout << "#########save heatMap to tiff" <<endl;
    cv::imwrite(heatMapFile, heatMap);
    Mat kern = (Mat_<char>(5,5) << 1, 1 ,1, 1, 1,
                                   1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1);
    cv::Mat binImg;
    cv::filter2D(heatMap, binImg, heatMap.depth(), kern);
    cv::Mat shrinkImg;
    cv::resize(binImg, shrinkImg, cv::Size(), imgDownRate, imgDownRate);
    string shrinkImgFile = heatMapFile + ".shrink.png";
    cv::imwrite(shrinkImgFile, shrinkImg);
}
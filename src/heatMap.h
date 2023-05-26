#ifndef HEATMAP_H
#define HEATMAP_H

#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;

class HeatMap{
public:
    HeatMap();
    HeatMap(long ColStart, long ColEnd, long RowStart, long RowEnd);
    ~HeatMap();
    void setPointValue(long col, long row, int value);
    void addPointValue(long col, long row, int value);
    void saveHeatMap(string heaMapFile);
    void init(long ColStart, long ColEnd, long RowStart, long RowEnd);
public:
    long colStart;
    long colEnd;
    long rowStart;
    long rowEnd;
    cv::Mat heatMap;
private:
    float imgDownRate = 0.03;
};



#endif
#include "IncludeDefine.h"
#include "Genome.h"

bool sjAlignSplit(uint a1,uint aLength, const Genome &mapGen, uint &a1D, uint &aLengthD, uint &a1A, uint &aLengthA, uint &isj) {
    uint sj1=(a1-mapGen.sjGstart)%mapGen.sjdbLength;
    if (sj1<mapGen.sjdbOverhang && sj1+aLength>mapGen.sjdbOverhang) {//align crosses the junctions
      // isj 是 sj数组的索引位置
      // aLengthD 和 aLengthA 是跨过sj两边的比对长度。D是供体Donor端，A是受体Accept端
      // a1D 和 a1A 是两边的比对位置
        isj=(a1-mapGen.sjGstart)/mapGen.sjdbLength;
        aLengthD=mapGen.sjdbOverhang-sj1;
        aLengthA=aLength-aLengthD;
        a1D=mapGen.sjDstart[isj]+sj1;
        a1A=mapGen.sjAstart[isj];
        return true;
    } else {
        return false;
    };
};

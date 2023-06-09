#include "ReadAlign.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"

uint ReadAlign::maxMappableLength2strands(uint pieceStartIn, uint pieceLengthIn, uint iDir, uint iSA1, uint iSA2, uint& maxLbest, uint iFrag) {
    //returns number of mappings, maxMappedLength=mapped length
    uint Nrep=0, indStartEnd[2], maxL;

    uint NrepAll[P.pGe.gSAsparseD], indStartEndAll[P.pGe.gSAsparseD][2], maxLall[P.pGe.gSAsparseD];
    maxLbest=0;

    bool dirR = iDir==0;

    // defaults:  (from genomeParameters.txt)
    // gSAsparseD = 1
    // gSAindexNbases = 14

    for (uint iDist=0; iDist<min(pieceLengthIn,P.pGe.gSAsparseD); iDist++) {//cycle through different distances
        uint pieceStart;
        uint pieceLength=pieceLengthIn-iDist;

        //calculate full index
        uint Lmax=min(P.pGe.gSAindexNbases,pieceLength);
        uint ind1=0;
        //将序列转换成0123编码的数值，存储到ind1
        if (dirR) {//forward search
            pieceStart=pieceStartIn+iDist;
            for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
                ind1 <<=2LLU;
                ind1 += ((uint) Read1[0][pieceStart+ii]);
            };
        } else {//reverse search
            pieceStart=pieceStartIn-iDist;
            for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
                ind1 <<=2LLU;
                ind1 += ( 3-((uint) Read1[0][pieceStart-ii]) );
            };
        };

        //find SA boundaries
        // Lind等于14或者剩下的序列长度
        uint Lind=Lmax;
        while (Lind>0) {//check the precense of the prefix for Lind
            iSA1=mapGen.SAi[mapGen.genomeSAindexStart[Lind-1]+ind1]; // starting point for suffix array search.
            if ((iSA1 & mapGen.SAiMarkAbsentMaskC) == 0) {//prefix exists
                break;
            } else {//this prefix does not exist, reduce Lind
                --Lind;
                ind1 = ind1 >> 2;
            };
        };

        // define lower bound for suffix array range search.
        if (mapGen.genomeSAindexStart[Lind-1]+ind1+1 < mapGen.genomeSAindexStart[Lind]) {//we are not at the end of the SA
		  if(mapGen.SAi[mapGen.genomeSAindexStart[Lind-1]+ind1+1]==0){
			iSA2=mapGen.nSA-1;
		  }else {
			iSA2 = ((mapGen.SAi[mapGen.genomeSAindexStart[Lind - 1] + ind1 + 1]&mapGen.SAiMarkNmask)&mapGen.SAiMarkAbsentMask) - 1;
		  }
        } else {
            iSA2=mapGen.nSA-1;
        };
		// iSA1为上界，iSA2为下界

    //#define SA_SEARCH_FULL

    #ifdef SA_SEARCH_FULL
        //full search of the array even if the index search gave maxL
        maxL=0;
        Nrep = maxMappableLength(mapGen, Read1, pieceStart, pieceLength, iSA1 & mapGen.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
    #else
        if (Lind < P.pGe.gSAindexNbases && (iSA1 & mapGen.SAiMarkNmaskC)==0 ) {//no need for SA search
            // very short seq, already found hits in suffix array w/o having to search the genome for extensions.
            // 可以看出indStartEnd存放的是比对位置的上下界，以及这一步骤需要产生如下几个数据给下游用：indStartEnd、Nrep和maxL
            indStartEnd[0]=iSA1;
            indStartEnd[1]=iSA2;
            //如果完全匹配上，则无需跟genome进行延伸对比，查询到的位置数量为上下界间的数量，即Nrep。maxL为序列长度
            Nrep=indStartEnd[1]-indStartEnd[0]+1;
            maxL=Lind;
        } else if (iSA1==iSA2) {//unique align already, just find maxL
            if ((iSA1 & mapGen.SAiMarkNmaskC)!=0) {
                ostringstream errOut;
                errOut  << "BUG: in ReadAlign::maxMappableLength2strands";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
            };
            indStartEnd[0]=indStartEnd[1]=iSA1;
            Nrep=1;
            bool comparRes;
            // 为unique比对时，对序列紧紧需要做一次延伸，跟genome序列对比找出MMP
            maxL=compareSeqToGenome(mapGen, Read1, pieceStart, pieceLength, Lind, iSA1, dirR, comparRes);
        } else {//SA search, pieceLength>maxL
            if ( (iSA1 & mapGen.SAiMarkNmaskC)==0 ) {//no N in the prefix
                maxL=Lind;
            } else {
                maxL=0;
            };
            // 大部分时间会走到这个分支，需要通过2分查找去在上下界间寻找最长的匹配序列，即mmp
            Nrep = maxMappableLength(mapGen, Read1, pieceStart, pieceLength, iSA1 & mapGen.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
        };
    #endif
		// iDist=0
        if (maxL+iDist > maxLbest) {//this idist is better
            maxLbest=maxL+iDist;
        };
        //maxLbest 就是mmp
        NrepAll[iDist]=Nrep;
        indStartEndAll[iDist][0]=indStartEnd[0];
        indStartEndAll[iDist][1]=indStartEnd[1];
        maxLall[iDist]=maxL;
    };
	// sparsity 表示稀疏性，1表示是连续的，如果为2呢，间隔1个？应该是代码开发时的测试参数
    for (uint iDist=0; iDist<min(pieceLengthIn,P.pGe.gSAsparseD); iDist++) {//cycle through different distances, store the ones with largest maxL
        if ( (maxLall[iDist]+iDist) == maxLbest) {
          //将mmp相关信息存起来，第二个参数表示在read中的起始位置，即偏移量，后续做stitch会用到
            storeAligns(iDir, (dirR ? pieceStartIn+iDist : pieceStartIn-iDist), NrepAll[iDist], maxLall[iDist], indStartEndAll[iDist], iFrag);
        };
    };
    return Nrep;
};

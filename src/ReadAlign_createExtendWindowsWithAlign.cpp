#include "IncludeDefine.h"
#include "Parameters.h"
//#include "Transcript.h"
#include "ReadAlign.h"
#include "SequenceFuns.h"

int ReadAlign::createExtendWindowsWithAlign(uint a1, uint aStr) {


    uint aBin = (a1 >> P.winBinNbits); //align's bin
    uint iBinLeft=aBin, iBinRight=aBin;
    uintWinBin* wB=winBin[aStr];
    uint iBin=-1, iWin=-1, iWinRight=-1;

#ifdef DEBUGOUT
//  cout<<"createExtWinWithAlign\t"<<aBin<<"\t";
#endif
    if (wB[aBin]==uintWinBinMax) {//proceed if there is no window at this bin
        //check neighboring bins
//		cout<<"no window\t";
        bool flagMergeLeft=false;
        if (aBin>0) {//merge left only if there are bins on the left
          // 往左看9个bin
          //winAnchorDistNbins默认为9，猜测是因为4个iBin大小为一个window，往左合并，即跨越2个window，8个iBin
            for (iBin=aBin-1;  iBin >= ( aBin>P.winAnchorDistNbins ? aBin-P.winAnchorDistNbins : 0 );  --iBin) {//go left, find windows in Anchor range
                if (wB[iBin]<uintWinBinMax) {
                    flagMergeLeft=true;
                    break;
                };
                if (iBin==0) break;
            };
            // P.winBinChrNbits =2 from debug
            // chrBin存储每个bin所属哪条染色体,每个bin大小为2^18（即genomeChrBinNbases），推测这里面有个公式 genomeChrBinNbases=winBinNbits+winBinChrNbits
            flagMergeLeft = flagMergeLeft && (mapGen.chrBin[iBin>>P.winBinChrNbits]==mapGen.chrBin[aBin>>P.winBinChrNbits]);
            if (flagMergeLeft) {//this align can be merged into the existing window
                iWin=wB[iBin];
                iBinLeft=WC[iWin][WC_gStart];
                for (uint ii=iBin+1; ii<=aBin; ii++) {//mark al bins with the existing windows ID
                    wB[ii]=iWin;
                };
            };
        };
        if(flagMergeLeft){
//          cout<<"merge left success\t"<<iBin<<"\t"<<iWin<<"\t"<<iBinLeft<<"\t";
        }
        bool flagMergeRight=false;
        if (aBin+1<P.winBinN) {//merge left only if there are bins on the right
          // 往右看9个bin
            for (iBin=aBin+1;  iBin<min(aBin+P.winAnchorDistNbins+1,P.winBinN);  ++iBin) {//go right, find windows in Anchor range
                if (wB[iBin]<uintWinBinMax) {
                    flagMergeRight=true;
                    break;
                };
            };

            flagMergeRight = flagMergeRight && (mapGen.chrBin[iBin>>P.winBinChrNbits]==mapGen.chrBin[aBin>>P.winBinChrNbits]);
            if (flagMergeRight) {//this align can be merged into the existing window
//			  cout<<"merge right success\t"<<iBin<<"\t";
                while (wB[iBin]==wB[iBin+1]) ++iBin; //extend through all bins of the right window
//                cout<<iBin<<"\t";
                iBinRight=iBin;
                iWinRight=wB[iBin];
//			  cout<<iBinRight<<"\t"<<iWinRight<<"\t";
//                 if (iBin!=WC[iWinRight][WC_gEnd]) {//debug, switch off!!!
//                     cerr <<"BUG in createWindows"<<endl<<flush;
//                     exit(0);
//                 };
                if (!flagMergeLeft) iWin=wB[iBin];//left window overwrites right
                for (uint ii=aBin; ii<=iBin; ii++) {//mark al bins with the existing windows ID
                    wB[ii]=iWin;
                };
            };
        };

        if (!flagMergeLeft && !flagMergeRight) {//no merging, a new window was added
//		  cout<<"no merge\t";
            wB[aBin]=iWin=nW; //add new window ID for now, may change it later
            WC[iWin][WC_Chr]=mapGen.chrBin[aBin >> P.winBinChrNbits];
            WC[iWin][WC_Str]=aStr;
            WC[iWin][WC_gEnd]=WC[iWin][WC_gStart]=aBin;
            ++nW;
//            cout<<aBin<<"\t"<<wB[aBin]<<"\t"<< WC[iWin][WC_Chr]<<"\t"<<WC[iWin][WC_Str]<<"\t"<<WC[iWin][WC_gStart]<<"\t"<<WC[iWin][WC_gEnd];
            if (nW>=P.alignWindowsPerReadNmax) {
                nW=P.alignWindowsPerReadNmax-1;
                return EXIT_createExtendWindowsWithAlign_TOO_MANY_WINDOWS; //too many windows, do not record TODO: record a marker
            };
        } else {//record windows after merging
            WC[iWin][WC_gStart]=iBinLeft;
            WC[iWin][WC_gEnd]=iBinRight;
            if (flagMergeLeft && flagMergeRight) {//kill right window, it was merged with the left one
                WC[iWinRight][WC_gStart]=1;
                WC[iWinRight][WC_gEnd]=0;
            };
//            cout<<"merge done\t"<<WC[iWin][WC_gStart]<<"\t"<<WC[iWin][WC_gEnd];
        };
//        cout<<endl;
    }else{
//      cout<<"has win\t"<<wB[aBin]<<endl;
    }
    return 0;
};


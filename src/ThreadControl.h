#ifndef THREAD_CONTROL_DEF
#define THREAD_CONTROL_DEF

#include "ReadAlignChunk.h"
#include <pthread.h>
#include "JoinedMultiThreadPara.h"

#define MAX_chunkOutBAMposition 100000

class ThreadControl {
public:
    bool threadBool;

    pthread_t *threadArray;
    pthread_mutex_t mutexInRead, mutexOutSAM, mutexOutBAM1, mutexOutChimSAM, mutexOutChimJunction, mutexOutUnmappedFastx, mutexOutFilterBySJout;
    pthread_mutex_t mutexStats, mutexLogMain, mutexBAMsortBins, mutexError;

    uint chunkInN,chunkOutN;

    ThreadControl();

//    static void* threadRAprocessChunks(void *RAchunk) {
//        ( (ReadAlignChunk*) RAchunk )->processChunks();
//        pthread_exit(0);
//        return NULL;
//    };
    static void* threadRAprocessChunks(void *JoinedPara) {
        JoinedMultiThreadPara* para=(JoinedMultiThreadPara*)JoinedPara;

        // do star
        para->starPara->processChunks(para->bmPara);
        pthread_exit(0);
        return NULL;
    };
};

#endif


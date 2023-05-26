all:BCSTAR
# or these user-set flags that will be added to standard flags
LDFLAGSextra ?=
CXXFLAGSextra ?=
src := ./src
obj := ./obj


# pre-defined flags
LDFLAGS_shared := -pthread -lz -lhts  


CXXFLAGS_common := -pipe -std=c++14 -Wall -Wextra -fopenmp -Isrc/folly/include
CXXFLAGS_main := -O3 $(CXXFLAGS_common) 

CFLAGS := -O3 -g -pipe -Wall -Wextra $(CFLAGS) 


##########################################################################################################

OBJECTS = ${obj}/ParametersChimeric_initialize.o ${obj}/ParametersSolo.o ${obj}/SoloRead.o ${obj}/SoloRead_record.o ${obj}/SoloReadBarcode.o ${obj}/SoloReadBarcode_getCBandUMI.o ${obj}/SoloReadFeature.o ${obj}/SoloReadFeature_record.o ${obj}/SoloReadFeature_inputRecords.o ${obj}/Solo.o ${obj}/SoloFeature.o ${obj}/SoloFeature_collapseUMI.o ${obj}/SoloFeature_outputResults.o ${obj}/SoloFeature_processRecords.o ${obj}/ReadAlign_outputAlignments.o ${obj}/ReadAlign.o ${obj}/STAR.o ${obj}/SharedMemory.o ${obj}/PackedArray.o ${obj}/SuffixArrayFuns.o ${obj}/Parameters.o ${obj}/InOutStreams.o ${obj}/SequenceFuns.o ${obj}/Genome.o ${obj}/Stats.o ${obj}/Transcript.o ${obj}/Transcript_alignScore.o ${obj}/Transcript_generateCigarP.o ${obj}/Chain.o ${obj}/Transcript_variationAdjust.o ${obj}/Variation.o ${obj}/ReadAlign_waspMap.o ${obj}/ReadAlign_storeAligns.o ${obj}/ReadAlign_stitchPieces.o ${obj}/ReadAlign_multMapSelect.o ${obj}/ReadAlign_mapOneRead.o ${obj}/readLoad.o ${obj}/ReadAlignChunk.o ${obj}/ReadAlignChunk_processChunks.o ${obj}/ReadAlignChunk_mapChunk.o ${obj}/OutSJ.o ${obj}/outputSJ.o ${obj}/blocksOverlap.o ${obj}/ThreadControl.o ${obj}/sysRemoveDir.o ${obj}/ReadAlign_maxMappableLength2strands.o ${obj}/binarySearch2.o ${obj}/ReadAlign_outputTranscriptSAM.o ${obj}/ReadAlign_outputTranscriptSJ.o ${obj}/ReadAlign_outputTranscriptCIGARp.o ${obj}/ReadAlign_createExtendWindowsWithAlign.o ${obj}/ReadAlign_assignAlignToWindow.o ${obj}/ReadAlign_oneRead.o ${obj}/ReadAlign_stitchWindowSeeds.o ${obj}/ReadAlign_peOverlapMergeMap.o ${obj}/ReadAlign_mappedFilter.o ${obj}/ReadAlign_chimericDetection.o ${obj}/ReadAlign_chimericDetectionOld.o ${obj}/ReadAlign_chimericDetectionOldOutput.o ${obj}/ChimericDetection.o ${obj}/ChimericDetection_chimericDetectionMult.o ${obj}/ReadAlign_chimericDetectionPEmerged.o ${obj}/ChimericAlign_chimericJunctionOutput.o ${obj}/ChimericAlign_chimericStitching.o ${obj}/stitchWindowAligns.o ${obj}/extendAlign.o ${obj}/stitchAlignToTranscript.o ${obj}/alignSmithWaterman.o ${obj}/Genome_genomeGenerate.o ${obj}/genomeParametersWrite.o ${obj}/genomeScanFastaFiles.o ${obj}/genomeSAindex.o ${obj}/Genome_insertSequences.o ${obj}/insertSeqSA.o ${obj}/funCompareUintAndSuffixes.o ${obj}/funCompareUintAndSuffixesMemcmp.o ${obj}/TimeFunctions.o ${obj}/ErrorWarning.o ${obj}/loadGTF.o ${obj}/streamFuns.o ${obj}/stringSubstituteAll.o ${obj}/Transcriptome.o ${obj}/Transcriptome_quantAlign.o ${obj}/Transcriptome_geneFullAlignOverlap.o ${obj}/ReadAlign_quantTranscriptome.o ${obj}/Quantifications.o ${obj}/Transcriptome_geneCountsAddAlign.o ${obj}/sjdbLoadFromFiles.o ${obj}/sjdbLoadFromStream.o ${obj}/sjdbPrepare.o ${obj}/sjdbBuildIndex.o ${obj}/sjdbInsertJunctions.o ${obj}/mapThreadsSpawn.o ${obj}/Parameters_readSAMheader.o ${obj}/bam_cat.o ${obj}/serviceFuns.o ${obj}/BAMoutput.o ${obj}/BAMfunctions.o ${obj}/ReadAlign_alignBAM.o ${obj}/BAMbinSortByCoordinate.o ${obj}/signalFromBAM.o ${obj}/bamRemoveDuplicates.o ${obj}/BAMbinSortUnmapped.o ${obj}/barcodeMap.o ${obj}/barcodeOverlap.o ${obj}/barcodePositionConfig.o ${obj}/barcodePositionMap.o ${obj}/barcodePositionMapGet.o ${obj}/barcodeProcessor.o ${obj}/barcodeStat.o ${obj}/barcodeToPositionMulti.o ${obj}/chipMaskHDF5.o ${obj}/fastqreader.o ${obj}/htmlreporter.o ${obj}/loadBarcodePositionMap.o ${obj}/options.o ${obj}/read.o ${obj}/result.o ${obj}/sequence.o ${obj}/slideMask.o ${obj}/writer.o ${obj}/writerThread.o ${obj}/bloomFilter.o ${obj}/mgzipParser.o ${obj}/GlobalVariables.o ${obj}/Parameters_closeReadsFiles.o ${obj}/Parameters_openReadsFiles.o ${obj}/ReadAlign_calcCIGAR.o ${obj}/ChimericSegment.o ${obj}/ChimericAlign.o

${obj}/%.o:${src}/%.cpp mk_dir
	$(CXX) $(CXXFLAGS) -c $< -o $@
${obj}/%.o:${src}/%.c
	$(CXX) $(CXXFLAGS) -c $< -o $@
mk_dir:
	@if test ! -d $(obj);\
	then\
		mkdir $(obj);\
	fi



.PHONY: clean
clean:
	rm -f ${obj}/*.o bcStar  Depend.list


.PHONY: cleanRelease
cleanRelease:
	rm -f *.o Depend.list
	$(MAKE) -C htslib clean

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleanRelease)
Depend.list: $(SOURCES) parametersDefault.xxd
	echo $(SOURCES)
	'rm' -f ./Depend.list
	$(CXX) $(CXXFLAGS_common) -MM $^ >> Depend.list
include Depend.list
endif
endif

htslib : htslib/libhts.a

htslib/libhts.a :
	$(MAKE) -C htslib lib-static

parametersDefault.xxd: parametersDefault
	xxd -i parametersDefault > parametersDefault.xxd

BCSTAR : CXXFLAGS := $(CXXFLAGSextra) $(CXXFLAGS_main) $(CXXFLAGS)
BCSTAR : LDFLAGS := $(LDFLAGSextra) $(LDFLAGS_shared) $(LDFLAGS)
BCSTAR : Depend.list parametersDefault.xxd $(OBJECTS)
	$(CXX) -o bcStar $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -lhdf5 -ldl src/folly/lib/libfolly.a src/folly/lib/libfmt.a -ldeflate src/mgzip/libmgzf.a src/zstd/libzstd.a
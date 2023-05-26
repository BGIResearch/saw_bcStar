/* Copyright (C) BGI-Reasearch - All Rights Reserved
* Unauthorized copying of this file, via any medium is strictly prohibited
* Proprietary and confidential
* Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
*/
#include "chipMaskHDF5.h"
#include <unistd.h>

ChipMaskHDF5::ChipMaskHDF5(std::string FileName){
    fileName = FileName;
}

ChipMaskHDF5::~ChipMaskHDF5(){

}

void ChipMaskHDF5::creatFile(){
    fileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //cout << "create hdf5 file successfully." << endl;
}

void ChipMaskHDF5::openFile(){
    fileID = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
}

herr_t ChipMaskHDF5::writeDataSet(std::string chipID, slideRange& sliderange, unordered_map<uint64_t, Position1>& bpMap, int barcodeLen, int segment, int slidePitch, uint compressionLevel, int index){
    //generate dataSet space
    hid_t dataspaceID;
    hsize_t dims[RANK];
    dims[0] = sliderange.rowEnd - sliderange.rowStart + 1;
    dims[1] = sliderange.colEnd - sliderange.colStart + 1;
    dims[2] = segment;
    dataspaceID = H5Screate_simple(RANK, dims, NULL);
    hsize_t memDims[1] = {dims[0]*dims[1]*dims[2]};
    hid_t memdataspaceID = H5Screate_simple(1, memDims, NULL);
    //transfer bpMap to bpMatrix
    bpMatrix = new uint64_t**[dims[0]];
    uint64_t* bpMatrix_buffer = new uint64_t[dims[0]*dims[1]*dims[2]]();
    for (uint32_t i = 0 ; i < dims[0] ; i++){
        bpMatrix[i] = new uint64_t*[dims[1]];
        for (uint32_t j = 0; j < dims[1]; j++){
            bpMatrix[i][j] = bpMatrix_buffer + i*dims[1]*dims[2] + j*dims[2];
        }
    }
    uint64_t* barcodes;
    for(auto mapIter = bpMap.begin(); mapIter!=bpMap.end(); mapIter++){
        int row = mapIter->second.y-sliderange.rowStart;
        int col = mapIter->second.x-sliderange.colStart;
        barcodes = bpMatrix[row][col];
        for(uint32_t i=0; i<(uint32_t)segment; i++){
            if (barcodes[i] == 0){
                barcodes[i] = mapIter->first;
                //cout << mapIter->second.y-sliderange.rowStart << " : " << mapIter->second.x-sliderange.colStart << " : " << mapIter->first << endl;
                break;
            }
        }
        //cout << mapIter->second.y-sliderange.rowStart << " : " << mapIter->second.x-sliderange.colStart << " : " << mapIter->first << endl;
    }

    //define dataSet chunk, compression need chunking data
    hid_t plistID;
    //create property list of dataSet generation
    plistID = H5Pcreate(H5P_DATASET_CREATE);

    hsize_t cdims[RANK] = {CDIM0, CDIM1, static_cast<hsize_t>(segment)};
    herr_t status;
    status=H5Pset_chunk(plistID, RANK, cdims);
    status=H5Pset_deflate(plistID, compressionLevel);

    //create dataset
    hid_t datasetID;
    std::string datasetName = DATASETNAME + std::to_string(index);
    datasetID = H5Dcreate2(fileID, datasetName.c_str(), H5T_NATIVE_UINT64, dataspaceID, H5P_DEFAULT, plistID, H5P_DEFAULT);

    //create dataset attribute [rowShift, colShift, barcodeLength, slidePitch]
    uint32_t attributeValues[ATTRIBUTEDIM] = {sliderange.rowStart, sliderange.colStart, static_cast<uint32_t>(barcodeLen), static_cast<uint32_t>(slidePitch)};
    //create attribute to store dnb infomation
    hsize_t dimsA[1] = {ATTRIBUTEDIM};
    hid_t attributeSpace = H5Screate_simple(1, dimsA, NULL);
    hid_t attributeID = H5Acreate(datasetID, ATTRIBUTENAME, H5T_NATIVE_UINT32, attributeSpace, H5P_DEFAULT, H5P_DEFAULT);

    //create attribute to sotre chip id
    hid_t attributeSpace1 = H5Screate(H5S_SCALAR);
    hid_t attributeID1 = H5Acreate(datasetID, ATTRIBUTENAME1, H5T_NATIVE_CHAR, attributeSpace1, H5P_DEFAULT, H5P_DEFAULT);

    //write attribute
    H5Awrite(attributeID, H5T_NATIVE_INT, &attributeValues[0]);
    H5Awrite(attributeID1, H5T_NATIVE_CHAR, chipID.c_str());

    //close attribute
    H5Sclose(attributeSpace);
    H5Sclose(attributeSpace1);
    H5Aclose(attributeID);
    H5Aclose(attributeID1);

    //write matrix data to hdf5 file dataset
    H5Dwrite(datasetID, H5T_NATIVE_UINT64, memdataspaceID, dataspaceID, H5P_DEFAULT, &bpMatrix_buffer[0]);
    H5Sclose(dataspaceID);
    H5Dclose(datasetID);
    H5Pclose(plistID);
    H5Fclose(fileID);

    delete[] bpMatrix;
    delete[] bpMatrix_buffer;
    return status;
}



void ChipMaskHDF5::readDataSet(folly::F14ValueMap<uint64_t, Position1>& bpMap,int index){
//  herr_t status;
  //open dataset with datasetName
  std::string datasetName = DATASETNAME + std::to_string(index);
  hid_t datasetID = H5Dopen2(fileID, datasetName.c_str(), H5P_DEFAULT);
  
  //read attribute of the dataset
  
  //uint32_t attributeValues[ATTRIBUTEDIM];
  //hid_t attributeID = H5Aopen_by_name(fileID, datasetName.c_str(), ATTRIBUTENAME, H5P_DEFAULT, H5P_DEFAULT);
  //H5Aread(attributeID, H5T_NATIVE_UINT32, &attributeValues[0]);
  //cout << "attribute values: " << attributeValues[0] << " "<< attributeValues[1] << endl;
  hid_t dspaceID = H5Dget_space(datasetID);
//  hid_t dtype_id = H5Dget_type(datasetID);
//  hid_t plistID = H5Dget_create_plist(datasetID);
  int rank = H5Sget_simple_extent_ndims(dspaceID);
  hsize_t dims[rank];
  H5Sget_simple_extent_dims(dspaceID, dims, NULL);
  
  uint64_t matrixLen = 1;
  for (int i = 0 ; i<rank; i++){
	matrixLen *= dims[i];
  }
  
  int segment = 1;
  if (rank>=3){
	segment = dims[2];
  }
//  pid_t pd=getpid();
//  cout<<"current memory use "<<get_proc_virtualmem(pd)<<"\t"<<get_proc_virtualmemPeak(pd)<<endl;
//cout<<"use hyperslab"<<endl;
  uint64_t cpRow=dims[0]/100==0?1:dims[0]/100;
  uint64_t* bpMatrix_buffer2 = new uint64_t[cpRow*dims[1]*dims[2]];
	for(uint64_t rowIdx=0;rowIdx<=dims[0];rowIdx+=cpRow){
	  uint64_t realCpRows=rowIdx+cpRow>dims[0]?dims[0]-rowIdx:cpRow;
	  uint64_t cpSize=realCpRows*dims[1]*dims[2];
	  
	  hsize_t dimsm2[1];
	  dimsm2[0]=cpSize;
	  hid_t memspace2=H5Screate_simple(1,dimsm2,NULL);
	  hsize_t offset[3],count[3];
	  offset[0] = rowIdx;
	  offset[1] = 0;
	  offset[2] = 0;
	  count[0]  = realCpRows;
	  count[1]  = dims[1];
	  count[2]  = 1;
	  hsize_t offset_out2[1];
	  hsize_t count_out2[1];
	  offset_out2[0]=0;
	  count_out2[0]=cpSize;
	  if(count[0]*count[1]!=cpSize){
		cerr<<"different size\t"<<rowIdx<<"\t"<<count[0]<<"\t"<<count[1]<<"\t"<<cpSize<<endl;
		exit(1);
	  }
	  H5Sselect_hyperslab(memspace2, H5S_SELECT_SET, offset_out2, NULL,
								   count_out2, NULL);
	  H5Sselect_hyperslab(dspaceID, H5S_SELECT_SET, offset, NULL,
								   count, NULL);
	  H5Dread(datasetID, H5T_NATIVE_ULONG, memspace2, dspaceID,
					   H5P_DEFAULT, bpMatrix_buffer2);
	  for (uint32_t r = 0; r < realCpRows; r++){
		//bpMatrix[r] = new uint64_t*[dims[1]];
		for (uint32_t c = 0; c< dims[1]; c++){
		  //bpMatrix[r][c] = bpMatrix_buffer + r*dims[1]*dims[2] + c*dims[2];
		  Position1 position;
		  position.x=c;
		  position.y=r+rowIdx;
		  position.count=0;
		  position.fullTag=0;
		  if (rank >= 3 ){
			segment = dims[2];
			for (int s = 0; s<segment; s++){
			  uint64_t barcodeInt = bpMatrix_buffer2[r*dims[1]*segment + c*segment + s];
			  if (barcodeInt == 0){
				continue;
			  }
			  bpMap[barcodeInt] = position;
			}
		  }else{
			uint64_t barcodeInt = bpMatrix_buffer2[r*dims[1]+c];
			if (barcodeInt == 0){
			  continue;
			}
			bpMap[barcodeInt] = position;
		  }
		}
	  }
	}
  delete[] bpMatrix_buffer2;
}

int ChipMaskHDF5::readDataSetParallelize(folly::F14ValueMap<uint64_t, Position1>& bpMap, BloomFilter *&bloomFilter,int index){
  herr_t status;
  //open dataset with datasetName
  std::string datasetName = DATASETNAME + std::to_string(index);
  hid_t datasetID = H5Dopen2(fileID, datasetName.c_str(), H5P_DEFAULT);
//#ifdef OPT_MEMORY_PLAN_B
//  cout<<"before reserve:\t"<<bpMap.getAllocatedMemorySize()<<endl;
//  uint64_t bcNum=5073598806;
//  bpMap.reserve(bcNum);
//  cout<<"after reserve:\t"<<bpMap.getAllocatedMemorySize()<<endl;
//#endif
  hid_t dspaceID = H5Dget_space(datasetID);
  bloomFilter=nullptr;
//  hid_t dtype_id = H5Dget_type(datasetID);
  
  H5Dget_create_plist(datasetID);
  
  int rank = H5Sget_simple_extent_ndims(dspaceID);
  
  hsize_t dims[rank];
  H5Sget_simple_extent_dims(dspaceID, dims, NULL);
  
  
  uint64 matrixLen = 1;
  for (int i = 0; i < rank; i++) {
	matrixLen *= dims[i];
  }
  

  H5Sselect_all(dspaceID);
  hsize_t nchunks;
  status=H5Dget_num_chunks(datasetID, dspaceID, &nchunks);
#ifdef D_PRINT
  cout<<"total chunks number:\t"<<nchunks<<endl;
#endif
  if(status<0){
    return status;
  }
  size_t chunk_dims[rank];
  int flag = 0;
  for (int r = 0; r < rank; r++) {
	if (dims[r] == 1) {
	  flag++;
	  chunk_dims[r] = 1;
	} else {
	  chunk_dims[r] = 0;
	}
	
  }
  hsize_t offset_next[rank];
  for (uint64_t cid = 1; cid < nchunks; cid++) {
	H5Dget_chunk_info(datasetID, dspaceID, cid, offset_next, NULL, NULL, NULL);
	for (int r = 0; r < rank; r++) {
	  if (!chunk_dims[r] && offset_next[r]) {
		chunk_dims[r] = offset_next[r];
		flag++;
	  }
	}
	if (flag >= rank) {
	  break;
	}
  }
  for (int r = 0; r < rank; r++){
	if(chunk_dims[r] == 0){
	  chunk_dims[r] = dims[r];
	}
  }
  hsize_t chunk_len = 1;
  for (int r = 0; r < rank; r++) {
	chunk_len *= chunk_dims[r];
  }
  
  uint64_t chunkBufSize=10;
  uint64 **compressed_buffer = new uint64 *[chunkBufSize];
  uint64 **buffer = new uint64 *[chunkBufSize];
  hsize_t **offset = new hsize_t *[chunkBufSize];
  uint64 **compressed_buffer2 = new uint64 *[chunkBufSize];
  uint64 **buffer2 = new uint64 *[chunkBufSize];
  hsize_t **offset2 = new hsize_t *[chunkBufSize];
  for (uint64_t i = 0; i < chunkBufSize; i++) {
	compressed_buffer[i] = new uint64[chunk_len];
	buffer[i] = new uint64[chunk_len];
	offset[i] = new hsize_t[rank];
	compressed_buffer2[i] = new uint64[chunk_len];
	buffer2[i] = new uint64[chunk_len];
	offset2[i] = new hsize_t[rank];
  }
  hsize_t *chunk_size = new hsize_t[nchunks];
  
  
//  uint32 bpmap_num = 0;
//  bloomFilter = new BloomFilter();
  
  uint64_t chunk_num = 0;
  uint64_t now_chunk[3];
#ifdef OPT_MEMORY_IN_H5
  string buf1status="empty";
  string buf2status="empty";
  uint64_t buf1size=0;
  uint64_t buf2size=0;
#endif
  int openMPthreads=3;
#ifdef OPT_MEMORY_IN_H5
  openMPthreads=2;
#endif
#pragma omp parallel num_threads(openMPthreads)
  {
	int num_id = omp_get_thread_num();
	if (num_id == 0) {
	  libdeflate_decompressor *decompressor = libdeflate_alloc_decompressor();
	  now_chunk[0]=0;
	  uint64_t chunk_start=0;
	  uint64_t chunk_end=chunk_start+chunkBufSize>nchunks?nchunks:chunk_start+chunkBufSize;
//	  cout<<"0\t"<<now_chunk[0]<<"\t"<<nchunks<<endl;
	  while(now_chunk[0]<nchunks){
	    if(buf1status=="empty"){
		  for(uint64_t chunk_index=chunk_start;chunk_index<chunk_end;chunk_index++){
		    uint64_t patchIdx=chunk_index-chunk_start;
			uint32_t filter = 0;
			size_t actual_out = 0;
			H5Dget_chunk_info(datasetID, dspaceID, chunk_index, offset[patchIdx], &filter, NULL,
							  &chunk_size[patchIdx]);
			H5Dread_chunk(datasetID, H5P_DEFAULT, offset[patchIdx], &filter,
								   compressed_buffer[patchIdx]);
			libdeflate_zlib_decompress(decompressor, (void *) compressed_buffer[patchIdx],
													 chunk_size[patchIdx], (void *) buffer[patchIdx],
													 chunk_len * sizeof(uint64), &actual_out);
			chunk_num++;
			now_chunk[0]++;
		  }
		  buf1size=chunk_end-chunk_start;
		  assert(buf1size<=chunkBufSize);
		  buf1status="full";
		  chunk_start=chunk_end;
		  chunk_end=chunk_start+chunkBufSize>nchunks?nchunks:chunk_start+chunkBufSize;
		}
		if(buf2status=="empty"){
		  if(buf1status=="full"){
		    uint64 **tmpPtr1=compressed_buffer2;
		    uint64 **tmpPtr2=buffer2;
		    uint64 **tmpPtr3=offset2;
		    compressed_buffer2=compressed_buffer;
		    buffer2=buffer;
		    offset2=offset;
		    compressed_buffer=tmpPtr1;
		    buffer=tmpPtr2;
		    offset=tmpPtr3;
		    buf1status="empty";
		    buf2status="full";
		    buf2size=buf1size;
		    assert(buf2size<=chunkBufSize);
		    buf1size=0;
		  }
		}else{
		  usleep(10);
		}
//		cout<<"0\t"<<now_chunk[0]<<"\t"<<nchunks<<endl;
	  }
	  while(1){
		if(buf2status=="empty"){
		  if(buf1status=="full"){
		    uint64 **tmpPtr1=compressed_buffer2;
		    uint64 **tmpPtr2=buffer2;
		    uint64 **tmpPtr3=offset2;
		    compressed_buffer2=compressed_buffer;
		    buffer2=buffer;
		    offset2=offset;
		    compressed_buffer=tmpPtr1;
		    buffer=tmpPtr2;
		    offset=tmpPtr3;
		    buf1status="empty";
		    buf2status="full";
		    buf2size=buf1size;
		    assert(buf2size<=chunkBufSize);
		    buf1size=0;
		  }else{
		    break;
		  }
		  usleep(100);
		}
	  }
#ifdef PRINT_INFO

	  printf("chunk num is %d   nchunks is %d\n", chunk_num, nchunks);
#endif
	  libdeflate_free_decompressor(decompressor);

	  H5Dclose(datasetID);
	  H5Fclose(fileID);
	}
	if (num_id == 1) {
	  // HashTable
	  now_chunk[1] = 0;
#ifdef OPT_MEMORY_IN_H5
#endif
	  Position1 position;
#ifdef OPT_MEMORY_PLAN_A
	  uint64_t coor=0;
#endif
	  while(now_chunk[1] < nchunks){
		if(buf2status=="full") {
		  for (uint64_t i = 0; i < buf2size; i++) {
			for (uint32_t y = offset2[i][0]; y < min(offset2[i][0] + chunk_dims[0], dims[0]); y++) {
			  for (uint32_t x = offset2[i][1]; x < min(offset2[i][1] + chunk_dims[1], dims[1]); x++) {
				position.x = x;
				position.y = y;
				position.count = 0;
				position.fullTag = 0;
				if(rank>=3) {
				  for (uint32_t z = offset2[i][2]; z < min(offset2[i][2] + chunk_dims[2], dims[2]); z++) {
					uint64_t barcodeInt = buffer2[i][(y - offset2[i][0])*chunk_dims[1]*dims[2] + (x - offset2[i][1])*dims[2] + z];
					if(barcodeInt==0) {
					  continue;
					}
#ifdef TMPD
				  if(position.x>=35500 && position.x<35700 && position.y>=20500 && position.y<20600){
				    cout<<position.x<<"\t"<<position.y<<endl;
				  }
#endif
					bpMap[barcodeInt] = position;
				  }
				}else{
				  uint64_t barcodeInt = buffer2[i][(y - offset2[i][0])*chunk_dims[1] + x - offset2[i][1]];
				  if(barcodeInt==0) {
					continue;
				  }
				  bpMap[barcodeInt] = position;
				}
			  }
			}
			now_chunk[1]++;
		  }
		  buf2status = "empty";
		} else {
		  usleep(10);
		}
		
//		cout<<"1\t"<<now_chunk[1]<<"\t"<<nchunks<<endl;
	  }
	}
  }
  uint64_t chunkSize=0;
#ifdef OPT_MEMORY_IN_H5
  chunkSize=chunkBufSize;
#else
  chunkSize=nchunks;
#endif
  for (uint32_t i = 0; i < chunkSize; i++) {
	delete[] compressed_buffer[i];
	delete[] buffer[i];
	delete[] offset[i];
#ifdef OPT_MEMORY_IN_H5
	delete[] compressed_buffer2[i];
	delete[] buffer2[i];
	delete[] offset2[i];
#endif
  }
  delete[] chunk_size;
  delete[] compressed_buffer;
  delete[] buffer;
  delete[] offset;
#ifdef OPT_MEMORY_IN_H5
  delete[] compressed_buffer2;
  delete[] buffer2;
  delete[] offset2;
#endif
  return 0;
//  mapSize = bpmap_num;
}

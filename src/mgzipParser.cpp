/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
//
// Created by berry on 2021/12/13.
//

#include "mgzipParser.h"

mgzipParser::mgzipParser(string f1, string f2, int threadIdx, int totalThreadsNum)
{
	this->threadIdx = threadIdx;
	if (totalThreadsNum <= 0)
	{
		cerr << "Error, threads number should be not less than 1" << endl;
		exit(1);
	}
	if (!f1.empty())
	{
		if (mgzf_is_mgzf(f1.c_str()))
		{
			r1 = mgzf_open(f1.c_str(), "r");
		}
		else
		{
			cerr << "Error, not a mgzip format," << f1 << endl;
			exit(1);
		}
	}
	else
	{
		r1 = nullptr;
	}
	if (!f2.empty())
	{
		if (mgzf_is_mgzf(f2.c_str()))
		{
			r2 = mgzf_open(f2.c_str(), "r");
		}
		else
		{
			cerr << "Error, not a mgzip format," << f2 << endl;
			exit(1);
		}
	}
	else
	{
		r2 = nullptr;
	}
	if (r1 != nullptr)
	{
		mgzf_index_build_init(r1);
		if (r2 != nullptr)
		{
			mgzf_index_build_init(r2);
		}
		if (mgzf_index_add(r1) == 0)
		{
			if (r2 != nullptr)
			{
				if (mgzf_index_add(r2) == 0)
				{
					if (r1->idx->noffs != r2->idx->noffs)
					{
						string errMsg="different blocks number in r1 and r2";
						sawErrCode(err_fileParse_failed,errMsg);
						cerr << "Error, different blocks number in r1 and r2" << endl;
						exit(1);
					}
				}
			}
			int blocks = r1->idx->noffs;
			int assignedBlocks = blocks % totalThreadsNum == 0 ? blocks / totalThreadsNum : blocks / totalThreadsNum + 1;
			startBlockIdx = threadIdx * assignedBlocks;
			endBlockIdx = startBlockIdx + (uint)assignedBlocks > (uint)blocks ? (uint)blocks : startBlockIdx + (uint)assignedBlocks;
			//	  posLock.lock();
			//	  cout<<threadIdx<<"\t"<<startBlockIdx<<"\t"<<endBlockIdx<<"\t"<<r1->idx->noffs<<endl;
			//	  posLock.unlock();
			if (startBlockIdx >= endBlockIdx)
			{
				r1 = nullptr;
				r2 = nullptr;
				return;
			}
			curBlockIdx = startBlockIdx;
			r1Block = block_read_init(r1->idx->offs[curBlockIdx].raw_size);
			if (r2 != nullptr)
			{
				r2Block = block_read_init(r2->idx->offs[curBlockIdx].raw_size);
			}
			bufIdx1 = 0;
			bufIdx2 = 0;
			if (moveToNextBlock() < 0)
			{
				string errMsg="cannot get block information from mgzip file";
				sawErrCode(err_fileParse_failed,errMsg);
				cerr << "Error, cannot get block information from mgzip file," << f1 << endl;
				exit(1);
			}
		}
		else
		{
			cerr << "Error, cannot get block information from mgzip file," << f1 << endl;
			exit(1);
		}
		//	mgzf_close(r1);
	}
}
int mgzipParser::readHeaderOfBlock(MGZF *fp)
{
	// ? means an option
	// 1f 8b method flag mtime[4 bytes) XFL OS XLEN[2 bytes] XETRAN[XLEN bytes] FNAME?[variable length, 00 terminated] COMMENT?[variable length, 00 terminated] CRC16?[2]
	/* flag
	  bit 0   FTEXT
	  bit 1   FHCRC
	  bit 2   FEXTRA
	  bit 3   FNAME
	  bit 4   FCOMMENT
	  bit 5   reserved
	  bit 6   reserved
	  bit 7   reserved
	 */
	int headerSize = 0;
	int count;
	uint8_t buf[10];
	count = _mgzf_read(fp->fp, buf, 10);
	headerSize += 10;
	if (count != 10)
	{ // no data read
		return -1;
	}

	if ((buf[3] & FEXTRA) == FEXTRA)
	{
		uint16_t xlen[1];
		count = _mgzf_read(fp->fp, xlen, 2);
		headerSize += 2;
		uint8_t *extraBuf = new uint8_t[xlen[0]];
		count = _mgzf_read(fp->fp, extraBuf, xlen[0]);
		headerSize += xlen[0];
		delete[] extraBuf;
	}
	if ((buf[3] & FNAME) == FNAME)
	{
		uint8_t tmpBuf[1];
		while (true)
		{
			_mgzf_read(fp->fp, tmpBuf, 1);
			headerSize++;
			if (tmpBuf[0] == 0)
			{
				break;
			}
		}
	}
	if ((buf[3] & FCOMMENT) == FCOMMENT)
	{
		uint8_t tmpBuf[1];
		while (true)
		{
			_mgzf_read(fp->fp, tmpBuf, 1);
			headerSize++;
			if (tmpBuf[0] == 0)
			{
				break;
			}
		}
	}
	if ((buf[3] & FHCRC) == FHCRC)
	{
		uint8_t tmpBuf[2];
		_mgzf_read(fp->fp, tmpBuf, 2);
		headerSize += 2;
	}
	return headerSize;
}
int mgzipParser::inflateFromBlock(MgzBlock *fp, uint8_t *compressedData, int block_length)
{
	z_stream zs;
	zs.zalloc = NULL;
	zs.zfree = NULL;
	zs.next_in = compressedData;
	zs.avail_in = block_length;
	zs.next_out = (Bytef *)fp->uncompressed_block;
	zs.avail_out = fp->raw_length;

	if (inflateInit2(&zs, -15) != Z_OK)
	{
		//	fp->fp->errcode |= MGZF_ERR_ZLIB;
		return -1;
	}
	if (inflate(&zs, Z_FINISH) != Z_STREAM_END)
	{
		inflateEnd(&zs);
		//	fp->errcode |= MGZF_ERR_ZLIB;
		return -1;
	}
	if (inflateEnd(&zs) != Z_OK)
	{
		//	fp->errcode |= MGZF_ERR_ZLIB;
		return -1;
	}
	return zs.total_out;
}
int mgzipParser::moveToNextBlock()
{

	//  posLock.lock();
	//  cout<<threadIdx<<"\t"<<curBlockIdx<<endl;
	//  posLock.unlock();
	if (curBlockIdx < endBlockIdx)
	{
		block_destroy(r1Block);
		r1->raw_length = r1->idx->offs[curBlockIdx].raw_size;
		r1->block_address = r1->idx->offs[curBlockIdx].start;
		r1->block_length = r1->idx->offs[curBlockIdx].block_size;
		if (_mgzf_seek(r1->fp, r1->block_address, SEEK_SET) < 0)
		{
			string errMsg = "mgzf seek error";
			sawErrCode(err_fileParse_failed, errMsg);
			cerr << "Error, mgzf seek error" << endl;
			exit(1);
		}
		if ((r1Block = readBlock(r1)) == nullptr)
		{
			string errMsg = "mgzf seek error";
			sawErrCode(err_fileParse_failed, errMsg);
			cerr << "Error, mgzf read error" << endl;
			exit(1);
		}
		if (r2 != nullptr)
		{
			block_destroy(r2Block);
			r2->raw_length = r2->idx->offs[curBlockIdx].raw_size;
			r2->block_address = r2->idx->offs[curBlockIdx].start;
			r2->block_length = r2->idx->offs[curBlockIdx].block_size;
			if (_mgzf_seek(r2->fp, r2->block_address, SEEK_SET) < 0)
			{
				string errMsg = "mgzf seek error";
				sawErrCode(err_fileParse_failed, errMsg);
				cerr << "Error, mgzf seek error" << endl;
				exit(1);
			}
			if ((r2Block = readBlock(r2)) == nullptr)
			{
				string errMsg = "mgzf read error";
				sawErrCode(err_fileParse_failed, errMsg);
				cerr << "Error, mgzf read error" << endl;
				exit(1);
			}
		}
		bufIdx1 = 0;
		bufIdx2 = 0;
		curBlockIdx++;
		return 0;
	}
	else
	{
		return -1;
	}
}
MgzBlock *mgzipParser::readBlock(MGZF *fp)
{
	MgzBlock *curBlock = block_read_init(fp->idx->offs[curBlockIdx].raw_size);
	// read header of a block
	int headerSize = readHeaderOfBlock(fp);
	//  uint8_t header[headerSize];
	uint64_t remaining = fp->idx->offs[curBlockIdx].block_size - headerSize;
	uint8_t *compressed_block = new uint8_t[remaining];

	size_t count = _mgzf_read(fp->fp, compressed_block, remaining);
	if (count != remaining)
	{
		fp->errcode |= MGZF_ERR_IO;
		return nullptr;
	}
	int count2 = inflateFromBlock(curBlock, compressed_block, remaining);
	if (count2 < 0)
		return nullptr;
	delete[] compressed_block;
	return curBlock;
}
int mgzipParser::getReadsFromBlock(ReadPair **readsBuf, int readsNumber)
{
	if (r1 == nullptr)
	{
		return 0;
	}
	int storedReadsNum1 = 0;
	int storedReadsNum2 = 0;
	string readName, readSeq, readStrand, readQuality;
	int newlineNum = 0;
	while (storedReadsNum1 < readsNumber)
	{
		if ((uint)bufIdx1 < r1Block->raw_length)
		{
			if (((char *)r1Block->uncompressed_block)[bufIdx1] != '\n')
			{
				switch (newlineNum)
				{
				case 0:
					readName += ((char *)r1Block->uncompressed_block)[bufIdx1];
					break;
				case 1:
					readSeq += ((char *)r1Block->uncompressed_block)[bufIdx1];
					break;
				case 2:
					readStrand += ((char *)r1Block->uncompressed_block)[bufIdx1];
					break;
				case 3:
					readQuality += ((char *)r1Block->uncompressed_block)[bufIdx1];
					break;
				default:
					break;
				}
				bufIdx1++;
			}
			else
			{
				bufIdx1++;
				newlineNum++;
				if (newlineNum == 4)
				{
					readsBuf[storedReadsNum1]->mLeft = new Read(readName, readSeq, readStrand, readQuality);
					readName = "";
					readSeq = "";
					readStrand = "";
					readQuality = "";
					storedReadsNum1++;
					newlineNum = 0;
					if (storedReadsNum1 == readsNumber)
					{
						break;
					}
				}
			}
		}
		else
		{
			if (r2 == nullptr)
			{
				if (moveToNextBlock() < 0)
				{
					break;
				}
			}
			break;
		}
	}
	readName = "";
	readSeq = "";
	readStrand = "";
	readQuality = "";
	if (r2 != nullptr)
	{
		while (storedReadsNum2 < readsNumber)
		{
			if ((uint)bufIdx2 < r2Block->raw_length)
			{
				if (((char *)r2Block->uncompressed_block)[bufIdx2] != '\n')
				{
					switch (newlineNum)
					{
					case 0:
						readName += ((char *)r2Block->uncompressed_block)[bufIdx2];
						break;
					case 1:
						readSeq += ((char *)r2Block->uncompressed_block)[bufIdx2];
						break;
					case 2:
						readStrand += ((char *)r2Block->uncompressed_block)[bufIdx2];
						break;
					case 3:
						readQuality += ((char *)r2Block->uncompressed_block)[bufIdx2];
						break;
					default:
						break;
					}
					bufIdx2++;
				}
				else
				{
					bufIdx2++;
					newlineNum++;
					if (newlineNum == 4)
					{
						readsBuf[storedReadsNum2]->mRight = new Read(readName, readSeq, readStrand, readQuality);
						readName = "";
						readSeq = "";
						readStrand = "";
						readQuality = "";
						storedReadsNum2++;
						newlineNum = 0;
						if (storedReadsNum2 == readsNumber)
						{
							if (storedReadsNum1 != storedReadsNum2)
							{
								cerr << "Error, reads number are different between r1 and r2" << endl;
								string errMsg="reads number are different between r1 and r2";
								sawErrCode(err_fileParse_failed,errMsg);
								exit(1);
							}
							else
							{
								if (bufIdx2 == r2->raw_length)
								{
									if (moveToNextBlock() < 0)
									{
										break;
									}
									break;
								}
								return storedReadsNum1;
							}
						}
					}
				}
			}
			else
			{
				if (moveToNextBlock() < 0)
				{
					break;
				}
				break;
			}
		}
		if (storedReadsNum1 != storedReadsNum2)
		{
			cerr << "Error, reads number are different between r1 and r2" << endl;
			string errMsg="reads number are different between r1 and r2";
			sawErrCode(err_fileParse_failed,errMsg);
			exit(1);
		}
	}
	else
	{
		if (bufIdx1 == r1->raw_length)
		{
			if (moveToNextBlock() < 0)
			{
				return 0;
			}
		}
		return storedReadsNum1;
	}
	return storedReadsNum1;
}
mgzipParser::~mgzipParser()
{
	if (r1 != nullptr)
	{
		mgzf_close(r1);
	}
	if (r2 != nullptr)
	{
		mgzf_close(r2);
	}
}

#ifndef MGZIP_MGZF_H
#define MGZIP_MGZF_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#ifdef _USE_KNETFILE
#include "knetfile.h"
typedef knetFile *_mgzf_file_t;
#define _mgzf_open(fn, mode) knet_open(fn, mode)
#define _mgzf_dopen(fp, mode) knet_dopen(fp, mode)
#define _mgzf_close(fp) knet_close(fp)
#define _mgzf_fileno(fp) ((fp)->fd)
#define _mgzf_tell(fp) knet_tell(fp)
#define _mgzf_seek(fp, offset, whence) knet_seek(fp, offset, whence)
#define _mgzf_read(fp, buf, len) knet_read(fp, buf, len)
#define _mgzf_write(fp, buf, len) knet_write(fp, buf, len)
#else // ~defined(_USE_KNETFILE)
#if defined(_WIN32) || defined(_MSC_VER)
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else // ~defined(_WIN32)
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif // ~defined(_WIN32)
typedef FILE *_mgzf_file_t;
#define _mgzf_open(fn, mode) fopen(fn, mode)
#define _mgzf_dopen(fp, mode) fdopen(fp, mode)
#define _mgzf_close(fp) fclose(fp)
#define _mgzf_fileno(fp) fileno(fp)
#define _mgzf_tell(fp) ftello(fp)
#define _mgzf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define _mgzf_read(fp, buf, len) fread(buf, 1, len, fp)
#define _mgzf_write(fp, buf, len) fwrite(buf, 1, len, fp)
#endif // ~define(_USE_KNETFILE)


#define MGZF_BLOCK_SIZE 0x20000000 // 512M
#define MGZF_BLOCK_SEQ_SIZE 0x100000 // 1M, the quantities of reads

#define BLOCK_HEADER_LENGTH 20
#define BLOCK_FOOTER_LENGTH 8

#define MGZF_ERR_ZLIB   1
#define MGZF_ERR_HEADER 2
#define MGZF_ERR_IO     4
#define MGZF_ERR_MISUSE 8

static const uint8_t end_block[32] = "\037\213\010\4\0\0\0\0\0\377\010\0\111\107\4\0\037\0\0\0\0\3\0\0\0\0\0\0\0\0\0";

typedef struct mgzidx_t mgzidx_t;
typedef struct
{
//    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t start;
    uint64_t block_size;
    uint64_t raw_size;
}mgzidx1_t;

struct mgzidx_t
{
    int noffs, moffs;       // the size of the index, n:used, m:allocated
    mgzidx1_t *offs;        // offsets
    uint64_t ublock_addr;   // offset of the current block (uncompressed data)
};

typedef struct {
    int open_mode:8, compress_level:8, errcode:16;
    int cache_size;
    int block_length, block_offset;
    int raw_length;
    int64_t block_address;
    void *uncompressed_block, *compressed_block;
    void *cache; // a pointer to a hash table
  	FILE *fp; // actual file handler; FILE* on writing; FILE* or knetFile* on reading
    mgzidx_t *idx;      // MGZF index
    int idx_build_otf;  // build index on the fly, set by mgzf_index_build_init()
} MGZF;
typedef struct{
//  int block_length;
  uint64_t raw_length;
  void *uncompressed_block;
} MgzBlock;


#ifndef KSTRING_T_MGZ
#define KSTRING_T_MGZ kstring_t_mgz
typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t_mgz;
#endif

#ifdef __cplusplus
extern "C" {
#endif

MgzBlock* block_read_init(uint64_t length);
int block_destroy(MgzBlock* b);
int inflate_FromBlock(MgzBlock* fp, void* compressedData,int block_length);
void mgz_error(const char *format, ...);

/******************
 * Basic routines *
 ******************/

/**
 * Open an existing file descriptor for reading or writing.
 *
 * @param fd    file descriptor
 * @param mode  mode matching /[rwu0-9]+/: 'r' for reading, 'w' for writing and a digit specifies
 *              the zlib compression level; if both 'r' and 'w' are present, 'w' is ignored.
 * @return      MGZF file handler; 0 on error
 */
MGZF* mgzf_dopen(int fd, const char *mode);

/**
 * Open the specified file for reading or writing.
 */
MGZF* mgzf_open(const char* path, const char *mode);

/**
 * Close the MGZF and free all associated resources.
 *
 * @param fp    MGZF file handler
 * @return      0 on success and -1 on error
 */
int mgzf_close(MGZF *fp);

/**
 * Read up to _length_ bytes from the file storing into _data_.
 *
 * @param fp     MGZF file handler
 * @param data   data array to read into
 * @param length size of data to read
 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
 */
ssize_t mgzf_read(MGZF *fp, void *data, ssize_t length);

/**
 * Write _length_ bytes from _data_ to the file.
 *
 * @param fp     MGZF file handler
 * @param data   data array to write
 * @param length size of data to write
 * @return       number of bytes actually written; -1 on error
 */
ssize_t mgzf_write(MGZF *fp, const void *data, ssize_t length);

/**
 * Write the data in the buffer to the file.
 */
int mgzf_flush(MGZF *fp);

/**
 * Return a virtual file pointer to the current location in the file.
 * No interpetation of the value should be made, other than a subsequent
 * call to mgzf_seek can be used to position the file at the same point.
 * Return value is non-negative on success.
 */
#define mgzf_tell(fp) ((fp->block_address << 16) | (fp->block_offset & 0xFFFF))

/**
 * Set the file to read from the location specified by _pos_.
 *
 * @param fp     MGZF file handler
 * @param pos    virtual file offset returned by mgzf_tell()
 * @param whence must be SEEK_SET
 * @return       0 on success and -1 on error
 */
int64_t mgzf_seek(MGZF *fp, int64_t pos, int whence);

/**
 * Check if the MGZF end-of-file (EOF) marker is present
 *
 * @param fp    MGZF file handler opened for reading
 * @return      1 if EOF is present; 0 if not or on I/O error
 */
int mgzf_check_EOF(MGZF *fp);

/**
 * Check if a file is in the MGZF format
 *
 * @param fn    file name
 * @return      1 if _fn_ is MGZF; 0 if not or on I/O error
 */
int mgzf_is_mgzf(const char *fn);

/*********************
 * Advanced routines *
 *********************/

/**
 * Flush the file if the remaining buffer size is smaller than _size_
 */
int mgzf_flush_try(MGZF *fp, ssize_t size);

/**
 * Read one byte from a MGZF file. It is faster than mgzf_read()
 * @param fp     MGZF file handler
 * @return       byte read; -1 on end-of-file or error
 */
int mgzf_getc(MGZF *fp);

/**
 * Read one line from a MGZF file. It is faster than mgzf_getc()
 *
 * @param fp     MGZF file handler
 * @param delim  delimitor
 * @param str    string to write to; must be initialized
 * @return       length of the string; 0 on end-of-file; negative on error
 */
int mgzf_getline(MGZF *fp, int delim, kstring_t_mgz *str);

/**
 * Read the next MGZF block.
 */
int mgzf_read_block(MGZF *fp);

int mgzf_read_block2(MGZF *fp);

int deflate_block(MGZF* fp, int block_length);

int inflate_block(MGZF* fp, int block_length);


/**
  * Tell MGZF to build index while compressing.
  *
  * @param fp          MGZF file handler; can be opened for reading or writing.
  *
  * Returns 0 on success and -1 on error.
  *
  * @note This function must be called before any data has been read or
  * written, and in particular before calling mgzf_mt() on the same
  * file handle (as threads may start reading data before the index
  * has been set up).
  */
int mgzf_index_build_init(MGZF *fp);

/// Load MGZF index
/**
 * @param fp          MGZF file handler
 * @param bname       base name
 * @param suffix      suffix to add to bname (can be NULL)
 * @return 0 on success and -1 on error.
 */
int mgzf_index_load(MGZF *fp, const char *bname, const char *suffix);

void mgzf_index_destroy(MGZF *fp);

/**
 * add index info to MGZF *fp
 * @param fp
 * @return 0 on success and -1 on error.
 */
int mgzf_index_add(MGZF *fp);

/**
 * write index to a file
 * @param fp          MGZF file handler
 * @param bname       base name
 * @param suffix      suffix to add to bname (can be NULL)
 * @return 0 on success and -1 on error.
 */
int mgzf_index_dump(MGZF *fp, const char *bname, const char *suffix);

/**
 * get index from a mgz file and write it to a file
 * @param mgz_file    MGZF file
 * @param bname       base name
 * @param suffix      suffix to add to bname (can be NULL)
 * @return 0 on success and -1 on error.
 */
int mgzf_index_dump_fromfile(const char *mgz_file, const char *bname, const char *suffix);

int check_header(const uint8_t *header);

MGZF *mgzf_read_init();

MGZF *mgzf_write_init(int compress_level);

int mgzf_destroy(MGZF* fp);

#ifdef __cplusplus
}
#endif

#endif //MGZIP_MGZF_H

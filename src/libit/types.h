#ifndef _TYPES_H_
#define _TYPES_H_

#ifndef _ITLIB_TYPES_
#define _ITLIB_TYPES_
typedef char       			int8;
typedef unsigned char       uint8;
typedef short int 			int16;
typedef unsigned short int	uint16;
typedef unsigned int       	uint32;
typedef unsigned long long 	uint64;
#endif

/**
    \brief The compression types.
 */
typedef enum {
    LOSS,		///Loss compression.
    LOSSLESS	///Lossless compression.
} Compression;

/**
    \brief The image format.
 */
typedef enum {
    CS420,
    CS444,
    RGB,
    RGBY,
    GREY,
    BAYER
} ColorSpace;

/**
    \brief The Decorrelation  method.
 */
typedef enum {
    WAVELET,
    PREDICT,
    RESIZE,
    VECTORIZE
} Decorrelation;

/**
    \brief The gamma correction method.
 */
typedef enum {
    LINEAR,
    BT709,
    sRGB
} Gamma;

/**
    \brief The frame type.
 */
typedef enum {
    I_FRAME = 0,
    P_FRAME = 1,
    B_FRAME = 2
} FrameType;

/**
    \brief The image transform states.
 */
typedef enum {
    FRAME_COPY 		=	1,
    BUFFER_READ		=	1<<1,
    DWT 			= 	1<<2,
    IDWT 			= 	1<<3,
    FILL_SUBBAND 	= 	1<<4,
    QUANTIZATION	= 	1<<5,
    RANGE_ENCODER	= 	1<<6,
    RANGE_DECODER	= 	1<<7,
    BITS_ALLOCATION	= 	1<<8,
    MEDIAN_FILTER	=	1<<9,
    COLOR_TRANSFORM = 	1<<10,
    FILL_HISTOGRAM 	= 	1<<11,
    WHITE_BALANCE	= 	1<<12,
    SEGMENTATION	= 	1<<13,
    MATCH			= 	1<<14
}CodecState;

/**
    \brief The bayer grid pattern.

    All RGB cameras use one of these Bayer grids:\n
*/
/*
    BGGR  0         GRBG 1          GBRG  2         RGGB 3
      0 1 2 3 4 5	  0 1 2 3 4 5	  0 1 2 3 4 5	  0 1 2 3 4 5
    0 B G B G B G	0 G R G R G R	0 G B G B G B	0 R G R G R G
    1 G R G R G R	1 B G B G B G	1 R G R G R G	1 G B G B G B
    2 B G B G B G	2 G R G R G R	2 G B G B G B	2 R G R G R G
    3 G R G R G R	3 B G B G B G	3 R G R G R G	3 G B G B G B
 */
typedef enum {
    BGGR = 0,
    GRBG = 1,
    GBRG = 2,
    RGGB = 3
} BayerGrid;

/**
    \brief Wavelet transform filter banks.
 */
typedef enum {
    FR_HAAR	= 0,
    FR_5_3	= 1
} FilterBank;

/**
    \brief The range coder type
 */
typedef enum {
    ADAP	= 0,	///Adaptive range coder.
    NADAP	= 1,	///Nonadaptive range coder.
    FAST	= 2     ///Nonadaptive fast range coder.
} RangeType;

#endif //_WALET_HH_

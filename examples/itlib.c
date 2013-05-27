#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libit/types.h"

static void Help(void) {
    printf("Usage: dwebp in_file [options] [-o out_file]\n\n"
           "Decodes the WebP image file to PNG format [Default]\n"
           "Use following options to convert into alternate image formats:\n"
           "  -pam ......... save the raw RGBA samples as a color PAM\n"
           "  -ppm ......... save the raw RGB samples as a color PPM\n"
           "  -pgm ......... save the raw YUV samples as a grayscale PGM\n"
           "                 file with IMC4 layout.\n"
           "  -yuv ......... save the raw YUV samples in flat layout.\n"
           "\n"
           " Other options are:\n"
           "  -version  .... print version number and exit.\n"
           "  -nofancy ..... don't use the fancy YUV420 upscaler.\n"
           "  -nofilter .... disable in-loop filtering.\n"
           "  -mt .......... use multi-threading\n"
           "  -crop <x> <y> <w> <h> ... crop output with the given rectangle\n"
           "  -scale <w> <h> .......... scale the output (*after* any cropping)\n"
           "  -alpha ....... only save the alpha plane.\n"
           "  -h     ....... this help message.\n"
           "  -v     ....... verbose (e.g. print encoding/decoding times)\n"
           );
}

int main(int argc, const char *argv[]) {

    const char *in_file = NULL;
    const char *out_file = NULL;
    int i;
    uint8 *in;
    int16 *out;

    for (i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
            Help();
            return 0;
        } else if (!strcmp(argv[i], "-o") && i < argc - 1) {
            out_file = argv[++i];
        } else if (!strcmp(argv[i], "-alpha")) {
            utils_image_copy_n(in, out, 320, 200, 8);

        }
    }

    if (in_file == NULL) {
        fprintf(stderr, "missing input file!!\n");
        Help();
        return -1;
    }


    if (out_file) {
        printf("Decoded %s. Now saving...\n", in_file);
        //SaveOutput(output_buffer, format, out_file);
    } else {
        printf("File %s can be decoded \n",in_file);
        printf("Nothing written; use -o flag to save the result as e.g. PNG.\n");
    }

    return 0;
}

//------------------------------------------------------------------------------

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <png.h>

#include "../libit/types.h"
#include "./itlib.h"

static int verb = 0;

static void get_first_pixels(uint8 * buff, int num, int bpp)
{
    int16 *pic;
    int i;

    if(bpp > 1){
        pic = (int16*)buff;
        for(i=0; i < num; i++){
            printf("%3d %3d %4d\n", (pic[i] & 0xFF00)>>8, pic[i] & 0x00FF,  pic[i]);
        }
    } else {
        for(i=0; i < num; i++){
            printf("%4d\n", buff[i]);
        }
    }
}

static int get_trans_status(const TransState* ts, const TransState* ts1)
{
    if(ts->colort == ts1->colort && ts->bpp == ts1->bpp) return 0;
    return 1;
}

static void PNGAPI error_function(png_structp png, png_const_charp dummy)
{
    (void)dummy;  // remove variable-unused warning
    longjmp(png_jmpbuf(png), 1);
}

/** \brief Read PNG file to the buffer.
    \param in_file  The input file descriptor.
    \param buff		The input buffer.
    \param w		The pointer to image width.
    \param h		The pointer to image height.
    \param bpp		The pointer to bits per pixel.
    \param colort	The pointer to color types.
    \retval         The return value, if 0 - OK.
*/
int readPNG(FILE* in_file, uint8** buff, int* const w, int* const h, int* const bpp, int* const colort)
{
    png_structp png;
    png_infop info = NULL, end_info = NULL;
    png_bytep row;
    int interlaced, has_alpha, num_passes, stride, p, ok=0;

    png_uint_32 y;

    png = png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    if (png == NULL) {
        fprintf(stderr, "Error! readPNG: Can't create png struct\n");
        ok = 1; goto End;
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error! readPNG: Can't create png_jmpbuf\n");
        ok = 1; goto End;
    }

    info = png_create_info_struct(png);
    if (info == NULL) {
        fprintf(stderr, "Error! readPNG: Can't create info struct\n");
        ok = 1; goto End;
    }

    end_info = png_create_info_struct(png);
    if (end_info == NULL) {
        fprintf(stderr, "Error! readPNG: Can't create end_info struct\n");
        ok = 1; goto End;
    }

    png_init_io(png, in_file);
    png_read_info(png, info);
    if (!png_get_IHDR(png, info, w, h, bpp, colort, &interlaced, NULL, NULL)) {
        fprintf(stderr, "Error! readPNG: Can't create png_get_IHDR\n");
        ok = 1; goto End;
    }
    //if(verb) printf("readPNG: w = %d h = %d bpp = %d colort = %d\n", *w, *h, *bpp, *colort);

    png_set_strip_16(png);
    png_set_packing(png);
    if (*colort == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
    if (*colort == PNG_COLOR_TYPE_GRAY || *colort == PNG_COLOR_TYPE_GRAY_ALPHA) {
        if (*bpp < 8) {
            png_set_expand_gray_1_2_4_to_8(png);
        }
        png_set_gray_to_rgb(png);
    }
    if (png_get_valid(png, info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(png);
        has_alpha = 1;
    } else {
        has_alpha = !!(*colort & PNG_COLOR_MASK_ALPHA);
    }

    num_passes = png_set_interlace_handling(png);
    png_read_update_info(png, info);
    //Bits per pixel to bytes per pixel
    *bpp = (*bpp > 8) ? 2 : 1;
    stride = (has_alpha ? 4 : 3) * (*bpp) * (*w);

    *buff = (uint8*)malloc(stride * (*h));
    if (*buff == NULL) {
        fprintf(stderr, "Error! readPNG: Can't allocate memory\n");
        ok = 1; goto End;
    }

    for (p = 0; p < num_passes; ++p) {
        for (y = 0; y < *h; y++) {
            row = (*buff) + y * stride;
            png_read_rows(png, &row, NULL, 1);
        }
    }
    //if(verb) printf("readPNG: Finish reading file\n");
    png_read_end(png, end_info);
    //if(verb) printf("readPNG: Finish png_read_end file\n");
    if(verb)  get_first_pixels(*buff, 10, *bpp);

End:
    png_destroy_read_struct(&png, &info, &end_info);
    //if(verb) printf("readPNG: Finish png_destroy_read_struct file\n");
    //if(*buff) free(*buff);
    return ok;
}

/** \brief Write buffer to PNG file.
    \param out_file  The output file descriptor.
    \param buff		The input buffer.
    \param w		The pointer to image width.
    \param h		The pointer to image height.
    \param bpp		The pointer to bits per pixel.
    \param colort	The pointer to color types.
    \retval         The return value, if 0 - OK.
*/
int writePNG(FILE* out_file, uint8* const buff, const int w, const int h, const int bpp, const int colort)
{
    const int stride =  bpp * w;
    png_structp png;
    png_infop info;
    png_uint_32 y;
    png_bytep row;

    png = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL, error_function, NULL);
    if (png == NULL) {
        fprintf(stderr, "Error! write PNG: Can't create png struct.");
        return 1;
    }
    info = png_create_info_struct(png);
    if (info == NULL) {
        fprintf(stderr, "Error! write PNG: Can't create inf struct.");
        png_destroy_write_struct(&png, NULL);
        return 1;
    }
    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error! write PNG: Can't create png_jmpbuf.");
        png_destroy_write_struct(&png, &info);
        return 1;
    }
    png_init_io(png, out_file);
    png_set_IHDR(png, info, w, h, bpp > 1 ? 8 : 8, colort,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);
    for (y = 0; y < h; y++) {
        row = buff + y * stride;
        png_write_rows(png, &row, 1);
    }
    png_write_end(png, info);
    png_destroy_write_struct(&png, &info);

    if(verb)  get_first_pixels(buff, 10, bpp);

    return 0;
}

/** \brief Read PGM file to the buffer.
    \param in_file  The input file descriptor.
    \param buff		The output buffer.
    \param w		The pointer to image width.
    \param h		The pointer to image height.
    \param bpp		The pointer to bits per pixel.
    \retval         The return value, if 0 - OK.
*/
int readPGM(FILE* in_file, uint8** buff, int* const w, int* const h, int* const bpp)
{
    uint8 line[10];
    int size;

    fscanf(in_file, "%s", line);
    if (strcmp(line, "P5") != 0) {
        fprintf(stderr, "Error! readPGM: It's not PPM file\n");
        return 1;
    }

    fscanf(in_file, "%d%d%d", w, h, bpp);
    *bpp = (*bpp == 255) ? 1 : 2;
    //if(verb) printf("w = %d h = %d bpp = %d\n", *w, *h, *bpp);


    size = (*w)*(*h)*(*bpp);
    //if(verb) printf("readPGM: Read file = %p size = %d\n",in_file, size);

    *buff = (uint8*)malloc(size);
    if (*buff == NULL) {
        fprintf(stderr, "Error! readPGM: Can't allocate memory\n");
        return 1;
    }

    if(fread(*buff, size, 1,  in_file) != 1){
        fprintf(stderr, "Error! readPGM: Image read error\n");
        free(*buff);
        return 1;
    }

    if(verb)  get_first_pixels(*buff, 10, *bpp);

    return 0;
}

/** \brief Write buffer to PGM file.
    \param out_file The output file descriptor.
    \param buff		The input buffer.
    \param w		The pointer to image width.
    \param h		The pointer to image height.
    \param bpp		The pointer to bits per pixel.
    \retval         The return value, if 0 - OK.
*/
int writePGM(FILE* out_file, uint8* buff, const int w, const int h, const int bpp)
{

    fprintf(out_file, "P5\n%d %d\n%d", w, h, (bpp > 1) ? 65535 : 255);

    //if(verb) printf("Write file = %p size = %d\n",out_file, w*h*bpp);

    if (fwrite(buff, w*h, bpp, out_file) != bpp) {
        fprintf(stderr, "Error! writePGM: Write PGM file problem\n");
        return 1;
    }

    if(verb)  get_first_pixels(buff, 10, bpp);

    return 0;
}

int main(int argc, const char *argv[]) {

    const char *in_file = NULL;
    const char *out_file = NULL;
    FILE *IN_FILE, *OUT_FILE;
    int i, ok = 0;
    uint8 *pic8;
    int16 *pic16;
    uint8* buff;
    int w, h, bpp, colort;

    for (i = 1; i < argc; i++)  if (!strcmp(argv[i], "-v")) verb = 1;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
            printf("Usage: itlib in_file [options] [-o out_file]\n\n"
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
            return 0;
        } else if (i == 1){
            //Open input file
            in_file = argv[i];
            IN_FILE = fopen(in_file, "rb");
            if (IN_FILE == NULL) {
                fprintf(stderr, "Error! Cannot open input file '%s'\n", in_file);
                return 1;
            }
            if(!strcmp(&in_file[strlen(in_file)-4],".pgm") || !strcmp(&in_file[strlen(in_file)-4],".PGM")){

                ok = readPGM(IN_FILE, &buff, &w, &h, &bpp);
                if(verb) printf("Read %s file w = %d h = %d bpp = %d\n", in_file, w, h, bpp);

            } else if(!strcmp(&in_file[strlen(in_file)-4],".png") || !strcmp(&in_file[strlen(in_file)-4],".PNG")){

                ok = readPNG(IN_FILE, &buff, &w, &h, &bpp, &colort);
                if(verb) printf("Read %s file w = %d h = %d bpp = %d colort = %d\n", in_file, w, h, bpp, colort);

            } else ok = 1;

            if(IN_FILE) fclose(IN_FILE);
            if(ok){
                fprintf(stderr, "Error! Read input file '%s'\n", in_file);
                return 1;
            }

        } else if (!strcmp(argv[i], "-o") && i < argc - 1) {
            //Open output file
            out_file = argv[++i];
            OUT_FILE = fopen(out_file, "wb");
            if (OUT_FILE == NULL) {
                fprintf(stderr, "Error! Cannot open input file '%s'\n", out_file);
                return 1;
            }
            if(!strcmp(&out_file[strlen(out_file)-4],".pgm") || !strcmp(&out_file[strlen(out_file)-4],".PGM")){

                ok = writePGM(OUT_FILE, buff, w, h, bpp);
                if(verb) printf("Write %s file\n", out_file);

            } else if(!strcmp(&out_file[strlen(out_file)-4],".png") || !strcmp(&out_file[strlen(out_file)-4],".PNG")){

                ok = writePNG(OUT_FILE, buff, w, h, bpp, GREY);
                if(verb) printf("Write %s file\n", out_file);

            } else ok = 1;

            if(OUT_FILE) fclose(OUT_FILE);
            if(ok){
                fprintf(stderr, "Error! Write output file '%s'\n", out_file);
                return 1;
            }

        } else if (!strcmp(argv[i], "-be")) {
            //utils_image_copy_n(in, out, 320, 200, 8);

        } else if (!strcmp(argv[i], "-le")) {
            //utils_image_copy_n(in, out, 320, 200, 8);

        }
    }

    if(argc == 1) fprintf(stderr, "Help: itlib -h\n");
    free(buff);

    return 0;
}

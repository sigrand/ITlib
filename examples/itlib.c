#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <png.h>

#include "../libit/types.h"
#include "./itlib.h"

static int verb = 0;

static void print_first_pixels(uint8 * buff, int num, int bpp, int colort, int w)
{
    int16 *pic;
    int i;

    if(colort == GREY || colort == BAYER){
        if(bpp > 8){
            pic = (int16*)buff;
            for(i=0; i < num; i++) printf("%4d ", pic[i]>>4);
            printf("\n");
            for(i=0; i < num; i++) printf("%4d ", pic[w+i]>>4);
            printf("\n");
        } else {
            for(i=0; i < num; i++) printf("%3d ", buff[i]);
            printf("\n");
            for(i=0; i < num; i++) printf("%3d ", buff[w+i]);
            printf("\n");
        }
    } else {
        for(i=0; i < num*3; i+=3) printf("%3d %3d %3d   ", buff[i], buff[i+1], buff[i+2]);
        printf("\n");
        for(i=0; i < num*3; i+=3) printf("%3d %3d %3d   ", buff[w*3+i], buff[w*3+i+1], buff[w*3+i+2]);
        printf("\n");
    }
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
    //*bpp = (*bpp > 8) ? 2 : 1;
    stride = (has_alpha ? 4 : 3) * 2 * (*w);

    //Create buffer up to 2 times more the bayer or gray image for the future used
    *buff = (uint8*)malloc(stride*(*h)*2);
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
    if(verb)  print_first_pixels(*buff, 10, *bpp, *colort, *w);

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
    int stride, color;
    png_structp png;
    png_infop info;
    png_uint_32 y;
    png_bytep row;

    if(colort == GREY || colort == BAYER) {
        color = PNG_COLOR_TYPE_GRAY;
        stride =  ((bpp > 8) ? 2 : 1) * w;
    } else {
        color = PNG_COLOR_TYPE_RGB;
        stride =  ((bpp > 8) ? 2 : 1) * w * 3;

    }
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

    //printf("color = %d  %d\n", color, PNG_COLOR_TYPE_RGB);
    png_set_IHDR(png, info, w, h, 8, color,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);
    for (y = 0; y < h; y++) {
        row = buff + y * stride;
        png_write_rows(png, &row, 1);
    }
    png_write_end(png, info);
    png_destroy_write_struct(&png, &info);

    if(verb)  print_first_pixels(buff, 10, bpp, colort, w);

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

    *bpp = (*bpp == 255) ? 8 : 16;
    //if(verb) printf("w = %d h = %d bpp = %d\n", *w, *h, *bpp);

    size = (*w)*(*h)*2;
    //if(verb) printf("readPGM: Read file = %p size = %d\n",in_file, size);

    //Create buffer up to 6 times more the bayer or gray image for the future used
    *buff = (uint8*)malloc(size*3*2);
    if (*buff == NULL) {
        fprintf(stderr, "Error! readPGM: Can't allocate memory\n");
        return 1;
    }

    fread(*buff, 1, 1, in_file);

    if(fread(*buff, size, 1,  in_file) != 1){
        fprintf(stderr, "Error! readPGM: Image read error\n");
        free(*buff);
        return 1;
    }

    if(verb)  print_first_pixels(*buff, 10, 8, GREY, *w);

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

    fprintf(out_file, "P5\n%d %d\n%d", w, h, (bpp > 8) ? 65535 : 255);

    //if(verb) printf("Write file = %p size = %d\n",out_file, w*h*bpp);
    if (fwrite(buff, w*h, (bpp > 8) ? 2 : 1, out_file) != ((bpp > 8) ? 2 : 1)) {
        fprintf(stderr, "Error! writePGM: Write PGM file problem\n");
        return 1;
    }

    if(verb)  print_first_pixels(buff, 10, bpp, GREY, w);

    return 0;
}

void copy_image16(int16 *in, int16 *out, const int w, const int h)
{
    int i, size = w*h;
    for(i=0; i < size; i++) out[i] = in[i];
}

void copy_image8(uint8 *in, uint8 *out, const int w, const int h)
{
    int i, size = w*h;
    for(i=0; i < size; i++) out[i] = in[i];
}

void copy_image8_to16(uint8 *in, int16 *out, const int w, const int h)
{
    int i, size = w*h;
    for(i=0; i < size; i++) out[i] = in[i];
}

int main(int argc, const char *argv[]) {

    const char *in_file = NULL;
    const char *out_file = NULL;
    FILE *IN_FILE, *OUT_FILE;
    int i, j, ok = 0, tr = 0, par, fc = 0;
    uint8 *buff = NULL, *buff1 = NULL, *tmpb = NULL;
    int min=0, max=0, size;
    TransState ts[2], *cr;
    void *tmp;

    cr = &ts[0];

    for (i = 1; i < argc; i++)  if (!strcmp(argv[i], "-v")) verb = 1;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
            printf("Usage: itlib -i[input options [bayer_option]] in_file [-t transform_options] [-o[output_options] out_file]\n\n"
                   "Input options:\n"
                   "  b  Input file is raw BAYER image.\n"
                   "  g  Input file is grayscale image.\n"
                   "  r  Input file is RGB image.\n"
                   "\n"
                   "bayer_option:\n"
                   "  bggr            grbg            gbrg            rggb  \n"
                   "    0 1 2 3 4 5	    0 1 2 3 4 5	    0 1 2 3 4 5	    0 1 2 3 4 5 \n"
                   "  0 B G B G B G	  0 G R G R G R	  0 G B G B G B	  0 R G R G R G \n"
                   "  1 G R G R G R	  1 B G B G B G	  1 R G R G R G	  1 G B G B G B \n"
                   "  2 B G B G B G	  2 G R G R G R	  2 G B G B G B	  2 R G R G R G \n"
                   "  3 G R G R G R	  3 B G B G B G	  3 R G R G R G	  3 G B G B G B \n"
                   "\n"
                   "Transform options:\n"
                   "  wb               White balancing.\n"
                   "  bay_rgb          Bayer to rgb bilinear interpolation algorithm.\n"
                   "  bay_grey         Bayer to grey bilinear interpolation algorithm.\n"
                   "  med <a>          3x3 median filter if a = 0 non adaptive, if a = 1 adaptive\n"
                   "  ace              Automatic Color Enhancement transform\n"
                   "  aver <x>         Averaging image with window of x radius \n"
                   "  sub              Subtract one image from another\n"
                   "\n"
                   "Output options:\n"
                   "  -h     ....... this help message.\n"
                   "  -v     ....... verbose \n"
                   );
            return 0;
        } else if (i == 1){
            //Read input file
            if(!strncmp(argv[i], "-i",2)){
                //Find input options
                if(!strncmp(&argv[i][2], "b",1)){
                    cr->colort = BAYER;
                    i++;
                    if(argc > 2) {
                        if     (!strcmp(argv[i], "bggr")) cr->bg = BGGR;
                        else if(!strcmp(argv[i], "grbg")) cr->bg = GRBG;
                        else if(!strcmp(argv[i], "gbrg")) cr->bg = GBRG;
                        else if(!strcmp(argv[i], "rggb")) cr->bg = RGGB;
                        else {
                            fprintf(stderr, "Error! Pls, write bayer grids pattern option\n");
                            ok = 1; goto End;
                        }
                    } else {
                        fprintf(stderr, "Error! Usage: itlib -h\n");
                        ok = 1; goto End;
                    }
                } else if(!strncmp(&argv[i][2], "g",1)){
                    cr->colort = GREY;
                } else if(!strncmp(&argv[i][2], "r",1)){
                    cr->colort = RGB;
                } else {
                    fprintf(stderr, "Error! Can't support '%s' options\n", argv[i]);
                    ok = 1; goto End;
                }

                //if(verb) printf("Color type of input file is %d\n", cr->colort);

                in_file = argv[++i];
                IN_FILE = fopen(in_file, "rb");
                if (IN_FILE == NULL) {
                    fprintf(stderr, "Error! Cannot open input file '%s'\n", in_file);
                    ok = 1; goto End;
                }
                if(!strcmp(&in_file[strlen(in_file)-4],".pgm") || !strcmp(&in_file[strlen(in_file)-4],".PGM")){

                    ok = readPGM(IN_FILE, &buff, &cr->w, &cr->h, &cr->bpp);

                    cr->pic[0] = buff;
                    cr->pic[1] = &buff[cr->w*cr->h*2*3]; //Temporary buffer


                    if(cr->bpp > 8) {
                        utils_cnange_bytes(cr->pic[0], cr->w, cr->h);
                        utils_get_stat(cr->pic[0], cr->w, cr->h, &cr->bpp, &min, &max);
                    } else {
                        copy_image8_to16(cr->pic[0], cr->pic[1], cr->w, cr->h);
                        tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
                    }

                    if(verb) printf("Read %s file w = %d h = %d bpp = %d max = %d min = %d\n", in_file, cr->w, cr->h, cr->bpp, max, min);

                } else if(!strcmp(&in_file[strlen(in_file)-4],".png") || !strcmp(&in_file[strlen(in_file)-4],".PNG")){

                    ok = readPNG(IN_FILE, &buff, &cr->w, &cr->h, &cr->bpp, &cr->colort);

                    cr->pic[0] = buff;
                    cr->pic[1] = &buff[cr->w*cr->h*2*3];
                    copy_image8_to16(cr->pic[0], cr->pic[1], cr->w, cr->h);
                    tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;

                    if(verb) printf("Read %s file w = %d h = %d bpp = %d colort = %d\n", in_file, cr->w, cr->h, cr->bpp, cr->colort);

                } else ok = 1;

                if(IN_FILE) fclose(IN_FILE);
                if(ok){
                    fprintf(stderr, "Error! Read input file '%s'\n", in_file);
                    ok = 1; goto End;
                }
                //Create tmp buffer
                tmpb = (uint8*)malloc(cr->w*cr->h*sizeof(uint32));
                if (tmpb == NULL) {
                    fprintf(stderr, "Error! Can't create tmpb buffer\n");
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! Can't find -i (input file)\n");
                ok = 1; goto End;
            }
        } else if (!strcmp(argv[i], "-o") && i < argc - 1) {
            //Write output file
            out_file = argv[++i];
            OUT_FILE = fopen(out_file, "wb");
            if (OUT_FILE == NULL) {
                fprintf(stderr, "Error! Cannot open input file '%s'\n", out_file);
                ok = 1; goto End;
            }
            if(!strcmp(&out_file[strlen(out_file)-4],".pgm") || !strcmp(&out_file[strlen(out_file)-4],".PGM")){

                ok = writePGM(OUT_FILE, cr->pic[0], cr->w, cr->h, cr->bpp);
                //if(verb) printf("Write %s file\n", out_file);

            } else if(!strcmp(&out_file[strlen(out_file)-4],".png") || !strcmp(&out_file[strlen(out_file)-4],".PNG")){

                if(cr->colort == GREY || cr->colort == BAYER){
                    utils_16_to_8(cr->pic[0], cr->pic[1], cr->w, cr->h, cr->bpp, 1);
                    tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; cr->bpp = 8;

                    ok = writePNG(OUT_FILE, cr->pic[0], cr->w, cr->h, cr->bpp, GREY);
                } else if (cr->colort == RGB){
                    utils_16_to_8(cr->pic[0], cr->pic[1], cr->w, cr->h, cr->bpp, 3);
                    tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; cr->bpp = 8;

                    ok = writePNG(OUT_FILE, cr->pic[0], cr->w, cr->h, cr->bpp, cr->colort);
                } else {
                    fprintf(stderr, "Error! Don't support color type %d\n", cr->colort);
                    ok = 1; goto End;
                }
                //if(verb) printf("Write %s file\n", out_file);

            } else ok = 1;

            if(OUT_FILE) fclose(OUT_FILE);
            if(ok){
                fprintf(stderr, "Error! Write output file '%s'\n", out_file);
                ok = 1; goto End;
            }

        } else if (!strcmp(argv[i], "-t")) {
            //The image transform
            tr = 1;
        } else if (!strcmp(argv[i], "!") && tr) {
            if(fc){
                 fc++;
                 cr = &ts[0];
            } else {
                fc++;
                //Clone image
                ts[1].w = ts[0].w; ts[1].h = ts[0].h; ts[1].bpp = ts[0].bpp; ts[1].bg = ts[0].bg; ts[1].colort = ts[0].colort;
                cr = &ts[1];
                size = ts[0].w*ts[0].h;
                buff1 = (uint8*)malloc(cr->h*cr->h*3*4);
                if (cr->pic[1] == NULL) {
                    fprintf(stderr, "Error! readPNG: Can't allocate memory buff1.\n");
                    ok = 1; goto End;
                }
                cr->pic[0] = buff1;
                cr->pic[1] = &buff1[cr->w*cr->h*2*3]; //Temporary buffer
                copy_image16(ts[0].pic[0], ts[1].pic[0], ts[0].w, ts[0].h);

            }
        } else if (!strcmp(argv[i], "wb") && tr) {
            //White balancing.........................................................................
            if(cr->colort == BAYER ){
                utils_wb_bayer(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, cr->bpp, cr->bg);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
            } else if (cr->colort == RGB){
                utils_wb(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, cr->bpp, cr->bg);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! wb: Input image should be in bayer or rgb format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("white balancing \n");

        } else if (!strcmp(argv[i], "bay_rgb") && tr) {
            if(cr->colort == BAYER){
                utils_bay_to_rgb_bi(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, cr->bg);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; cr->colort = RGB;
            } else {
                fprintf(stderr, "Error! bay_to_rgb_bi: Input image should be in bayer format .\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("bay_to_rgb_bi transform\n");

        }  else if (!strcmp(argv[i], "bay_grey") && tr) {
            if(cr->colort == BAYER){
                utils_bay_to_grey_bi(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, cr->bg);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; cr->colort = GREY;
            } else {
                fprintf(stderr, "Error! bay_to_grey_bi: Input image should be in bayer format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("bay_to_grey_bi transform\n");

        } else if (!strcmp(argv[i], "med") && tr) {
            par = strtol(argv[i+1], NULL, 0);

            if(cr->colort == BAYER){
                filters_median_bayer(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, par);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
            } else if(cr->colort == GREY){
                filters_median(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, par);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! median: Input image should be in bayer format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("median filter\n");
        } else if (!strcmp(argv[i], "ace") && tr) {
            if(cr->colort == BAYER || cr->colort == GREY){
                hdr_ace(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, cr->bpp, 8);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; cr->bpp = 8;
            } else {
                fprintf(stderr, "Error! ace: Input image should be in bayer or grey.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("ace filter\n");
        } else if (!strcmp(argv[i], "ace_local") && tr) {
            if(cr->colort == BAYER || cr->colort == GREY){
                hdr_ace_local(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, cr->bpp);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; cr->bpp = 8;
            } else {
                fprintf(stderr, "Error! ace_local: Input image should be in bayer or grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("ace_local filter\n");
        } else if (!strcmp(argv[i], "aver") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(cr->colort == GREY){
                utils_average(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, par);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
            } else if(cr->colort == BAYER){
                utils_average_bayer(cr->pic[0], cr->pic[1], tmpb, cr->w, cr->h, par);
                tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! average: Input image should be in bayer or grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("average filter\n");
        } else if (!strcmp(argv[i], "sub") && tr ) {
            if(fc == 3){
                if(cr->colort == GREY || cr->colort == BAYER){
                    utils_subtract(cr->pic[0], ts[1].pic[0], cr->pic[1], cr->w, cr->h, cr->bpp);
                    tmp = cr->pic[0]; cr->pic[0] = cr->pic[1]; cr->pic[1] = tmp; //cr->bpp = 8;
                } else {
                    fprintf(stderr, "Error! subtract: Input image should be in bayer or grey format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! subtract: Should be two image in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("subtract two image\n");
        }
    }

    if(argc == 1) printf("Usage: itlib -h\n");
End:
    if(buff) free(buff);
    if(buff1) free(buff1);
    if(tmpb) free(tmpb);
    return ok;
}

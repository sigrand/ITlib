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
    stride = (has_alpha ? 4 : 3) * (*w);

    //Create buffer up to 2 times more the bayer or gray image for the future used
    *buff = (uint8*)malloc(stride*(*h)*4);
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

    //if(verb) printf("readPNG: Finish reading file h = %d\n", *h);
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

void copy_image8_to16(uint8 *in, int16 *out, const int size)
{
    int i;
    for(i=0; i < size; i++) {
        out[i] = in[i];
    }
}

int main(int argc, const char *argv[]) {

    char *in_file[2], cb[2][1000];
    char *out_file = NULL;
    FILE *IN_FILE, *OUT_FILE;
    int i, j, n = 0, ok = 0, tr = 0, par, fc = 0, f, nf, ster = 0;
    uint8 *buff[2], *tmpb = NULL;
    int min=0, max=0, size;
    TransState ts[2];
    void *tmp;
    in_file[0] = NULL; in_file[1] = NULL;
    buff[0] = NULL; buff[1] = NULL;

    //cr[0] = &ts[0];

    for (i = 1; i < argc; i++)  if (!strcmp(argv[i], "-v")) verb = 1;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help")) {
            printf("Usage:\n"
                   "For single image input:\n"
                   "itlib -i[input options] [bayer_option] in_file [-t transform_options] [-o[output_options] out_file]\n"
                   "For stereo pair input:\n"
                   "itlib -is in_dir in_file1 in_file2 [-t transform_options] [-o[output_options] out_file]\n\n"
                   "Input options:\n"
                   "  b  Input file is raw BAYER image.\n"
                   "  s  Input stereo pair.\n"
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
                   "  !  !  !          Image copy construction\n"
                   "  in_image ! in_image transforms ! copy_in_image transforms ! in_image and copy_in_image transforms\n\n"
                   "  wb               White balancing.\n"
                   "  bay_rgb_bi       Bayer to rgb bilinear interpolation algorithm.\n"
                   "  bay_grey_bi      Bayer to grey bilinear interpolation algorithm.\n"
                   "  bay_grey_s       Bayer to rgb bicubic B-spline aproximation \n"
                   "  med <a>          3x3 median filter if a = 0 non adaptive, if a = 1 adaptive\n"
                   "  ace <b>          Automatic Color Enhancement transform b - the bits per pixel for output image\n"
                   "  aver <x>         Averaging image with window of x radius \n"
                   "  sub              Subtract one image from another\n"
                   "  add              Add two images\n"
                   "  dnois_nlm <x>    Non-Local means denoise algorithm, x - radius around pixels for matching\n"
                   "  dnois_reg <x>    The mean square error (MSE) regression of plane denoise filter\n"
                   "                   x - radius around pixels for matching\n"
                   "  hess             Calculate the determinant of Hessian of grey image\n"
                   "  mse <x>          Calculate the MSE, x - radius around pixels for matching \n"
                   "  dnois            Denoise algorithm base on sum of difference\n"
                   "  rgb_yuv444       rgb to yuv444 transform\n"
                   "  rgb_yuv420       rgb to yuv420 transform\n"
                   "  res_down_2       Two times downsampling image\n"
                   "\n"
                   "  grad <x>         The image gradient, x - gardient threshould, if less than th,  = 0 \n"
                   "  lmax             Find local maximums \n"
                   "  end_edge         Find end of edges \n"
                   "  edge             Edge detection \n"
                   "  canny <x>        Canny edge detection, x - gardient threshould, if less than th,  = 0 \n"
                   "  corner <x>       Corners detection, x - threshould, if less than th,  = 0 \n"
                    "\n"
                   "  Stereo options\n"
                   "  s_disp           Calculate disparity"
                   "Output options:\n"
                   "  -h               This help message.\n"
                   "  -v               Verbose \n"
                   );
            return 0;
        } else if (i == 1){
            //Read input file
            if(!strncmp(argv[i], "-i",2)){
                //Find input options
                j = 2; nf = 1;
                if (!strncmp(&argv[i][j], "s",1)){
                    //For two stereo imput images
                    j++;
                    nf = 2;
                    ster = 1;
                }

                if(!strncmp(&argv[i][j], "b",1) ){
                    //Bore bayer input image.
                    for(f=0; f < nf; f++) ts[n].colort = BAYER;
                    i++;
                    if(argc > 2) {
                        if     (!strcmp(argv[i], "bggr")) ts[n].bg = BGGR;
                        else if(!strcmp(argv[i], "grbg")) ts[n].bg = GRBG;
                        else if(!strcmp(argv[i], "gbrg")) ts[n].bg = GBRG;
                        else if(!strcmp(argv[i], "rggb")) ts[n].bg = RGGB;
                        else {
                            fprintf(stderr, "Error! Pls, write bayer grids pattern option\n");
                            ok = 1; goto End;
                        }
                    } else {
                        fprintf(stderr, "Error! Usage: itlib -h\n");
                        ok = 1; goto End;
                    }
                }

                if(nf==1) in_file[0] = (char*)argv[i+1];
                else {
                    strcpy (cb[0], argv[++i]);
                    strcpy (cb[1], argv[i]);
                    strcat (cb[0], argv[++i]);
                    strcat (cb[1], argv[++i]);
                    in_file[0] = cb[0];
                    in_file[1] = cb[1];
                }

                //if(verb) printf("Color type of input file is %d\n", ts[n].colort);
                for(f=0; f < nf; f++) {

                    IN_FILE = fopen(in_file[f], "rb");
                    if (IN_FILE == NULL) {
                        fprintf(stderr, "Error! Cannot open input file '%s'\n", in_file[f]);
                        ok = 1; goto End;
                    }
                    if(!strcmp(&in_file[f][strlen(in_file[f])-4],".pgm") || !strcmp(&in_file[f][strlen(in_file[f])-4],".PGM")){

                        ok = readPGM(IN_FILE, &buff[f], &ts[f].w, &ts[f].h, &ts[f].bpp);

                        ts[f].pic[0] = buff[f];
                        ts[f].pic[1] = &buff[f][ts[f].w*ts[f].h*2*3]; //Temporary buffer


                        if(ts[f].bpp > 8) {
                            utils_cnange_bytes(ts[f].pic[0], ts[f].w, ts[f].h);
                            utils_get_stat(ts[f].pic[0], ts[f].w, ts[f].h, &ts[f].bpp, &min, &max);
                        } else {
                            copy_image8_to16(ts[f].pic[0], ts[f].pic[1], ts[f].w*ts[f].h);
                            tmp = ts[f].pic[0]; ts[f].pic[0] = ts[f].pic[1]; ts[f].pic[1] = tmp;
                        }

                        if(verb) printf("Read %s file w = %d h = %d bpp = %d max = %d min = %d\n", in_file[f], ts[f].w, ts[f].h, ts[f].bpp, max, min);

                    } else if(!strcmp(&in_file[f][strlen(in_file[f])-4],".png") || !strcmp(&in_file[f][strlen(in_file[f])-4],".PNG")){

                        ok = readPNG(IN_FILE, &buff[f], &ts[f].w, &ts[f].h, &ts[f].bpp, &ts[f].colort);

                        ts[f].pic[0] = buff[f];
                        ts[f].pic[1] = &buff[f][ts[f].w*ts[f].h*2*3];
                        if(ts[f].colort == RGB) copy_image8_to16(buff[f], ts[f].pic[1], ts[f].w*ts[f].h*3);
                        else copy_image8_to16(buff[f], ts[f].pic[1], ts[f].w*ts[f].h);
                        tmp = ts[f].pic[0]; ts[f].pic[0] = ts[f].pic[1]; ts[f].pic[1] = tmp; ts[f].bpp = 16;

                        if(verb) printf("Read %s file w = %d h = %d bpp = %d colort = %d\n", in_file[f], ts[f].w, ts[f].h, ts[f].bpp, ts[f].colort);

                    } else ok = 1;

                    if(IN_FILE) fclose(IN_FILE);
                    if(ok){
                        fprintf(stderr, "Error! Read input file '%s'\n", in_file[f]);
                        ok = 1; goto End;
                    }
                }

                //Create tmp buffer
                //printf("f = %d w = %d h = %d\n", f, ts[f].w, ts[f].h);
                tmpb = (uint8*)malloc(ts[0].w*ts[0].h*sizeof(int16)*5);
                if (tmpb == NULL) {
                    fprintf(stderr, "Error! Can't create tmpb buffer\n");
                    ok = 1; goto End;
                }

            }
        }

        if (!strcmp(argv[i], "-o") && i < argc - 1) {
            //Write output file
            out_file = (char*)argv[++i];
            OUT_FILE = fopen(out_file, "wb");
            if (OUT_FILE == NULL) {
                fprintf(stderr, "Error! Cannot open input file '%s'\n", out_file);
                ok = 1; goto End;
            }


            if(!strcmp(&out_file[strlen(out_file)-4],".pgm") || !strcmp(&out_file[strlen(out_file)-4],".PGM")){

                ok = writePGM(OUT_FILE, ts[n].pic[0], ts[n].w, ts[n].h, ts[n].bpp);
                //if(verb) printf("Write %s file\n", out_file);

            } else if(!strcmp(&out_file[strlen(out_file)-4],".png") || !strcmp(&out_file[strlen(out_file)-4],".PNG")){
                printf("n = %d color = %d\n", n, ts[n].colort);
                if(ts[n].colort == GREY || ts[n].colort == BAYER || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                    if(ts[n].bpp != 8) {
                        utils_16_to_8(ts[n].pic[0], ts[n].pic[1], ts[n].w, ts[n].h, ts[n].bpp, 0);
                        tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp; ts[n].bpp = 8;
                    }

                    ok = writePNG(OUT_FILE, ts[n].pic[0], ts[n].w, ts[n].h, ts[n].bpp, GREY);
                } else if (ts[n].colort == RGB){
                    if(ts[n].bpp != 8) {
                        utils_16_to_8(ts[n].pic[0], ts[n].pic[1], ts[n].w, ts[n].h, ts[n].bpp, 2);
                        tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp; ts[n].bpp = 8;
                    }

                    ok = writePNG(OUT_FILE, ts[n].pic[0], ts[n].w, ts[n].h, ts[n].bpp, RGB);
                } else {
                    fprintf(stderr, "Error! Don't support color type %d\n", ts[n].colort);
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
                 if(fc == 2) n = 1;
                 else n = 0;
                 //cr[0] = &ts[0];
            } else {
                fc++;
                //Clone image
                if(buff[1] == NULL){
                    //printf("main buff = %p\n", buff[1]);
                    ts[1].w = ts[n].w; ts[1].h = ts[n].h; ts[1].bpp = ts[n].bpp; ts[1].bg = ts[n].bg; ts[1].colort = ts[n].colort;
                    size = ts[n].w*ts[n].h;
                    buff[1] = (uint8*)malloc(ts[1].w*ts[1].h*3*4);
                    if (ts[1].pic[1] == NULL) {
                        fprintf(stderr, "Error! Can't allocate memory buff[1].\n");
                        ok = 1; goto End;
                    }
                    ts[1].pic[0] = buff[1];
                    ts[1].pic[1] = &buff[1][ts[1].w*ts[1].h*2*3]; //Temporary buffer
                    copy_image16(ts[n].pic[0], ts[1].pic[0], ts[n].w, ts[n].h);
                }
            }
            printf("%d ", n);
        } else if (!strcmp(argv[i], "wb") && tr) {
            //White balancing.........................................................................
            if(ts[n].colort == BAYER ){
                utils_wb_bayer(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, ts[n].bpp, ts[n].bg);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else if (ts[n].colort == RGB){
                utils_wb_rgb(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, ts[n].bpp);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! wb: Input image should be in bayer or rgb format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("white balancing \n");

        } else if (!strcmp(argv[i], "bay_rgb_bi") && tr) {
            if(ts[n].colort == BAYER){
                trans_bay_to_rgb_bi(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, ts[n].bg);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp; ts[n].colort = RGB;
            } else {
                fprintf(stderr, "Error! bay_to_rgb_bi: Input image should be in bayer format .\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("bay_to_rgb_bi transform\n");

        }  else if (!strcmp(argv[i], "bay_grey_bi") && tr) {
            if(ts[n].colort == BAYER){
                trans_bay_to_grey_bi(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, ts[n].bg);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp; ts[n].colort = GREY;
            } else {
                fprintf(stderr, "Error! bay_to_grey_bi: Input image should be in bayer format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("bay_to_grey_bi transform\n");

        } else if (!strcmp(argv[i], "med") && tr) {
            par = strtol(argv[i+1], NULL, 0);

            if(ts[n].colort == BAYER){
                filters_median_bayer(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else if(ts[n].colort == GREY){
                filters_median(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! median: Input image should be in bayer format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("median filter\n");
        } else if (!strcmp(argv[i], "ace") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(!par) par = 8;
            //printf("par = %d\n", par);
            if(ts[n].colort < RGBA){
                hdr_ace(ts[n].pic[0], ts[n].pic[1], (int*)tmpb, ts[n].w, ts[n].h, ts[n].bpp, par, ts[n].colort);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp; //ts[n].bpp = par;
            } else {
                fprintf(stderr, "Error! ace: Input image should be in bayer or grey.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("ace filter\n");
        } else if (!strcmp(argv[i], "ace_local") && tr) {
            if(ts[n].colort == BAYER || ts[n].colort == GREY){
                hdr_ace_local(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, ts[n].bpp);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1];  ts[n].pic[1] = tmp; //ts[n].bpp = 8;
            } else {
                fprintf(stderr, "Error! ace_local: Input image should be in bayer or grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("ace_local filter\n");
        } else if (!strcmp(argv[i], "aver") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                utils_average(ts[n].pic[0], ts[n].pic[1], (uint32*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;

            } else if(ts[n].colort == BAYER){
                utils_average_bayer(ts[n].pic[0], ts[n].pic[1], (uint32*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! average: Input image should be in bayer or grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("average filter\n");
        }  else if (!strcmp(argv[i], "hess") && tr) {
            if(ts[n].colort == GREY || ts[n].colort == BAYER || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                filters_hessian(ts[n].pic[0], ts[n].pic[1], (uint32*)tmpb, ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! hess: Input image should be in bayer or grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("hessian filter\n");
        } else if (!strcmp(argv[i], "sub") && tr ) {
            if(fc == 3 || ster){
                fc = 0;
                if(ts[n].colort == GREY || ts[n].colort == BAYER){
                    //hdr_diff(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], ts[n].w, ts[n].h, ts[n].bpp);
                    utils_subtract(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], ts[n].w, ts[n].h, ts[n].bpp);
                    tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;

                } else {
                    fprintf(stderr, "Error! subtract: Input image should be in bayer or grey format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! subtract: Should be two image in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("subtract two image\n");
        } else if (!strcmp(argv[i], "add") && tr ) {
            if(fc == 3 || ster){
                fc = 0;
                if(ts[n].colort == GREY || ts[n].colort == BAYER || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                    utils_add(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], ts[n].w, ts[n].h, ts[n].bpp, ts[1].bpp);
                    tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1];  ts[n].bpp = ts[1].bpp; ts[n].pic[1] = tmp;

                } else {
                    fprintf(stderr, "Error! add: Input image should be in bayer or grey format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! add: Should be two image in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("add two image\n");
        } else if (!strcmp(argv[i], "dnois_nlm") && tr) {
            if(fc == 3 || ster){
                fc = 0;
                par = strtol(argv[i+1], NULL, 0);
                if(ts[n].colort == GREY || ts[n].colort == BAYER){
                    filters_NLM_denoise_bayer(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], (int16*)tmpb, par, 50, ts[n].w, ts[n].h);
                    tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1];  ts[n].bpp = ts[1].bpp; ts[n].pic[1] = tmp;

                } else {
                    fprintf(stderr, "Error! dnois_nlm: Input image should be in bayer or grey format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! dnois_nlm: Should be two image in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("denoise filter\n");
        } else if (!strcmp(argv[i], "mse") && tr) {
            if(fc == 3 || ster){
                fc = 0;
                par = strtol(argv[i+1], NULL, 0);
                if(ts[n].colort == GREY || ts[n].colort == BAYER || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                    filters_MSE_bayer(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], (int16*)tmpb, par, ts[n].bpp, ts[n].w, ts[n].h);
                    tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].bpp = ts[1].bpp; ts[n].pic[1] = tmp;

                } else {
                    fprintf(stderr, "Error! mse: Input image should be in bayer or grey format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! mse: Should be two image in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("mse\n");
        } else if (!strcmp(argv[i], "dnois") && tr) {
            if(fc == 3 || ster){
                fc = 0;
                //par = strtol(argv[i+1], NULL, 0);
                if(ts[n].colort == GREY || ts[n].colort == BAYER){
                    filters_denoise(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], (int*)tmpb, ts[n].bpp, ts[n].w, ts[n].h);
                    tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].bpp = ts[1].bpp; ts[n].pic[1] = tmp;

                } else {
                    fprintf(stderr, "Error! dnois: Input image should be in bayer or grey format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! dnois: Should be two image in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("dnois\n");
        } else if (!strcmp(argv[i], "dnois_reg") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == BAYER){
                filters_denoise_regression_bayer(ts[n].pic[0], ts[n].pic[1], (int*)tmpb, par, ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! dnois_reg: Input image should be in bayer or grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("denoise  filter\n");
        } else if (!strcmp(argv[i], "bay_rgb_s") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == BAYER){
                trans_bay_to_rgb_b_spline(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, ts[n].bg);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].colort = RGB; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! bay_rgb_s: Input image should be in bayer format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("bay_rgb_s spline approximation \n");

        } else if (!strcmp(argv[i], "rgb_yuv444") && tr) {
            if(ts[n].colort == RGB){
                trans_rgb_to_yuv444(ts[n].pic[0], ts[n].pic[1], &((int16*)ts[n].pic[1])[ts[n].w*ts[n].h], &((int16*)ts[n].pic[1])[ts[n].w*ts[n].h<<1], ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].colort = YUV444; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! rgb_yuv444: Input image should be in rgb format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("rgb_yuv444 transform \n");
        } else if (!strcmp(argv[i], "rgb_yuv420") && tr) {
            if(ts[n].colort == RGB){
                trans_rgb_to_yuv420(ts[n].pic[0], ts[n].pic[1], &((int16*)ts[n].pic[1])[ts[n].w*ts[n].h], &((int16*)ts[n].pic[1])[ts[n].w*ts[n].h<<1], ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].colort = YUV420; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! rgb_yuv420: Input image should be in rgb format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("rgb_yuv420 transform \n");
        } else if (!strcmp(argv[i], "grad") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                seg_gradient(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;

            } else {
                fprintf(stderr, "Error! grad: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("gradient transform \n");
        } else if (!strcmp(argv[i], "lmax") && tr) {
            //par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                seg_local_max(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! lmax: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Find local maximum\n");
        } else if (!strcmp(argv[i], "edge") && tr) {
            //par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                seg_edge_detection(ts[n].pic[0], ts[n].pic[1], ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;  ts[n].bpp = 8;
            } else {
                fprintf(stderr, "Error! edge: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Edge detection\n");
        } else if (!strcmp(argv[i], "canny") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                seg_canny_edge(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! canny: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Canny edge detection\n");
        } else if (!strcmp(argv[i], "res_down_2") && tr) {
            //par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                utils_resize_down_2(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
                ts[n].w = ts[n].w>>1; ts[n].h = ts[n].h>>1;
            } else {
                fprintf(stderr, "Error! res_down_2: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Two times downsampling image\n");
        } else if (!strcmp(argv[i], "end_edge") && tr) {
            //par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                seg_end_of_edges(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! end_of_edges: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Found end of edges\n");
        } else if (!strcmp(argv[i], "corner") && tr) {
            par = strtol(argv[i+1], NULL, 0);
            if(ts[n].colort == GREY || ts[n].colort == YUV444 || ts[n].colort == YUV420){
                seg_corners(ts[n].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h, par);
                tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp;
            } else {
                fprintf(stderr, "Error! corners: Input image should be in grey format.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Corners detection\n");
        } else if (!strcmp(argv[i], "s_disp") && tr ) {
            if(ster){
                n = 0;
                if(ts[n].colort == RGB){
                    stereo_disparity(ts[n].pic[0], ts[1].pic[0], ts[n].pic[1], (int16*)tmpb, ts[n].w, ts[n].h);
                    tmp = ts[n].pic[0]; ts[n].pic[0] = ts[n].pic[1]; ts[n].pic[1] = tmp; ts[n].colort = GREY;

                } else {
                    fprintf(stderr, "Error! s_disp: Input image should be RGB format.\n", out_file);
                    ok = 1; goto End;
                }
            } else {
                fprintf(stderr, "Error! s_disp: Should be stereo pair in input.\n", out_file);
                ok = 1; goto End;
            }
            if(verb) printf("Calculate disparity\n");
        }
    }

    if(argc == 1) printf("Usage: itlib -h\n");
End:
    if(buff[0]) free(buff[0]);
    if(buff[1]) free(buff[1]);
    if(tmpb) free(tmpb);
    return ok;
}

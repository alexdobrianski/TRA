/***********************************************************************
     
    2009-2011 (C) Alex Dobrianski tra.cpp oprbit calculation app

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
************************************************************************/
// file modified from sample source code from Jpeg library:
/*
 * jpeglib.h
 *
 * Copyright (C) 1991-1995, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file defines the application interface for the JPEG library.
 * Most applications using the library need only include this file,
 * and perhaps jerror.h if they want to know the exact error codes.
 */



#include "stdafx.h"
#define _USE_MATH_DEFINES 1
#include <math.h>
#include <stdio.h>
extern "C"
{
#include "JPEGLIB.H"
}
void
write_JPEG_file (char * filename, 
				 int quality, // 
				 int SizeW, // Width
				 int SizeH, // hight
				 int SizeB, // color components
				 unsigned char *bArray, 
				 J_COLOR_SPACE ColorCode
				 )
{
  /* This struct contains the JPEG compression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   * It is possible to have several such structures, representing multiple
   * compression/decompression processes, in existence at once.  We refer
   * to any one struct (and its associated working data) as a "JPEG object".
   */
  struct jpeg_compress_struct cinfo;
  struct jpeg_compress_struct *pcinfo;
  /* This struct represents a JPEG error handler.  It is declared separately
   * because applications often want to supply a specialized error handler
   * (see the second half of this file for an example).  But here we just
   * take the easy way out and use the standard error handler, which will
   * print a message on stderr and call exit() if compression fails.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct jpeg_error_mgr jerr;
  struct jpeg_error_mgr *pjerr;
  /* More stuff */
  FILE * outfile;		/* target file */
  JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
  int row_stride;		/* physical row width in image buffer */
  int iRealH;
  JSAMPLE * im_buffer = (JSAMPLE*)bArray;
  int CountOfSets = 0;

  /* Step 1: allocate and initialize JPEG compression object */

  /* We have to set up the error handler first, in case the initialization
   * step fails.  (Unlikely, but it could happen if you are out of memory.)
   * This routine fills in the contents of struct jerr, and returns jerr's
   * address which we place into the link field in cinfo.
   */
  cinfo.err = jpeg_std_error(&jerr);

  /* Now we can initialize the JPEG compression object. */

  /* Step 2: specify data destination (eg, a file) */
  /* Note: steps 2 and 3 can be done in either order. */

  /* Here we use the library-supplied code to send compressed data to a
   * stdio stream.  You can also write your own code to do something else.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to write binary files.
   */
   if (((outfile) = fopen(filename, "wb")) == NULL) 
   {
	return ;
   }
    
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);


  /* Step 3: set parameters for compression */

  /* First we supply a description of the input image.
   * Four fields of the cinfo struct must be filled in:
   */
  cinfo.image_width = SizeW; 	/* image width and height, in pixels */
  cinfo.image_height = SizeH;
  cinfo.input_components = SizeB;		/* # of color components per pixel */
  cinfo.in_color_space = ColorCode;//JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_FASTEST;//JDCT_FLOAT;

//  pcinfo->image_width = SizeW; 	/* image width and height, in pixels */
//  pcinfo->image_height = SizeH;
//  pcinfo->input_components = SizeB;		/* # of color components per pixel */
//  pcinfo->in_color_space = ColorCode;//JCS_RGB; 	/* colorspace of input image */
//  pcinfo->dct_method = JDCT_FASTEST;//JDCT_FLOAT;

  
  //cinfo.optimize_coding = TRUE;
  //cinfo.restart_in_rows = 240;
		cinfo.dct_method = JDCT_FLOAT;//JDCT_FASTEST;//JDCT_FLOAT;
		cinfo.optimize_coding = TRUE;
			//cinfo.smoothing_factor = 100;
			cinfo.smoothing_factor = 0;
  /* Now use the library's routine to set default compression parameters.
   * (You must set at least cinfo.in_color_space before calling this,
   * since the defaults depend on the source color space.)
   */
  jpeg_set_defaults(&cinfo);
  //jpeg_set_defaults(pcinfo);
  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);
       //jpeg_set_quality(pcinfo, quality, TRUE/*FALSE*/ );
	//	jpeg_simple_progression (pcinfo);
	//  jpeg_set_quality(pcinfo, quality, TRUE);


  /* Step 4: Start compressor */

  /* TRUE ensures that we will write a complete interchange-JPEG file.
   * Pass TRUE unless you are very sure of what you're doing.
   */
  jpeg_start_compress(&cinfo, TRUE);
  //jpeg_start_compress(pcinfo, TRUE);

  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */

  /* Here we use the library's state variable cinfo.next_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * To keep things simple, we pass one scanline per call; you can pass
   * more if you wish, though.
   */
  row_stride = SizeW * SizeB;	/* JSAMPLEs per row in image_buffer */
 // *prow_stride = SizeW * SizeB;	/* JSAMPLEs per row in image_buffer */

	
//  while (cinfo.next_scanline < cinfo.image_height) {
//    row_pointer[0] = & im_buffer[cinfo.next_scanline * row_stride];
//    (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
//  }
		  iRealH = cinfo.image_height;



	  while (cinfo.next_scanline < iRealH) 
	  {
		/* jpeg_write_scanlines expects an array of pointers to scanlines.
		* Here the array is only one element long, but you could pass
		* more than one scanline at a time if that's more convenient.
		*/
		row_pointer[0] = & im_buffer[cinfo.next_scanline * (row_stride)];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	  }


  /* Step 6: Finish compression */

  jpeg_finish_compress(&cinfo);
  //jpeg_finish_compress(pcinfo);
  /* After finish_compress, we can close the output file. */
  if (filename)
  {
	fclose(outfile);
    //fclose(*poutfile);
  }

  /* Step 7: release JPEG compression object */

  /* This is an important step since it will release a good deal of memory. */
  jpeg_destroy_compress(&cinfo);
  //jpeg_destroy_compress(pcinfo);

  /* And we're done! */
}
///////////////////////////////////////////////////


/******************************************/
/*
  reverb.cpp
  by Gary P. Scavone, 2006-2007.
  by Robin Larvor and Inoë ANDRE 2016-2017

  This program opens a duplex stream and passes
  input and put a reverb effect directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "include/somefunc.h"



/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8
*/

typedef signed short MY_TYPE;
#define FORMAT RTAUDIO_SINT16

/*
typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64
*/

void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device to use (default = 0),\n";
  std::cout << "    oDevice = optional output device to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
  exit( 0 );
}


int conv_add(double* h, double* x,double* prec, unsigned int L, long M)
{
  int i=0,j=0;
  int tmp=0;
  int kmin=0;
  int kmax=0;
  double * conv= (double *)malloc(sizeof(conv)*(L+M-1));

  for(i=0;i<L+M-1;i++)
    {
      tmp=0;
      if(i>=M){
	kmin=i-M+1;
      } else {
	kmin= 1;
      }

      if(i<L){
	kmax=i;
      } else {
	kmax = L;
      }


      for(j=kmin;j<=kmax;j++){
	tmp+=x[j]*h[i-j+1];
      }
      conv[i]=tmp;
    }
  
  //bloc L
  for(i=0;i<L;i++){
    x[i]=conv[i];
  }
  //bloc M
  for(i=0;i<M-1;i++){
    prec[i]=conv[i+L+1];
  }
  
  return 0;
}


int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data_stream )
{
  
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;
  conv_add( ((data)data_stream)->h, (double *)inputBuffer,((data)data_stream)->buffer_prec, ((data)data_stream)->L,((data)data_stream)->M);
  unsigned int bytes = ((data)data_stream)->L ;
  memcpy( outputBuffer, inputBuffer, bytes );
  return 0;
}


int main( int argc, char *argv[] )
{

//read the binary file
  FILE * pFile;
  long lSize;
  double * h_filter;
  size_t nb_buffer;

 

  //ouverture du fichier contenant de la réponse impulsionnelle du filtre
  pFile = fopen ( "../../ressources_tstr_v1_1/c/impres" , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  h_filter = (double*) malloc (sizeof(char)*lSize);
  if (h_filter == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  nb_buffer = fread (h_filter,sizeof(double),lSize,pFile);
  if ((long)nb_buffer*sizeof(double) != lSize) {printf("nb_buffer:%u lSize:%u \n", (long)nb_buffer, lSize); fputs ("Reading error\n ",stderr); exit (3);}
/* the whole file is now loaded in the memory buffer. */

  // int i = 0;
  // for(i=1000;i<1010;i++){
  // std::cout << "h_filter :" << h_filter[i] << std::endl;
  // }




  unsigned int channels, fs, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

  // Minimal command-line checking
  if (argc < 3 || argc > 7 ) usage();

  RtAudio adac;
  if ( adac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  if ( argc > 3 )
    iDevice = (unsigned int) atoi(argv[3]);
  if ( argc > 4 )
    oDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    iOffset = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    oOffset = (unsigned int) atoi(argv[6]);

  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 512;
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = iDevice;
  iParams.nChannels = channels;
  iParams.firstChannel = iOffset;
  oParams.deviceId = oDevice;
  oParams.nChannels = channels;
  oParams.firstChannel = oOffset;

  if ( iDevice == 0 )
    iParams.deviceId = adac.getDefaultInputDevice();
  if ( oDevice == 0 )
    oParams.deviceId = adac.getDefaultOutputDevice();

  RtAudio::StreamOptions options;
  //options.flags |= RTAUDIO_NONINTERLEAVED;



   //allocation de la structure pour la stream
  
  data data_stream = (data)malloc(sizeof(*data_stream));
  
  data_stream->h = (double*)calloc(nb_buffer,sizeof(*data_stream->h));
  data_stream->fft_h = (double*)calloc(nb_buffer,sizeof(*data_stream->fft_h));
  data_stream->M= nb_buffer;
  data_stream->L= 512;//bufferFrames * channels * sizeof( MY_TYPE );
  data_stream->buffer_prec = (double*)calloc( (data_stream->M + data_stream->L - 1) , sizeof(*data_stream->buffer_prec));
  
  //remplir la structure
  data_stream->h=h_filter;


					     
  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)data_stream, &options );
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    exit( 1 );
  }

  data_stream->L = 512;//bufferFrames * channels * sizeof( MY_TYPE );

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  try {
    adac.startStream();

    char input;
    std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
    std::cin.get(input);

    // Stop the stream.
    adac.stopStream();
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    goto cleanup;
  }

 cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  // terminate
  fclose (pFile);
  free (h_filter);

  //faire des free dans data_stream

  free(data_stream);
  return 0;
}

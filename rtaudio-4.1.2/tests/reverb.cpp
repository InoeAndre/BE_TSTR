/******************************************/
/*
  reverb.cpp
  by Gary P. Scavone, 2006-2007.
  by Robin Larvor and Inoë ANDRE 2016-2017

  This program opens a duplex stream and passes
  input and is supposed to put a reverb effect directly through to the output either in FREQUENCY or TIME DOMAIN.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include "somefunc.h"


#define SHORT_FILTER //mettre en commentaire pour avoir la réponse impulsionnelle complète
#define DOMAINE_FREQUENTIEL //mettre en commentaire pour avoir la convolution add dans le domaine temporelle


/*
  typedef char MY_TYPE;
  #define FORMAT RTAUDIO_SINT8


  typedef signed short MY_TYPE;
  #define FORMAT RTAUDIO_SINT16


  typedef S24 MY_TYPE;
  #define FORMAT RTAUDIO_SINT24

  typedef signed long MY_TYPE;
  #define FORMAT RTAUDIO_SINT32

  typedef float MY_TYPE;
  #define FORMAT RTAUDIO_FLOAT32
*/
typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64


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

/**************************************************************************************
                 fonction copier du fichier somefunc.cpp
***************************************************************************************/

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


double *_sintbl = 0;
int maxfftsize = 0;


///////////////////////////////
// FFT functions
size_t get_nextpow2(size_t n)
{
  size_t k = 1;
  while (k < n){
    k *= 2;
  }

  return k;
}


static int checkm(const int m)
{
  int k;

  for (k = 4; k <= m; k <<= 1) {
    if (k == m)
      return (0);
  }
  fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n",m);

  return (-1);
}

int fft(double *x, double *y, const int m)
{
  int j, lmx, li;
  double *xp, *yp;
  double *sinp, *cosp;
  int lf, lix, tblsize;
  int mv2, mm1;
  double t1, t2;
  double arg;
  int checkm(const int);

  /**************
   * RADIX-2 FFT *
   **************/

  if (checkm(m))
    return (-1);

  /***********************
   * SIN table generation *
   ***********************/

  if ((_sintbl == 0) || (maxfftsize < m)) {
    tblsize = m - m / 4 + 1;
    arg = M_PI / m * 2;
    if (_sintbl != 0)
      free(_sintbl);
    _sintbl = sinp = dgetmem(tblsize);
    *sinp++ = 0;
    for (j = 1; j < tblsize; j++)
      *sinp++ = sin(arg * (double) j);
    _sintbl[m / 2] = 0;
    maxfftsize = m;
  }

  lf = maxfftsize / m;
  lmx = m;

  for (;;) {
    lix = lmx;
    lmx /= 2;
    if (lmx <= 1)
      break;
    sinp = _sintbl;
    cosp = _sintbl + maxfftsize / 4;
    for (j = 0; j < lmx; j++) {
      xp = &x[j];
      yp = &y[j];
      for (li = lix; li <= m; li += lix) {
	t1 = *(xp) - *(xp + lmx);
	t2 = *(yp) - *(yp + lmx);
	*(xp) += *(xp + lmx);
	*(yp) += *(yp + lmx);
	*(xp + lmx) = *cosp * t1 + *sinp * t2;
	*(yp + lmx) = *cosp * t2 - *sinp * t1;
	xp += lix;
	yp += lix;
      }
      sinp += lf;
      cosp += lf;
    }
    lf += lf;
  }

  xp = x;
  yp = y;
  for (li = m / 2; li--; xp += 2, yp += 2) {
    t1 = *(xp) - *(xp + 1);
    t2 = *(yp) - *(yp + 1);
    *(xp) += *(xp + 1);
    *(yp) += *(yp + 1);
    *(xp + 1) = t1;
    *(yp + 1) = t2;
  }

  /***************
   * bit reversal *
   ***************/
  j = 0;
  xp = x;
  yp = y;
  mv2 = m / 2;
  mm1 = m - 1;
  for (lmx = 0; lmx < mm1; lmx++) {
    if ((li = lmx - j) < 0) {
      t1 = *(xp);
      t2 = *(yp);
      *(xp) = *(xp + li);
      *(yp) = *(yp + li);
      *(xp + li) = t1;
      *(yp + li) = t2;
    }
    li = mv2;
    while (li <= j) {
      j -= li;
      li /= 2;
    }
    j += li;
    xp = x + j;
    yp = y + j;
  }

  return (0);
}

int fftr(double *x, double *y, const int m)
{
  int i, j;
  double *xp, *yp, *xq;
  double *yq;
  int mv2, n, tblsize;
  double xt, yt, *sinp, *cosp;
  double arg;

  mv2 = m / 2;

  /* separate even and odd  */
  xq = xp = x;
  yp = y;
  for (i = mv2; --i >= 0;) {
    *xp++ = *xq++;
    *yp++ = *xq++;
  }

  if (fft(x, y, mv2) == -1)    /* m / 2 point fft */
    return (-1);


  /***********************
   * SIN table generation *
   ***********************/

  if ((_sintbl == 0) || (maxfftsize < m)) {
    tblsize = m - m / 4 + 1;
    arg = M_PI / m * 2;
    if (_sintbl != 0)
      free(_sintbl);
    _sintbl = sinp = dgetmem(tblsize);
    *sinp++ = 0;
    for (j = 1; j < tblsize; j++)
      *sinp++ = sin(arg * (double) j);
    _sintbl[m / 2] = 0;
    maxfftsize = m;
  }
  //printf("Debug: m=%i, maxfftsize=%i\n",m,maxfftsize);

  n = maxfftsize / m;
  sinp = _sintbl;
  cosp = _sintbl + maxfftsize / 4;

  xp = x;
  yp = y;
  xq = xp + m;
  yq = yp + m;
  *(xp + mv2) = *xp - *yp;
  *xp = *xp + *yp;
  *(yp + mv2) = *yp = 0;

  for (i = mv2, j = mv2 - 2; --i; j -= 2) {
    ++xp;
    ++yp;
    sinp += n;
    cosp += n;
    yt = *yp + *(yp + j);
    xt = *xp - *(xp + j);
    *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
    *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
  }

  xp = x + 1;
  yp = y + 1;
  xq = x + m;
  yq = y + m;

  for (i = mv2; --i;) {
    *xp++ = *(--xq);
    *yp++ = -(*(--yq));
  }

  return (0);
}

int ifft(double *x, double *y, const int m)
{
  int i;

  if (fft(y, x, m) == -1)
    return (-1);

  for (i = m; --i >= 0; ++x, ++y) {
    *x /= m;
    *y /= m;
  }

  return (0);
}

double *dgetmem(int leng)
{
  return ( (double *)getmem(leng, sizeof(double)) );
}

char *getmem(int leng, unsigned size)
{
  char *p = NULL;

  if ((p = (char *)calloc(leng, size)) == NULL){
    fprintf(stderr, "Memory allocation error !\n");
    exit(3);
  }
  return (p);
}

double get_process_time() {
  struct rusage usage;
  if( 0 == getrusage(RUSAGE_SELF, &usage) ) {
    return (double)(usage.ru_utime.tv_sec + usage.ru_stime.tv_sec) +
      (double)(usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1.0e6;
  }
  return 0;
}




/**************************************************************************************
                 fonctions de convolution add
***************************************************************************************/



/**
 * \fn int conv_add(double* h, double* x,double* conv, double* prec, unsigned int L, long int M)
 *  \brief fonction qui permet de faire une convolution add dans le domaine temporelle
 *  \param h    = filtre impulsionnel, 
 *  \param x    = signal à traiter et à mettre dans le buffer
 *  \param conv = tableau de temporaire pour calculer les nouvelles valeurs du buffer 
 *  \param prec = la portion qui ne peut pas rentrer dans le buffer du précédent calcul
 *  \param L    = longueur du buffer
 *  \param M    = longueur du filtre impulsionnel
 *  \return 0 quand le calcul est fini 
 */
int conv_add(double* h, double* x,double* conv, double* prec, unsigned int L, long int M)
{
  long int i=0;
  int j=0;
  double tmp=0.0;// buffer du résultat d'un indice de la convolution
  long int kmin=0;
  long int kmax=0;

  //calcul de la convolution    				
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

      for(j=kmin-1;j<=kmax-1;j++){
	tmp+=x[j]*h[i-j+1];
  					
      }
      conv[i]=tmp;
    }


  //bloc L : Partie envoyé dans le stream
  for(i=0;i<L;i++){
    x[i]=conv[i]+prec[i];
    //printf("x[%d]= %lf",i,x[i]);
  }

  //bloc M : Le reste est a stocker pour la prochaine convolution. Sachant que M-1 est plus grand que L, il y a encore M-1-L valeur à stocker de l'ancien prec dans le nouveau prec.
  for(i=0;i<M-1;i++){
    if(i+L<M-1)  prec[i]=conv[i+L]+prec[i+L];
    else prec[i]=conv[i+L];
  }
  					
  return 0;
}


/** \fn int conv_add( double* x,double* fft_x,double * tmp_x,double* fft_h, double* prec,double * conv, double* res_ifft, double* tmp3,unsigned int L, long int M,long int FFTOrder)
 *  \brief fonction qui permet de faire une convolution add dans le domaine frequentiel
 *  \param x        = Signal à traiter et à mettre dans le buffer
 *  \param fft_x    = TF du signal du stream  
 *  \param tmp_x    = Pointeur entre la partie imaginaire et la partie réelle pour la fft du signal x
 *  \param fft_h    = Filtre impulsionnel dans le domaine frequentielle,
 *  \param prec     = La portion qui ne peut pas rentrer dans le buffer du précédent calcul
 *  \param conv     = Multiplication des fft
 *  \param res_ifft = TF inverse de conv
 *  \param tmp3     = Pointeur entre la partie imaginaire et la partie réelle pour la fft inverse
 *  \param L        = Longueur du buffer
 *  \param M        = Longueur du filtre impulsionnel
 *  \param FFTOrder = Longueur des FFT : puissance de 2 supérieur la plus proche de L+M-1
 *  \return 0 quand le calcul est fini 
 */
int conv_add_freq( double* x,double* fft_x,double * tmp_x,double* fft_h, double* prec,double * conv, double* res_ifft, double* tmp3,unsigned int L, long int M,long int FFTOrder)
{
  long int i=0;

  //réinitialisation des buffers 
  for(i=0;i<2*FFTOrder;i++){
    fft_x[i]=0.0;
    conv[i]=0.0;
    res_ifft[i]=0.0;
  }

  //Calcul de la FFT du signal du stream
  memcpy(fft_x,x, L);
  fftr(fft_x,tmp_x,FFTOrder); // res_fft contient maintenant le résultat de la FFT

  //multiplication des fft 
  for(i=0;i<FFTOrder;i++){
    //parite réelle
    conv[i]=fft_x[i]*fft_h[i]-fft_x[FFTOrder+i]*fft_h[FFTOrder+i];
    //partie imaginaire
    conv[FFTOrder+i]=fft_x[i]*fft_h[FFTOrder+i]+fft_x[FFTOrder+i]*fft_h[i];
  }
  
  //Calcul de la FFT inverse
  memcpy(res_ifft,conv, 2*FFTOrder*sizeof(double));  
  ifft(res_ifft,tmp3,FFTOrder); // res_ifft contient maintenant la transformée de Fourier inverse

  //bloc L : Partie envoyé dans le stream
  for(i=0;i<L;i++){
    x[i]=res_ifft[i]+prec[i];
    // DEBUG
    // printf("x[%zu]= %lf\nres_ifft[%zu]= %lf\nprec[%zu]= %lf\n",i,x[i],i,res_ifft[i],i,prec[i]);
  }

  //bloc M : Le reste est a stocker pour la prochaine convolution. Sachant que M-1 est plus grand que L, il y a encore M-1-L valeur à stocker de l'ancien prec dans le nouveau prec.
  for(i=0;i<M-1;i++){
    if(i+L<M-1) prec[i]=res_ifft[i+L]+prec[i+L];
    else prec[i]=res_ifft[i+L];
  }

  return 0;
}

int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void *data_stream )
{

  //MESURE DU TEMPS
  double t1 = get_process_time();
  
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

#ifndef DOMAINE_FREQUENTIEL
  conv_add( ((data)data_stream)->h,
	    (double *)inputBuffer,
	    ((data)data_stream)->conv,
	    ((data)data_stream)->buffer_prec,
	    ((data)data_stream)->L,
	    ((data)data_stream)->M);
#else
  conv_add_freq((double *)inputBuffer,
		((data)data_stream)->fft_x,
		((data)data_stream)->tmp_x,
		((data)data_stream)->fft_h,
		((data)data_stream)->buffer_prec,
		((data)data_stream)->conv,
		((data)data_stream)->res_ifft,
		((data)data_stream)->tmp3,
		((data)data_stream)->L,
		((data)data_stream)->M,
		((data)data_stream)->dFFTOrder);
#endif
  
  unsigned int bytes = ((data)data_stream)->L ;
  memcpy( outputBuffer, inputBuffer, sizeof(double)*bytes );

  //MESURE DU TEMPS DIFFERENTIEL
  double t2 = get_process_time();
  printf("time = %lf s\n",t2-t1);

  return 0;
}






int main( int argc, char *argv[] )
{

  //read the binary file
  FILE * pFile;
  long int lSize;
  double * h_filter;
  long int nb_buffer;
 

  //ouverture du fichier contenant de la réponse impulsionnelle du filtre
  pFile = fopen ( "../../ressources_tstr_v1_1/c/impres" , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  h_filter = (double*) malloc (sizeof(char)*lSize);
  if (h_filter == NULL) {
    fputs ("Memory error",stderr);
    exit (2);
  }

  // copy the file into the buffer:
  nb_buffer = fread (h_filter,sizeof(double),lSize,pFile);
  if ((long int)(nb_buffer*sizeof(double)) != lSize) {
    printf("nb_buffer:%zu lSize:%zu \n", nb_buffer, lSize);
    fputs ("Reading error\n ",stderr);
    exit (3);
  }
  /* the whole file is now loaded in the memory buffer. */

  // //DEBUG POUR VOIR SI LE FILTRE LU EST COHERENT
  // int i =0;
  // for(i=100;i<20000;i++){
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


  /**************************************************************************************
                         allocation de la structure pour la stream 
  ****************************************************************************************/

  data data_stream = (data)malloc(sizeof(*data_stream));


  
#ifdef SHORT_FILTER
  //taille racourcie de la fenetre
  size_t nb_buffer_short=40000;
  //nb_buffer = taille de la réponse impulsionnelle
  nb_buffer = nb_buffer_short;
  double * h_short = (double*)calloc(nb_buffer,sizeof(double));
  memcpy(h_short,h_filter,nb_buffer);

  //remplissage de la structure
  (data_stream->h)=h_short;
#else
  (data_stream->h)=h_filter;
#endif

  //suite remplissage de la structure
  data_stream->M= nb_buffer;
  data_stream->L= 512;
  data_stream->buffer_prec = (double*)calloc( (data_stream->M - 1) , sizeof(*data_stream->buffer_prec));
  
#ifndef DOMAINE_FREQUENTIEL
  
  data_stream->conv = (double*)calloc( (data_stream->L + data_stream->M - 1),sizeof(*data_stream->conv));
  
#else
  
  //ordre des fft
  long int FFTOrder = get_nextpow2(data_stream->L+data_stream->M-1);
  data_stream->dFFTOrder=FFTOrder;
  
  //buffer de la multiplication des FFT
  data_stream->conv = (double*)calloc( 2*FFTOrder,sizeof(*data_stream->conv));

  //calcul et stockage de la FFT de h la RI 
  double* res_fft=(double *)calloc(2*FFTOrder,sizeof(double));
  memcpy(res_fft,data_stream->h,data_stream->M);  
  double* tmp = (res_fft + FFTOrder);
  fftr(res_fft,tmp,FFTOrder); 
  //res_fft contient maintenant le résultat de la FFT de h

  //Remplissage de la structure
  data_stream->fft_h = res_fft;

  //allocation des buffers utilisés dans le stream
  data_stream->fft_x = (double*)calloc( 2*FFTOrder,sizeof(*data_stream->fft_x));
  data_stream->res_ifft = (double *)calloc(2*FFTOrder,sizeof(double));
  //tmp pour le signal
  data_stream->tmp_x = (data_stream->fft_x + FFTOrder);
  data_stream->tmp3 = (data_stream->res_ifft + FFTOrder);
  
#endif  
  
  std::cout << "init OK" << std::endl;
					     
  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void *)data_stream, &options );
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    exit( 1 );
  }

  //data_stream->L = 512;//bufferFrames * channels * sizeof( MY_TYPE );

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

  //Libération de mémoire
  free (h_filter);

  //Liberer la struc data_stream
  free(data_stream);
  return 0;
}


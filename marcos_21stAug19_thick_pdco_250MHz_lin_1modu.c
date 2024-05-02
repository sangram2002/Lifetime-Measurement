#include <fcntl.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#define PMODE 0666
#define EVENT_INFO_LENGTH             68      // Information length for each event
//#define BIG 88888888
//#define BIG 17179869184 works april 8
#define BIG 1500000
//#define Odsize 262144
//16 channels 16 * 2* 8192 + 4 addback spectra
//#define Odsize 344064
//#define Odsize 327680
//3 modules 16 channels in each 3 * 16 * 2 * 8192 + 12 addback spectra (16bit)
//#define Odsize 983040
//4 modules 16 channels in each 4 * 16 * 2 * 8192 + 16 addback spectra (16bit)
//4*16*2*8192+16*2*8192=1310720 ! 16-bit
//#define Odsize 1310720  
//4*16*4*8192+16*4*8192=2621440 ! 32-bit
//#define Odsize 2621440
//4*16*4*8192+16*4*8192+8*4*8192=2621440+262144=2883584 ! 32-bit
//#define Odsize 2883584
//6*16*4*8192+24*4*8192+8*4*8192=4194304 ! 32-bit
#define Odsize 4194304
#define Odsizecr 32768
//#define TOdsize 16384
//#define TOdsize 262144
//#define TOdsize 786432 //48 Channels
//#define TOdsize 1048576 //64 Channels (16 bit)
//#define TOdsize 2097152 //64 Channels (32 bit)
//#define TOdsize 2097152 //64 Channels (32 bit)
#define TOdsize 3145728 //96 Channels (32 bit)
#define dim 4096

//Modified on 4th Feb 2011 for incorporation of DCO Matrix by R. Palit
//Modified on 15th Nov 2011 for incorporation of 6 modules by R. Palit

/////////////////////////////////////////////////////////////////////
float Offset[256],Slope[256],Quad[256],cube[256],TOffset[256];
int OptGainMatch;
void GainMatch(int Param);
//////////////////////////////////////////////////////////////////////
void GainMatch(int Param)
{
int i,Option;
char GainFile[40];
long int rand_number;
time_t *t;
FILE *Gp;

printf("Do you want to gain match?(Yes=1,No=0)"); scanf("%d",&OptGainMatch);
if (OptGainMatch)
   {
   rand_number = time(NULL);
   srand(rand_number);                             /*Random Number Seed*/
   printf("Random number seed is generated\n");
   printf("Gain match parameters from Keyboard or File(0 or 1)?");
   scanf("%d",&Option);
   if (Option)
      {
      printf("Name of Gain matching parameter file?"); scanf("%s",GainFile);
          if (access(GainFile,0))
             {
             printf("File: '%s' does not exist!\n",GainFile);
             printf("Try another gain file name !\n");
             printf("Gain file name?"); scanf("%s",GainFile);
             }
      Gp=fopen(GainFile,"r");
      for (i=0;i<Param;i++)
        {
        fscanf(Gp," %g %g %g %g %g",&Offset[i],&Slope[i],&Quad[i],&cube[i],&TOffset[i]);
        printf(" Offset= %g Slope= %g  Quad= %g cube= %g TOffset= %g\n",Offset[i],Slope[i],Quad[i],cube[i],TOffset[i]);
        }

      }
     else
      {
      printf("Enter the gain match parameters:\n");
      for (i=0;i<Param;i++)
         {
         printf("Give gain offset, slope, quad, cube anf TOffset for detector %d ", i+1);
         scanf("%g %g %g %g %g",&Offset[i],&Slope[i],&Quad[i],&cube[i], &TOffset[i]);
         printf(" Offset= %g Slope= %g  Quad= %g cube= %g TOffset= %g\n",Offset[i],Slope[i],Quad[i],cube[i],TOffset[i]);
         }
      }
   }
if (!OptGainMatch)
   {
   for (i=0;i<Param;i++) { Offset[i]=0.; Slope[i]=1.; Quad[i]=0.; cube[i]=0.; TOffset[i]=0.;}
   }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/





main()
{

int i,Rt,Odf,TOdf,Tdf,Tdf1;
int Rtt;
int j,l,k;
int jnum;
int Param;
unsigned long t1,t2,t3,t4,t5,t6,t7;
double x,y,tfpga,tcfd;
unsigned long timespec;
double timespec1;
unsigned z;
unsigned long z1;
char ListFileName[150];
char Bfile1d[20];
unsigned long eventdata, headerlength, eventlength;
unsigned long EVTTIME_LO, EVTTIME_HI, CFD_FRAC_TIME;
double CFD_TIME, tlow;
unsigned long TotalWords, TotalSkippedWords, NumEvents,NumEvents1;
unsigned long TotNumEvents;
//unsigned long *EventInformation;
unsigned long *Cha1,*En1;
int FileNumber2bScanned;
long int rand_number;

unsigned long delta12;
double time0,deltat,CTm[96];
double timewin;
unsigned long CEvent,CEn[96];
unsigned long EAdd[24],GEnD[96],x1,x2;
float of[96],sl[96];
double *tcfd1;



//unsigned  short *Twod;
int *Twod;

int cla[24], clsum, clov1f,clov2f,clov3f;
int clov4f,clov5f,clov6f,clov7f,clov8f,clov9f,clov10f,clov11f,clov12f,clov13f;
int MergeFileNum;

int a[96],FP;
int Tdsize; 

//unsigned short *sp;
int *sp;
//unsigned short *tsp;
int *tsp;
unsigned long En,GEn,GEn0,GEn1,GEn2,GEn3,GEnA,Channel,Slot,Crate;
unsigned long GEn01,GEn11,GEn21,GEn31,GEnA1;
FILE *ListModeFile = NULL;
FILE *Lf;

//Lamps spectrum
int *spcr;
int Odfcr,Rtcr;



//Cube parameters
int TDF;
unsigned short CubBuf[8192];
int BufSize,CubBufEvent;
static short sbuf[32768];

//TwodDCO matrix parameters
int m1,m2;
unsigned long EAddx[3],EAddy[4];
int *TwodDCO;
int Tdfdco;

//TwodPolarization Matrix parameters
unsigned long EPar[10],EPer[10];
unsigned long EAdddum[24];
int *TwodPol;
int *TwodPol1;
int Tdfpol;
int Tdfpol1;

///////////////////////////////////////////////////////////////////////
//For 16-bit 4k x 4k 2d-spectrum
//Tdsize=2* (dim) * (dim) ;
//Twod = (unsigned short *)malloc(Tdsize);
//for(i=0;i<(Tdsize/2);i++) {Twod[i]=0;}
////////////////////////////////////////////////////////////////////
////For 32-bit 4k x 4k 2d-spectrum
Tdsize=4* (dim) * (dim) ;

///////////////////////////////////////////////////////////////////

Twod = (unsigned int *)malloc(Tdsize);
for(i=0;i<(Tdsize/4);i++) {Twod[i]=0;}
///////////////////////////////////////////////////////////////////

TwodDCO = (unsigned int *)malloc(Tdsize);
for(i=0;i<(Tdsize/4);i++) {TwodDCO[i]=0;}

///////////////////////////////////////////////////////////////////

TwodPol = (unsigned int *)malloc(Tdsize);
TwodPol1 = (unsigned int *)malloc(Tdsize);
for(i=0;i<(Tdsize/4);i++) {TwodPol[i]=0; TwodPol1[i]=0;}

///////////////////////////////////////////////////////////////////
//
//
//EventInformation = (unsigned long *)malloc(BIG);
//if (EventInformation == NULL) { printf("memory allocation problem\n"); exit(0);}
En1 = (unsigned long *)malloc(BIG);
if (En1 == NULL) { printf("memory allocation problem\n"); exit(0);}
Cha1 = (unsigned long *)malloc(BIG);
if (Cha1 == NULL) { printf("memory allocation problem\n"); exit(0);}
tcfd1 = (double *)malloc(BIG);
if (tcfd1 == NULL) { printf("memory allocation problem\n"); exit(0);}
//sp = (unsigned short *)malloc(Odsize);
//for(i=0;i<(Odsize/2);i++) {sp[i]=0;}
sp = (unsigned int *)malloc(Odsize);
for(i=0;i<(Odsize/4);i++) {sp[i]=0;}


spcr = (unsigned int *)malloc(Odsizecr);
for(i=0;i<(Odsizecr/4);i++) {spcr[i]=0;}



//tsp = (unsigned short *)malloc(TOdsize);
//for(i=0;i<(TOdsize/2);i++) {tsp[i]=0;}
tsp = (unsigned int *)malloc(TOdsize);
for(i=0;i<(TOdsize/4);i++) {tsp[i]=0;}

Param = 96;

	printf("set0\n");
system("clear");
printf("Multi pARameter time stamped based COincidence Search program\n");
printf("MARCOS --- R. Palit written: 18 Dec 2009: Sorting of XIA list files\n");

printf("Total number of Files to be scanned?"); scanf("%d",&FileNumber2bScanned);
printf("Give the time window?");scanf("%lf",&timewin);
printf("Time window = %lf \n",timewin);
TDF=creat("xx.dmp",0666);


printf("*.his 1d-file?"); scanf("%s",Bfile1d);


      Lf=fopen("Listfile.log","w");


MergeFileNum=0;
do
{
	while(1)
	{
	printf("XIA List File Name?"); scanf("%s",ListFileName);
	printf("%s\n",ListFileName);
        	if (access(ListFileName,0))
         	{
		printf("File: '%s' doesn't exist!\n",ListFileName);
		printf("Try another XIA List File name\n");
		printf("XIA List File Name?"); scanf("%s",ListFileName);
		}
		if (!(access(ListFileName,0))) break; 
         }

ListModeFile = fopen(ListFileName, "rb");

        fprintf(Lf," %s\n",ListFileName);


rand_number = time(NULL);
   srand(rand_number);                             /*Random Number Seed*/
   printf("Random number seed is generated\n");

if(ListModeFile == NULL) break;

if(ListModeFile != NULL)
        {

	printf("set1\n");
                // Get file length
                fseek(ListModeFile, 0, SEEK_END);
                TotalWords = (ftell(ListModeFile) + 1) / 4;
                printf("TotalWords = %d\n", TotalWords);                                                                                                                                       
                fseek(ListModeFile, 0, SEEK_SET);
//                fseek(ListModeFile, 24915520, SEEK_SET);
                                                                                                                                                       
                TotalSkippedWords = 0;
		NumEvents = 0;

GainMatch(Param);

////////////////////////////////

clov1f = 0;
clov2f = 0;
clov3f = 0;
clov4f = 0;
clov5f = 0;
clov6f = 0;
clov7f = 0;
clov8f = 0;
clov9f = 0;
clov10f = 0;
clov11f = 0;
clov12f = 0;
clov13f = 0;

////////////////////////////////


jnum = 0;
do
{

printf("jnum = %d\n", jnum);
	printf("TotalSkippedWords = %d\n",TotalSkippedWords);
NumEvents=0;
for (i=0;i<15000;i++)
{
tcfd1[i] = 0;
En1[i] = 0;
Cha1[i] = 0;
}

	   do
	   {

Rt =0;
      Rt=fread(&eventdata, 4, 1, ListModeFile);

	if (Rt<1) { printf("Rt = %d\n",Rt); break;}
	// Event #
	Channel = (eventdata & 0xF);
        // Slot #
        Slot = (eventdata & 0xF0) >> 4;
        // Crate #
  //      Crate = (eventdata & 0xF00) >>8;
      //  printf("Crate = %lu\n",Crate);
        // Header length
        headerlength = (eventdata & 0x1F000) >> 12;
       // Event length
        eventlength = (eventdata & 0x7FFE0000) >> 17;
	if (eventlength>4) printf("eventlength = %d\n",eventlength);
      // Finish code

	         	fread(&eventdata, 4, 1, ListModeFile);
                        // EventTime_Low
			EVTTIME_LO = 0;
			EVTTIME_LO = eventdata;
t1 = eventdata & 0xFFFF;
t2 = (eventdata & 0xFFFF0000) >> 16;


                        fread(&eventdata, 4, 1, ListModeFile);
                        // EventTime_High
			EVTTIME_HI = (eventdata & 0xFFFF);
	t3 = (eventdata & 0xFFFF);
//			printf("EVTTIME_HI = %d\n", EVTTIME_HI);
			// CFD TIME
		//	CFD_FRAC_TIME = (eventdata & 0xFFFF0000)>>16;  100MHz
			 
//t4=(eventdata & 0xFFFF0000)>>16; 100MHz
//tcfd = (t3 * pow (x,z) + t2 * pow(x,y) + t1 + t4/65536.)*10.; 100
//250MHZ

	CFD_FRAC_TIME = (eventdata & 0x3FFF0000)>>16;
    t4=(eventdata & 0x3FFF0000)>>16;
    t5=(eventdata & 0x40000000)>>30;
    t6=(eventdata & 0x80000000)>>31;

x = 2.0;
y = 16.0;
z = 32.0;


        if(t6 == 0)
			{
			tcfd = ((t3 * pow (x,z) + t2 * pow(x,y)+t1)*2 - t5 + t4/16384.0)*4.;	
			}	
		else if(t6 == 1)
		    {
			tcfd = ((t3 * pow (x,z) + t2 * pow(x,y) + t1)*2)*4.; 	
			}	

                        fread(&eventdata, 4, 1, ListModeFile);
                        // Event Energy
			En = (eventdata & 0xFFFF);

			TotalSkippedWords += eventlength;
			NumEvents ++;
//			fseek(ListModeFile, 0, SEEK_CUR);
			fseek(ListModeFile, (eventlength - 4) * 4, SEEK_CUR);
//if ((Slot == 2) && (Channel == 3))tcfd = tcfd - 400.;
//if ((Slot == 4) && (Channel == 0))tcfd = tcfd - 400.;
//if ((Slot == 4) && (Channel == 2))tcfd = tcfd - 400.;
tcfd1[NumEvents-1] = tcfd;
En1[NumEvents-1] = En;
Cha1[NumEvents-1] = Channel + (Slot -2)*16;
//printf("value of tcfd1 %f, En1 %lu and Cha1 %lu \n",tcfd1[NumEvents-1],En1[NumEvents-1] ,Cha1[NumEvents-1]); 


//if (TotalSkippedWords%1000000) {printf("TotalSkippedWords = %d\n",TotalSkippedWords);} 
//printf("try1 %d %d %f\n",NumEvents-1,En,(tcfd-tcfd1[0])/1000);
//if ((Cha1[NumEvents-1] ==5) && (En > 30) && (En< 8000)){GEn0 = En; sp[GEn0]++;}
	   } while (NumEvents < 15000 );

NumEvents1 = NumEvents;
//while( NumEvents < 150000 );

//TotalSkippedWords < TotalWords
   printf("Random number seed is generated11 %d\n", NumEvents1);
//////////////////////////////////////////////////////////////
//
//
//GainMatch(Param);
//
//
//////////////////////////////////////////////////////////////

NumEvents = 0;
time0 = tcfd1[0];
CEvent = 0;
GEn =0;
for(i=0;i<96;i++)
{a[i]=0;CEn[i]=0;CTm[i]=0.0;}
i = Cha1[0];
CEn[i] = En1[0];
CTm[i] = tcfd1[0];
a[i]=1;
/*
for(i=0;i<96;i++)
{
printf("Before loop: value of Channel %d and CTm %f \n",i,CTm[i]); 
}
*/
NumEvents = 1;

	printf("loop\n");
	printf("Num %d %f\n", NumEvents, tcfd1[0]);


CubBufEvent=0;

for (i=0;i<32768;i++) {sbuf[i]=0;}



do
{

 deltat = fabs(tcfd1[NumEvents]-time0);	

// if(deltat<1500.) 
 if(deltat<timewin) 
// if(deltat<200.) 
 {

//
//	z1 = timespec;
//	
//        if (deltat > 0) printf("try2 %d %f\n",NumEvents, tcfd1[NumEvents]-time0+2000);
	timespec = tcfd1[NumEvents]-time0+2000 ;
//        if (deltat > 0) printf("try3 %d %f\n",NumEvents, timespec);
	i=Cha1[NumEvents];
	a[i] = 1;
	CEn[i]=En1[NumEvents];
	CTm[i]=tcfd1[NumEvents];
//	printf("value of Channel %d and CTm %f \n",i,CTm[i]); // for time testing with pulsar
//	printf("value of Channel 0 and CTm[0] %f \n",CTm[0]); // for time testing with pulsar
	if ((deltat > 0) && (CEn[i] > 20) ) tsp[timespec]++;
 }

 else 
 {
 FP =0;
 for(i=0;i<96;i++) {FP = FP + a[i];}
// printf("CEvent = %d FP = %d\n",CEvent, FP);
//
// 
//
/////////////////////////////////////////

for (i=0;i<24;i++) cla[i] = 0;

for (i=0;i<24;i++)
{
cla[i]=0;
j = 4*i;
if ((a[j]==1) || (a[j+1]==1) || (a[j+2]==1) || (a[j+3]==1)) cla[i] = 1;
}

clsum=0;
for (i=0;i<24;i++) clsum = clsum + cla[i];

if (clsum ==1) clov1f = clov1f++ ;
if (clsum ==2) clov2f = clov2f++ ;
if (clsum ==3) clov3f = clov3f++ ;
if (clsum ==4) clov4f = clov4f++ ;
if (clsum ==5) clov5f = clov5f++ ;
if (clsum ==6) clov6f = clov6f++ ;
if (clsum ==7) clov7f = clov7f++ ;
if (clsum ==8) clov8f = clov8f++ ;
if (clsum ==9) clov9f = clov9f++ ;
if (clsum ==10) clov10f = clov10f++ ;
if (clsum ==11) clov11f = clov11f++ ;
if (clsum ==12) clov12f = clov12f++ ;
if (clsum ==13) clov13f = clov13f++ ;


/////////////////////////////////////////
//if ( (a[i] == 1) && (CEn[i] > 100.) && (CEn[i]<28000.) )
//{
/*
CEn[0]=CEn[0]*1.0114;CEn[1]=CEn[1]*1.0114;CEn[2]=CEn[2]*1.0114;CEn[3]=CEn[3]*1.0114;
CEn[4]=CEn[4]*1.0114;CEn[5]=CEn[5]*1.0114;CEn[6]=CEn[6]*1.0114;CEn[7]=CEn[7]*1.0114;
CEn[8]=CEn[8]*1.0114;CEn[9]=CEn[9]*1.0114;CEn[10]=CEn[10]*1.0114;CEn[11]=CEn[11]*1.0114;
CEn[12]=CEn[12]*1.0104;CEn[13]=CEn[13]*1.0104;CEn[14]=CEn[14]*1.0104;CEn[15]=CEn[15]*1.0104;
CEn[16]=CEn[16]*1.010;CEn[17]=CEn[17]*1.010;CEn[18]=CEn[18]*1.010;CEn[19]=CEn[19]*1.010;
CEn[20]=CEn[20]*1.0109;CEn[21]=CEn[21]*1.0109;CEn[22]=CEn[22]*1.0109;CEn[23]=CEn[23]*1.0109;
CEn[24]=CEn[24]*1.0067;CEn[25]=CEn[25]*1.0067;CEn[26]=CEn[26]*1.0067;CEn[27]=CEn[27]*1.0067;
CEn[28]=CEn[28]*1.0073;CEn[29]=CEn[29]*1.0073;CEn[30]=CEn[30]*1.0073;CEn[31]=CEn[31]*1.0073;
CEn[32]=CEn[32]*1.0078;CEn[33]=CEn[33]*1.0078;CEn[34]=CEn[34]*1.0078;CEn[35]=CEn[35]*1.0078;
CEn[36]=CEn[36]*.99234;CEn[37]=CEn[37]*.99234;CEn[38]=CEn[38]*.99234;CEn[39]=CEn[39]*.99234;
CEn[40]=CEn[40]*1.0015;CEn[41]=CEn[41]*1.0015;CEn[42]=CEn[42]*1.0015;CEn[43]=CEn[43]*1.0015;
CEn[44]=CEn[44]*1.0021;CEn[45]=CEn[45]*1.0021;CEn[46]=CEn[46]*1.0021;CEn[47]=CEn[47]*1.0021;
//CEn[48]=CEn[48]*1.0;CEn[49]=CEn[49]*1.0;CEn[50]=CEn[50]*1.0;CEn[51]=CEn[51]*1.0;
CEn[52]=CEn[52]*1.0015;CEn[53]=CEn[53]*1.0015;CEn[54]=CEn[54]*1.0015;CEn[55]=CEn[55]*1.0015;
CEn[56]=CEn[56]*1.0021;CEn[57]=CEn[57]*1.0021;CEn[58]=CEn[58]*1.0021;CEn[59]=CEn[59]*1.0021;
//CEn[60]=CEn[60]*1.0;CEn[61]=CEn[61]*1.0;CEn[62]=CEn[62]*1.0;CEn[63]=CEn[63]*1.0;
CEn[64]=CEn[64]*0.9969;CEn[65]=CEn[65]*0.9969;CEn[66]=CEn[66]*0.9969;CEn[67]=CEn[67]*0.9969;
CEn[68]=CEn[68]*.99641;CEn[69]=CEn[69]*.99641;CEn[70]=CEn[70]*.99641;CEn[71]=CEn[71]*.99641;
CEn[72]=CEn[72]*.99234;CEn[73]=CEn[73]*.99234;CEn[74]=CEn[74]*.99234;CEn[75]=CEn[75]*.99234;
CEn[76]=CEn[76]*.99234;CEn[77]=CEn[77]*.99234;CEn[78]=CEn[78]*.99234;CEn[79]=CEn[79]*.99234;
CEn[80]=CEn[80]*.99285;CEn[81]=CEn[81]*.99285;CEn[82]=CEn[82]*.99285;CEn[83]=CEn[83]*.99285;
*/
//}

for(i=0;i<96;i++) {GEnD[i]=0.;}

 for(i=0;i<96;i++)
{
if ( (a[i] == 1) && (CEn[i] > 50.) && (CEn[i]<64000.) ) 
{ 
//CEn[i]=0.5*CEn[i]; 

if (i!=6 && i!=11){
CEn[i]=0.5*CEn[i];}
else {CEn[i]=0.2*CEn[i];}
GEn = 0; GEn = Offset[i] + Slope[i]*CEn[i] + Quad[i]*CEn[i]*CEn[i] + cube[i]*CEn[i]*CEn[i]*CEn[i]
               + (float) rand()/(float) RAND_MAX;



 

//  GEn = 0.25*GEn;
  sp[i*8192 + GEn]++; 
  GEnD[i] = GEn;
}
}

//printf("value of CTm[0] %f and CTm[4] %f \n",CTm[0],CTm[4]); 

//Normal timing
//if (CEn[8] > 100.)

for(j=0;j<16;j++){
if ((GEnD[j] > 73) && (GEnD[j] < 89) ) //&& (CEn[4] > 50.0) )
{
for(i=j+1;i<16;i++)
{
timespec1 = CTm[j]-CTm[i]; 
timespec1 = 5.*(CTm[j]-CTm[i]); 
timespec = fabs(CTm[j]-CTm[i]); 

 
if ( (timespec > -1.) && (GEnD[i] > 347) && (GEnD[i]<365) ) 
{
	timespec = 8192*1 + timespec1+ 2000 - TOffset[i] + TOffset[j]; tsp[timespec]++;
	//printf("value of timespec1 %d \n",TOffset[i]); 
}
}
}
}


for(j=0;j<16;j++){
if ((GEnD[j] > 73) && (GEnD[j] < 89) ) //&& (CEn[4] > 50.0) )
{
for(i=j+1;i<16;i++)
{
timespec1 = CTm[j]-CTm[i]; 
timespec1 = 5.*(CTm[j]-CTm[i]); 
timespec = fabs(CTm[j]-CTm[i]); 

 
if ( ((timespec > -1.) && (GEnD[i] > 323) && (GEnD[i]<327))  || ((timespec > -1.) && (GEnD[i] > 421) && (GEnD[i]<433) )) 
{
	timespec = 8192*2 + timespec1+ 2000 - TOffset[i] + TOffset[j]; tsp[timespec]++;
	//printf("value of timespec1 %d \n",TOffset[i]); 
}
}
}
}


for(j=0;j<16;j++){
if (((GEnD[j] > 54) && (GEnD[j] < 63)) || ((GEnD[j] > 100) && (GEnD[j] < 106))) //&& (CEn[4] > 50.0) )
{
for(i=j+1;i<16;i++)
{
timespec1 = CTm[j]-CTm[i]; 
timespec1 = 5.*(CTm[j]-CTm[i]); 
timespec = fabs(CTm[j]-CTm[i]); 

 
if ( (timespec > -1.) && (GEnD[i] > 347) && (GEnD[i]<365) ) 
{
	timespec = 8192*3 + timespec1+ 2000 - TOffset[i] + TOffset[j]; tsp[timespec]++;
	//printf("value of timespec1 %d \n",TOffset[i]); 
}
}
}
}



for(j=0;j<16;j++){
if (((GEnD[j] > 54) && (GEnD[j] < 63)) || ((GEnD[j] > 100) && (GEnD[j] < 106))) //&& (CEn[4] > 50.0) )
{
for(i=j+1;i<16;i++)
{
timespec1 = CTm[j]-CTm[i]; 
timespec1 = 5.*(CTm[j]-CTm[i]); 
timespec = fabs(CTm[j]-CTm[i]); 

 
if ( ((timespec > -1.) && (GEnD[i] > 323) && (GEnD[i]<327))  || ((timespec > -1.) && (GEnD[i] > 421) && (GEnD[i]<433) )) 
{
	timespec = 8192*4 + timespec1+ 2000 - TOffset[i] + TOffset[j]; tsp[timespec]++;
	//printf("value of timespec1 %d \n",TOffset[i]); 
}
}
}
}


/*
// 60Co gated
if ((CEn[2] >2000.)&&(CEn[2]<3400.)) //1333 keV
{
for(i=3;i<96;i++)
{
timespec1 = CTm[2]-CTm[i]; 
timespec1 = 5.*(CTm[2]-CTm[i]); 
timespec = fabs(CTm[2]-CTm[i]); 
if ((timespec > -1.) && (CEn[i] > 2000.) && (CEn[i]<3400.) ) {timespec = 8192*i + timespec1+ 2000; tsp[timespec]++;}
//printf("xxxxx %f\n", timespec);
//if (timespec > 0) {timespec = 8192*i + timespec+ 2000.; tsp[timespec]++;}
}
}
*/

/*
// 60Co or 133Ba or 152 Eu gated with calibration file
if ((GEnD[2]> 60.)&&(GEnD[2]<6000.)) //1333 keV
{
for(i=3;i<96;i++)
{
timespec1 = CTm[2]-CTm[i]; 
timespec1 = 5*(CTm[2]-CTm[i]); 
timespec = fabs(CTm[2]-CTm[i]); 
if ((timespec > -1.) && (GEnD[i] > 300.) && (GEnD[i]<340.) ) {timespec = 8192*i + timespec1+ 2000; tsp[timespec]++;}
//printf("xxxxx %f\n", timespec);
//if (timespec > 0) {timespec = 8192*i + timespec+ 2000.; tsp[timespec]++;}
}
}
*/

/*
// Ba gated
if ((CEn[0] > 650.)&&(CEn[0]<729)) //355 keV
{
for(i=1;i<96;i++)
{
timespec1 = CTm[0]-CTm[i]; 
timespec1 = 10*(CTm[0]-CTm[i]); 
timespec = fabs(CTm[0]-CTm[i]); 
if ((timespec > -1.) && (CEn[i] > 150) && (CEn[i]<200) ) {timespec = 8192*i + timespec1+ 2000; tsp[timespec]++;}
//printf("xxxxx %f\n", timespec);
//if (timespec > 0) {timespec = 8192*i + timespec+ 2000.; tsp[timespec]++;}
}
}

*/
/*
// Eu gated
if ((CEn[0] > 450.)&&(CEn[0]<500.)) //244 keV
{
for(i=1;i<96;i++)
{
timespec1 = CTm[0]-CTm[i]; 
timespec1 = 10*(CTm[0]-CTm[i]); 
timespec = fabs(CTm[0]-CTm[i]); 
if ((timespec > -1.) && (CEn[i] >230.) && (CEn[i]<290.) ) {timespec = 8192*i + timespec1+ 2000; tsp[timespec]++;}
//printf("xxxxx %f\n", timespec);
//if (timespec > 0) {timespec = 8192*i + timespec+ 2000.; tsp[timespec]++;}
}
}
*/


EAdd[0]=0;
EAdd[1]=0;
EAdd[2]=0;
EAdd[3]=0;
EAdd[4]=0;
EAdd[5]=0;
EAdd[6]=0;
EAdd[7]=0;
EAdd[8]=0;
EAdd[9]=0;
EAdd[10]=0;
EAdd[11]=0;
EAdd[12]=0;
EAdd[13]=0;
EAdd[14]=0;
EAdd[15]=0;
EAdd[16]=0;
EAdd[17]=0;
EAdd[18]=0;
EAdd[19]=0;
EAdd[20]=0;
EAdd[21]=0;
EAdd[22]=0;
EAdd[23]=0;
/*
EAdd[0]=GEnD[0]+GEnD[1]+GEnD[2]+GEnD[3];      // -23 deg PC1  24 Nov 11 - 23 Dec 11
EAdd[1]=GEnD[4]+GEnD[5]+GEnD[6]+GEnD[7];      // -23 deg PC2
EAdd[2]=GEnD[8]+GEnD[9]+GEnD[10]+GEnD[11];    // -23 deg PC3
EAdd[3]=GEnD[12]+GEnD[13]+GEnD[14]+GEnD[15];  // -40 deg PC4
EAdd[4]=GEnD[16]+GEnD[17]+GEnD[18]+GEnD[19];  // -40 deg PC5
EAdd[5]=GEnD[20]+GEnD[21]+GEnD[22]+GEnD[23];  // -40 deg PC6
EAdd[6]=GEnD[24]+GEnD[25]+GEnD[26]+GEnD[27];  // -65 deg PC7
EAdd[7]=GEnD[28]+GEnD[29]+GEnD[30]+GEnD[31];  // -65 deg PC8
EAdd[8]=GEnD[32]+GEnD[33]+GEnD[34]+GEnD[35];  // -65 deg PC9
EAdd[9]=GEnD[36]+GEnD[37]+GEnD[38]+GEnD[39];  // +40 deg PC21
EAdd[10]=GEnD[40]+GEnD[41]+GEnD[42]+GEnD[43]; // 90 deg PC11  par 40+41, 42+43 per 40+42, 41+43
EAdd[11]=GEnD[44]+GEnD[45]+GEnD[46]+GEnD[47]; // 90 deg PC12  par 44+45, 46+47 per 44+47, 45+46
EAdd[12]=GEnD[48]+GEnD[49]+GEnD[50]+GEnD[51]; // Empty  
EAdd[13]=GEnD[52]+GEnD[53]+GEnD[54]+GEnD[55]; // 90 deg PC14  par 52+53, 54+55 per 52+55, 53+54
EAdd[14]=GEnD[56]+GEnD[57]+GEnD[58]+GEnD[59]; // 90 deg PC15  par 56+57, 58+59 per 56+59, 57+58
EAdd[15]=GEnD[60]+GEnD[61]+GEnD[62]+GEnD[63]; // LaBr2 in 2 channels  
EAdd[16]=GEnD[64]+GEnD[65]+GEnD[66]+GEnD[67]; // +65 deg PC17
EAdd[17]=GEnD[68]+GEnD[69]+GEnD[70]+GEnD[71]; // +65 deg PC16
EAdd[18]=GEnD[72]+GEnD[73]+GEnD[74]+GEnD[75]; // +40 deg PC19
EAdd[19]=GEnD[76]+GEnD[77]+GEnD[78]+GEnD[79]; // +40 deg PC20
EAdd[20]=GEnD[80]+GEnD[81]+GEnD[82]+GEnD[83]; // Not connected
EAdd[21]=GEnD[84]+GEnD[85]+GEnD[86]+GEnD[87]; // Not connected
EAdd[22]=GEnD[88]+GEnD[89]+GEnD[90]+GEnD[91]; // Not connected
EAdd[23]=GEnD[92]+GEnD[93]+GEnD[94]+GEnD[95]; // Not connected   
*/

EAdd[0]=GEnD[0]; 
EAdd[1]=GEnD[1]; 
EAdd[2]=GEnD[2]; 
EAdd[3]=GEnD[3]; 
EAdd[4]=GEnD[4]; 
EAdd[5]=GEnD[5]; 
EAdd[6]=GEnD[6]; 
EAdd[7]=GEnD[7]; 
EAdd[8]=GEnD[8]; 
EAdd[9]=GEnD[9]; 
EAdd[10]=GEnD[10]; 
EAdd[11]=GEnD[11]; 
EAdd[12]=GEnD[12]; 
EAdd[13]=GEnD[13]; 
EAdd[14]=GEnD[14]; 
EAdd[14]=GEnD[14]; 
EAdd[15]=GEnD[15]; 

for(i=0;i<24;i++) {if(EAdd[i] > 50) sp[(Param+i)*8192+EAdd[i]]++;}
for(i=0;i<24;i++) {if((EAdd[i] > 50)&&(EAdd[i] < 8100) ) spcr[EAdd[i]]++;}


for(l=0;l<16;l++)
              {
                for(k=l+1;k<16;k++)
                {
                  x=0; y=0; z=0;
if((EAdd[l]>50) && (EAdd[l]<dim) && (EAdd[k]>50) && (EAdd[k]<dim))
                 {x=EAdd[l]; y=EAdd[k]; z= x+4096*y; Twod[z]++;}
                }
              }

//===============================DCO======================================
for(m1=0;m1<3;m1++){EAddx[m1]=0;}
for(m2=0;m2<4;m2++){EAddy[m2]=0;}

if (EAdd[0]>0) EAddx[0]=EAdd[0];
if (EAdd[1]>0) EAddx[1]=EAdd[1];
if (EAdd[2]>0) EAddx[2]=EAdd[2];

if (EAdd[10]>0) EAddy[0]=EAdd[10];
if (EAdd[11]>0) EAddy[1]=EAdd[11];
if (EAdd[13]>0) EAddy[2]=EAdd[13];
if (EAdd[14]>0) EAddy[3]=EAdd[14];


for(m1=0;m1<3;m1++)
              {
                for(m2=0;m2<4;m2++)
                {
                  x=0; y=0; z=0;
if((EAddx[m1]>50) && (EAddx[m1]<dim) && (EAddy[m2]>50) && (EAddy[m2]<dim))
                 {x=EAddx[m1]; y=EAddy[m2]; z= x+4096*y; TwodDCO[z]++;}
                }
              }

//============================Polarization=====================================

EPar[0]=0.;
if ((GEnD[40]>20.)&&(GEnD[41]>20.)&&(GEnD[42]<20.)&&(GEnD[43]<20.))
{EPar[0] = GEnD[40]+GEnD[41];}
if ((GEnD[40]<20.)&&(GEnD[41]<20.)&&(GEnD[42]>20.)&&(GEnD[43]>20.))
{EPar[0] = GEnD[42]+GEnD[43];}

EPer[0]=0.;
if ((GEnD[40]>20.)&&(GEnD[41]<20.)&&(GEnD[42]>20.)&&(GEnD[43]<20.))
{EPer[0] = GEnD[40]+GEnD[42];}
if ((GEnD[40]<20.)&&(GEnD[41]>20.)&&(GEnD[42]<20.)&&(GEnD[43]>20.))
{EPer[0] = GEnD[41]+GEnD[43];}

EPar[1]=0.;
if ((GEnD[44]>20.)&&(GEnD[45]>20.)&&(GEnD[46]<20.)&&(GEnD[47]<20.))
{EPar[1] = GEnD[44]+GEnD[45];}
if ((GEnD[44]<20.)&&(GEnD[45]<20.)&&(GEnD[46]>20.)&&(GEnD[47]>20.))
{EPar[1] = GEnD[46]+GEnD[47];}

EPer[1]=0.;
if ((GEnD[44]>20.)&&(GEnD[45]<20.)&&(GEnD[46]<20.)&&(GEnD[47]>20.))
{EPer[1] = GEnD[44]+GEnD[47];}
if ((GEnD[44]<20.)&&(GEnD[45]>20.)&&(GEnD[46]>20.)&&(GEnD[47]<20.))
{EPer[1] = GEnD[45]+GEnD[46];}

EPar[2]=0.;
if ((GEnD[52]>20.)&&(GEnD[53]>20.)&&(GEnD[54]<20.)&&(GEnD[55]<20.))
{EPar[2] = GEnD[52]+GEnD[53];}
if ((GEnD[52]<20.)&&(GEnD[53]<20.)&&(GEnD[54]>20.)&&(GEnD[55]>20.))
{EPar[2] = GEnD[54]+GEnD[55];}

EPer[2]=0.;
if ((GEnD[52]>20.)&&(GEnD[53]<20.)&&(GEnD[54]<20.)&&(GEnD[55]>20.))
{EPer[2] = GEnD[52]+GEnD[55];}
if ((GEnD[52]<20.)&&(GEnD[53]>20.)&&(GEnD[54]>20.)&&(GEnD[55]<20.))
{EPer[2] = GEnD[53]+GEnD[54];}

EPar[3]=0.;
if ((GEnD[56]>20.)&&(GEnD[57]>20.)&&(GEnD[58]<20.)&&(GEnD[59]<20.))
{EPar[3] = GEnD[56]+GEnD[57];}
if ((GEnD[56]<20.)&&(GEnD[57]<20.)&&(GEnD[58]>20.)&&(GEnD[59]>20.))
{EPar[3] = GEnD[58]+GEnD[59];}

EPer[3]=0.;
if ((GEnD[56]>20.)&&(GEnD[57]<20.)&&(GEnD[58]<20.)&&(GEnD[59]>20.))
{EPer[3] = GEnD[56]+GEnD[59];}
if ((GEnD[56]<20.)&&(GEnD[57]>20.)&&(GEnD[58]>20.)&&(GEnD[59]<20.))
{EPer[3] = GEnD[57]+GEnD[58];}


for(i=0;i<4;i++){if(EPar[i]>20) sp[(Param+24+i)*8192+EPar[i]]++;  }
for(i=0;i<4;i++){if(EPer[i]>20) sp[(Param+28+i)*8192+EPer[i]]++;  }

for(m2=0;m2<24;m2++){EAdddum[m2]=EAdd[m2];}

for(m1=0;m1<4;m1++)
              {
if (m1==0) EAdddum[10] = 0;
if (m1==1) EAdddum[11] = 0;
if (m1==2) EAdddum[13] = 0;
if (m1==3) EAdddum[14] = 0;
                for(m2=0;m2<24;m2++)
                {
                  x=0; y=0; z=0;
if((EPar[m1]>50) && (EPar[m1]<dim) && (EAdddum[m2]>50) && (EAdddum[m2]<dim))
                 {x=EPar[m1]; y=EAdddum[m2]; z= x+4096*y; TwodPol[z]++;}
                  x=0; y=0; z=0;
if((EPer[m1]>50) && (EPer[m1]<dim) && (EAdddum[m2]>50) && (EAdddum[m2]<dim))
                 {x=EPer[m1]; y=EAdddum[m2]; z= x+4096*y; TwodPol1[z]++;}
                }

for(m2=0;m2<24;m2++){EAdddum[m2]=EAdd[m2];}

              }






//========================================================================

for (j=0;j<4;j++) CubBuf[j]=0;
if (clsum> 2)
{
l=CubBufEvent;
j=0;
for (i=0;i<24;i++)
    {
            if (EAdd[i]>50.) {CubBuf[j]=EAdd[i]; j++;}
    }

                if ((j>2) && (j<5))
                {
                for (j=0;j<4;j++) sbuf[4*l+4+j]=CubBuf[j];
                CubBufEvent++;
                }

// if ((j>1) && (j<4)) write(TDF,CubBuf,8);

}

//========================================================================





 CEvent++;
 time0 = tcfd1[NumEvents] ;

 for(i=0;i<96;i++){a[i]=0;CEn[i]=0;CTm[i]=0.0;}
 i=Cha1[NumEvents];
 a[i]=1;
 CEn[i] = En1[NumEvents];
 CTm[i] = tcfd1[NumEvents];
 }

 NumEvents ++;
} while (NumEvents < NumEvents1);

//======================Cube Buffer Update=================================


BufSize = 8*CubBufEvent;

printf("BufSize 1111 %d\n", BufSize);
write(TDF,&BufSize,4);
sbuf[0]=-8;
sbuf[1]=BufSize/2;
sbuf[2]=4;
sbuf[3]=4*CubBufEvent;

//for (i=0;i<CubBufEvent;i++)
//{
//printf("Blocksize : %d %d %d %d\n",sbuf[4*i+4],sbuf[4*i+5],sbuf[4*i+6],sbuf[4*i+7]); 
//}


write(TDF,sbuf,BufSize);

write(TDF,&BufSize,4);
//break;
////=======================End of Cube Buffer Update=========================



if (NumEvents1 < 15000) {printf("EOF File\n"); break;}

	printf("loop1\n");
jnum++;
TotNumEvents = jnum*NumEvents;

////////////////////////////////////////////////////////////////
//} while (TotNumEvents < 15000000);
}while (TotNumEvents < TotalWords);
//}while (TotalSkippedWords < TotalWords);

	printf("Final TotalSkippedWords = %d\n",TotalSkippedWords);
printf("clov fold %d %d %d %d %d %d %d %d %d %d %d %d %d\n", clov1f,clov2f,clov3f,clov4f,clov5f,clov6f,clov7f,clov8f,clov9f,clov10f,clov11f,clov12f,clov13f);
fprintf(Lf, "clov fold %d %d %d %d %d %d %d %d %d %d %d %d %d\n", clov1f,clov2f,clov3f,clov4f,clov5f,clov6f,clov7f,clov8f,clov9f,clov10f,clov11f,clov12f,clov13f);

        }
fclose(ListModeFile);

MergeFileNum++;


Tdf1=creat("2d_dum.his",0666);
Rt=write(Tdf1,Twod,Tdsize);
close(Tdf1);


} while (MergeFileNum<FileNumber2bScanned);



Odf=creat(Bfile1d,0666);
Rt=write(Odf,sp,Odsize);
close(Odf);

Odfcr=creat("add1.z1d",0666);
Rtcr=write(Odfcr,spcr,Odsizecr);
close(Odfcr);

TOdf=creat("t1d.his",0666);
Rtt=write(TOdf,tsp,TOdsize);
close(TOdf);

Tdfdco=creat("2d_dco.his",0666);
Rt=write(Tdfdco,TwodDCO,Tdsize);
close(Tdfdco);

Tdfpol=creat("2d_pol.his",0666);
Rt=write(Tdfpol,TwodPol,Tdsize);
close(Tdfpol);

Tdfpol1=creat("2d_pol1.his",0666);
Rt=write(Tdfpol1,TwodPol1,Tdsize);
close(Tdfpol1);

Tdf=creat("2d.his",0666);
Rt=write(Tdf,Twod,Tdsize);
close(Tdf);

close(TDF);
fclose(Lf);

return(0);
	
}

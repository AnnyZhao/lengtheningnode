
#include <math.h>
#include <stdio.h>
#include "header.h"

/******************************************/
/*  Functions to read and write header    */
/*    data for analysis files             */
/*                                        */
/*  Rob Maher    April 1987               */
/*                                        */
/*    use:                                */
/*           HEADER headptr;              */
/*              int fd;                   */
/*                                        */
/*        fd=open("file",0);              */
/*          --or--                        */
/*        fd=creat("file",0644);          */
/*                                        */
/*     rdat(fd,&headptr);                 */
/*     wdat(fd,&headptr);                 */
/*                                        */
/*Changes:                                */ 
/*  08/12/02 mert bay insert byteswap4()  */
/*                    for little endian   */
/*                    machines.(eg., PC)  */
/******************************************/


/********************************************/
/** function to acquire string from a file **/
/** -->returns pointer to allocated buffer **/
/********************************************/

 char *sgstring(fd)
     int fd;
  {
   int i;
   char buf;
   char *string;

   /** count chars in string **/
    i=0;
    do{
      
      i++;
      read(fd,&buf,sizeof(char));
      }while(buf!=0);
    
    /* create string buffer */
    string=(char *)calloc(i,sizeof(char));

    /* seek to head of string */
    lseek(fd,-1*i,1);
    read(fd,string,i*sizeof(char));
    return(string);
  }
  /*---*/


 /* function to fill header structure with data from file */

 void rdat(fd,head_data)
    int fd;         /*  file descriptor ( from open() )  */
    HEADER *head_data;
     
  {
  int i, *data;                                              /*mb 08/12/02*/
  head_data->performer= sgstring(fd);
  head_data->instrument= sgstring(fd);
  head_data->date= sgstring(fd);
  head_data->pitch= sgstring(fd);
  head_data->dyn= sgstring(fd);
  head_data->vibra= sgstring(fd);
  head_data->part= sgstring(fd);
  head_data->type= sgstring(fd);
  head_data->comments= sgstring(fd);
  head_data->andate= sgstring(fd);
  read(fd,&head_data->interpval,sizeof(float));
  read(fd,&head_data->sr,sizeof(float));
  read(fd,&head_data->tl,sizeof(float));
  read(fd,&head_data->smax,sizeof(float));
  read(fd,&head_data->fa,sizeof(float));
  read(fd,&head_data->dt,sizeof(float));
  read(fd,&head_data->fftlen,sizeof(int));
  read(fd,&head_data->nhar,sizeof(int));
  read(fd,&head_data->nchans,sizeof(int));
  read(fd,&head_data->npts,sizeof(int));
  if(byte_order)                                               /*mb 08/12/02*/
  {
  data=(int *)head_data;           /*if the machine byte order is little endian*/ 
  for(i=10;i<20;i++)               /*this loop swaps the bytes  mb 08/12/02*/
    byteswap4(data+i);                
  }                                                         
   if((head_data->npts)<=0 ){
     printf("File read error in rdat()\n");
     fflush(stdout);
     exit(1);
    }

  }
  /*---*/

 /** function to write header structure data to file **/

 void wdat(fd,head_data)
   int fd;
   HEADER *head_data;
  {
  /* include null terminators */
  write(fd,head_data->performer,1+strlen(head_data->performer));
  write(fd,head_data->instrument,1+strlen(head_data->instrument));
  write(fd,head_data->date,1+strlen(head_data->date));
  write(fd,head_data->pitch,1+strlen(head_data->pitch));
  write(fd,head_data->dyn,1+strlen(head_data->dyn));
  write(fd,head_data->vibra,1+strlen(head_data->vibra));
  write(fd,head_data->part,1+strlen(head_data->part));
  write(fd,head_data->type,1+strlen(head_data->type));
  write(fd,head_data->comments,1+strlen(head_data->comments));
  write(fd,head_data->andate,1+strlen(head_data->andate));
  write(fd,&head_data->interpval,sizeof(float));
  write(fd,&head_data->sr,sizeof(float));
  write(fd,&head_data->tl,sizeof(float));
  write(fd,&head_data->smax,sizeof(float));
  write(fd,&head_data->fa,sizeof(float));
  write(fd,&head_data->dt,sizeof(float));
  write(fd,&head_data->fftlen,sizeof(int));
  write(fd,&head_data->nhar,sizeof(int));
  write(fd,&head_data->nchans,sizeof(int));
  if( write(fd,&head_data->npts,sizeof(int)) <=0 ){
     printf("File write error in wdat()\n");
     exit(1);
    }
 }

/************* END OF FILE ****************/



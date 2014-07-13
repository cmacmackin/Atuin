/* posixwrapper.c - wrapper for posix functions, called from fortranposix.f90 */

/*     This file is a part of the fortranposix library. This makes system calls on Linux and Unix like systems available to Fortran programs. */
/*     Copyright (C) 2004  Madhusudan Singh */

/*     This library is free software; you can redistribute it and/or */
/*     modify it under the terms of the GNU Lesser General Public */
/*     License as published by the Free Software Foundation; either */
/*     version 2.1 of the License, or (at your option) any later version. */

/*     This library is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/*     Lesser General Public License for more details. */

/*     You should have received a copy of the GNU Lesser General Public */
/*     License along with this library; if not, write to the Free Software */
/*     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

/*     I would strongly welcome bug reports, feature enhancement requests, compliments, donations and job offers :) My email address : msc@ieee.org */


#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<sys/types.h>
#include<sys/wait.h>
#include<sys/stat.h>
#include<errno.h>
#include<pwd.h>

#ifdef __linux__
#include<linux/limits.h>
#endif

#ifdef __sun__
#include<limits.h>
#endif

FILE* fileptr[100];

char *appendtostring(char *str,int len)
{
  char *res; /* C arrays are from 0:len-1 */
  if((res=(char*)malloc(len+1)))
    {
      strncpy(res,str,len);
      res[len]='\0';
    }
  return res;
}

void fortrangetpipebufsize_(int *pipesize)
{
  *pipesize=PIPE_BUF;
}

void fortranpopen_(int *fd, char *cmd, char *typ,  int cmdlength, int typlength)
{
  char *command, *type;
  FILE* cmdpipe;
  extern int errno;

  if(!(command=appendtostring(cmd,cmdlength)))
    {
      perror("Failed : appendtostring (command) in posixwrapper:fortranpopen");
      exit(1);
    }

  if(!(type=appendtostring(typ,typlength)))
    {
      perror("Failed : appendtostring (type) in posixwrapper:fortranpopen");
      exit(1);
    }

  printf("Executing %s (%d) with option %s (%d) ... \n",command,strlen(command),type,strlen(type));

  errno=0;
  cmdpipe=popen(command,type);
  if(cmdpipe==NULL)
    {
      perror("Failed : popen in posixwrapper:fortranpopen");
      exit(1);
    }

  if((*fd=fileno(cmdpipe))==-1)
    {
      perror("Failed : fileno in posixwrapper:fortranpopen");
      exit(1);
    }

  fileptr[*fd]=cmdpipe;
  printf("Opened pipe for GNUPlot-Fortran95 IPC : file descriptor %d\n",*fd);
}

void fortranpclose_(int *fd, int *stat)
{
  int status;

  if((status=pclose(fileptr[*fd]))==-1)
    {
      perror("Failed : pclose in posixwrapper:fortranpclose");
      exit(1);
    }
  *stat=status;
}

void fortranaccess_(char *pathnm, int *mode, int *status, int pathnmlength)
{
  char *pathname;

  if(!(pathname=appendtostring(pathnm,pathnmlength)))
    {
      perror("Failed : appendtostring (pathname) in posixwrapper:fortranaccess");
      exit(1);
    }

  switch (*mode)
    {
      case 4 : // Check for read permissions
	*status=access(pathname,R_OK);
	break;
      case 2 : // Check for write permissions
	*status=access(pathname,W_OK);
	break;
      case 1 : // Check for execute permissions
	*status=access(pathname,X_OK);
	break;
      case 0 : // Check for existence
	*status=access(pathname,F_OK);
	break;
      default : // Invalid input
	*status=10;
	break;
    }
}

void fortranfflush_(int *fd, int *stat)
{
  extern int errno;

  errno=0;
  *stat=-1;
  if(fflush(fileptr[*fd])==EOF)
    {
      perror("Failed : fflush in posixwrapper:fortranfflush");
      exit(1);
    }
  *stat=0;
}

void fortrangetcwd_(char *pathnm, int *stat, int *lenstr, int pathnmlength)
{
  char *buf="getcwd";

  if(!(buf=getcwd(pathnm,pathnmlength)))
    {
      perror("Failed : getcwd in posixwrapper:fortrangetcwd");
      *stat=-1;
      exit(1);
    }
  *stat=0;
  *lenstr=strlen(pathnm);
}

void fortrangethostname_(char *hostnmstring, int *stat, int *lenstr, int hostnmlength)
{

  if((*stat=gethostname(hostnmstring,hostnmlength))==-1)
    {
      perror("Failed : gethostname in posixwrapper:fortrangethostname");
      exit(1);
    }
  *lenstr=strlen(hostnmstring);
  *stat=0;
}

void fortrangetlogin_(char *usernm, int *stat, int *lenstr, int usernmlength)
{
  char *username="getpwuid(geteuid)";
  struct passwd *p1;
  uid_t uid;

  *stat=-1;

  if((uid=getuid())==-1)
    {
      perror("Failed : getuid in posixwrapper:fortrangetlogin");
      exit(1);
    }

  if((p1=getpwuid(uid))==NULL)
    {
      perror("Failed : getpwuid in posixwrapper:fortranlogin");
      exit(1);
    }
  
  username=p1->pw_name;
  *lenstr=strlen(username);
  strncpy(usernm,username,*lenstr);
  *stat=0;
}

void fortrangetenv_(char *varnm, char *envnm, int *lenstr, int *status, int varnmlength,int envnmlength)
{
  char *varname="getenv";
  char *envname;
  
  if(!(varname=appendtostring(varnm,varnmlength)))
    {
      perror("Failed : appendtostring (varname) in posixwrapper:fortrangetenv");
      exit(1);
    }

  if(!(envname=getenv(varname)))
    {
      perror("Failed : getenv in posixwrapper:fortrangetenv");
      exit(1);
    }
  
  *lenstr=strlen(envname);
  strncpy(envnm,envname,*lenstr);
  *status=0;
}

void fortranchdir_(char *pathnm, int *stat, int pathnmlength)
{
  int status;
  char *pathname;

  if(!(pathname=appendtostring(pathnm,pathnmlength)))
    {
      perror("Failed : appendtostring (pathname) in posixwrapper:fortranchdir");
      exit(1);
    }
  
  if((status=chdir(pathname))==-1)
    {
      perror("Failed : chdir in posixwrapper:fortranchdir");
      exit(1);
    }
  *stat=status;
}

void fortranfputs_(char* buf, int *fd, int* status, int buflength)
{
  int stat;
  char *buffer="fputs";
  extern int errno;

  *status=-1;

  if(!(buffer=appendtostring(buf,buflength)))
    {
      perror("Failed : appendtostring (buf) in posixwrapper:fortranfputs");
      exit(1);
    }

  errno=0;

  if(fputs(buffer,fileptr[*fd])==EOF)
    {
      perror("Failed : fputs in posixwrapper:fortranfputs");
      exit(1);
    }

  *status=0;
}

void fortranfgets_(char *buf, int *fd, int *lenstr, int *status, int buflength)
{
  *status=-1;
  
  if(fgets(buf,buflength,fileptr[*fd])==NULL)
    {
      perror("Failed : fgets in posixwrapper:fortranfgets");
      exit(1);
    }
 
  *lenstr=strlen(buf);
  *status=0;
}

void fortranmkdir_(char *buf, char *mode, int *status,int buflength, int modelength)
{
  char *buffer;
  char *modestring;

  *status=-1;

  if(!(buffer=appendtostring(buf,buflength)))
    {
      perror("Failed : appendtostring (buf) in posixwrapper:fortranmkdir");
      exit(1);
    }

  if(!(modestring=appendtostring(mode,modelength)))
    {
      perror("Failed : appendtostring (mode) in posixwrapper:fortranmkdir");
      exit(1);
    }

  if(strcmp(modestring,"400")==0)
    *status=mkdir(buffer,S_IRUSR);
  else if(strcmp(modestring,"200")==0)
    *status=mkdir(buffer,S_IWUSR);
  else if(strcmp(modestring,"100")==0)
    *status=mkdir(buffer,S_IXUSR);
  else if(strcmp(modestring,"700")==0)
    *status=mkdir(buffer,S_IRWXU);
  else if(strcmp(modestring,"040")==0)
    *status=mkdir(buffer,S_IRGRP);
  else if(strcmp(modestring,"020")==0)
    *status=mkdir(buffer,S_IWGRP);
  else if(strcmp(modestring,"010")==0)
    *status=mkdir(buffer,S_IXGRP);
  else if(strcmp(modestring,"070")==0)
    *status=mkdir(buffer,S_IRWXG);
  else if(strcmp(modestring,"004")==0)
    *status=mkdir(buffer,S_IROTH);
  else if(strcmp(modestring,"002")==0)
    *status=mkdir(buffer,S_IWOTH);
  else if(strcmp(modestring,"001")==0)
    *status=mkdir(buffer,S_IXOTH);
  else if(strcmp(modestring,"007")==0)
    *status=mkdir(buffer,S_IRWXO);
  else
    {
      *status=-1;
      printf("Invalid octal string specified : mkdir in posixwrapper:fortranmkdir\n");
    }
      
  if(*status==-1)
    {
      perror("Failed : mkdir in posixwrapper:fortranmkdir");
      exit(1);
    }

  *status=0;
}


void fortranrmdir_(char *buf, int *status,int buflength)
{
  char *buffer;
  extern int errno;

  *status=-1;
  
  if(!(buffer=appendtostring(buf,buflength)))
    {
      perror("Failed : appendtostring (buf) in posixwrapper:fortranrmdir");
      exit(1);
    }

  errno=0;
  *status=rmdir(buffer);
    
  if(*status==-1)
    {
      printf("%d %s\n",errno,strerror(errno));
      perror("Failed : rmdir in posixwrapper:fortranrmdir");
      exit(1);
    }
  
  *status=0;
}

void fortrangetpid_(int *pidnow)
{
  pid_t pid;

  pid=getpid();
  *pidnow=pid;
}

void fortrangetppid_(int *pidnow)
{
  pid_t ppid;

  ppid=getppid();
  *pidnow=ppid;
}

void fortrangetchar_(char* readchar, int* status)
{
  *status=-1;

  if((*readchar=getchar())==EOF)
    {
      perror("Failed : getchar in posixwrapper:fortrangetchar");
    }

  *status=0;  
}


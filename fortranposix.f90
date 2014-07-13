!!$ fortranposix.f90 - Fortran POSIX routines that call wrapper routines in posixwrapper.c
!!$
!!$     This file is a part of the fortranposix library. This makes system calls on Linux and Unix like systems available to Fortran programs.
!!$     Copyright (C) 2004  Madhusudan Singh 
!!$
!!$     This library is free software; you can redistribute it and/or 
!!$     modify it under the terms of the GNU Lesser General Public 
!!$     License as published by the Free Software Foundation; either 
!!$     version 2.1 of the License, or (at your option) any later version.
!!$
!!$     This library is distributed in the hope that it will be useful,
!!$     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!!$     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
!!$     Lesser General Public License for more details. 
!!$
!!$     You should have received a copy of the GNU Lesser General Public 
!!$     License along with this library; if not, write to the Free Software 
!!$     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
!!$
!!$       I would strongly welcome bug reports, feature enhancement requests, compliments, donations and job offers :) My email address : msc@ieee.org
!!$
!!$ ADDITIONAL NOTES :
!!$
!!$ Uses calls to underlying C routines. Keep in mind the usual bromides about Fortran 95 / C interop and underscores, etc.
!!$




!!$          Function : fortran_pipebufsize
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. pipe_buf_size : integer : size of the buffer that can be atomically written to a pipe
!!$
!!$          Purpose : Finds the size of the buffer that can be atomically written to a pipe
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_pipebufsize() result(pipe_buf_size)
!!$                   end function fortran_pipebufsize
!!$                 end interface

          integer function fortran_pipebufsize() result(pipe_buf_size)

            call fortrangetpipebufsize(pipe_buf_size)
            
          end function fortran_pipebufsize

!!$          Function : chopnull
!!$
!!$          Inputs :
!!$
!!$                 1. str : character string : String to be processed
!!$                 2. len : integer : length of the string without any trailing blanks
!!$                 3. maxlen : integer : maximum length allocated for the string
!!$
!!$          Outputs :
!!$
!!$                 1. str : character string : Processed string
!!$
!!$          Purpose : Chops off the null character and the spaces from a string upon return from F95/C interop
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Interface :
!!$
!!$                 None

          
          subroutine chopnull(str,len,maxlen)
            
            character(len=*), intent(inout) :: str
            integer, intent(in) :: len,maxlen
            
            str(len+1:maxlen)=''
          end subroutine chopnull

!!$          Function : fortran_popen
!!$
!!$          Inputs :
!!$
!!$                 1. cmd : character string : Command to be piped
!!$                 2. len : character string (len = 1) : type of pipe to be opened
!!$
!!$          Outputs :
!!$
!!$                 1. fd : integer : file descriptor for the pipe
!!$
!!$          Purpose : Mimics the popen() system call
!!$
!!$          Diagnostics (incomplete) :
!!$
!!$                 1. fd = -1 : failed to initiate the popen() call
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_popen(cmd,type) result(fd)
!!$                     character(len=*), intent(in) :: cmd
!!$                     character(len=*), intent(in) :: type
!!$                   end function fortran_popen
!!$                 end interface
          
          integer function fortran_popen(cmd,type) result(fd)
            character(len=*), intent(in) :: cmd
            character(len=*), intent(in) :: type
            
            fd=-1
            
            if(len_trim(type).ne.1) then
               print*, 'Failed : fortran_popen : Invalid specification for pipe type'
               return
            end if
            
            call fortranpopen(fd,cmd,type)
          end function fortran_popen

!!$          Function : fortran_pclose
!!$
!!$          Inputs :
!!$
!!$                 1. fd : integer : file descriptor for the pipe to be closed
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the pclose() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to close pipe
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_pclose(fd) result(status)
!!$                     integer, intent(in) :: fd
!!$                   end function fortran_pclose
!!$                 end interface
          
          integer function fortran_pclose(fd) result(status)
            integer, intent(in) :: fd
            
            status=-1
            
            call fortranpclose(fd,status)
            
          end function fortran_pclose

!!$          Function : fortran_access
!!$
!!$          Inputs :
!!$
!!$                 1. pathname : character string : pathname for the file to be checked
!!$                 2. mode : integer : 
!!$                                     4 - check for read permissions
!!$                                     2 - check for write permissions
!!$                                     1 - check for execute permissions
!!$                                     0 - check for existence
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the access() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to check for permissions
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_access(pathname,mode) result(status)
!!$                     character(len=*), intent(in) :: pathname
!!$                     integer, intent(in) :: mode
!!$                   end function fortran_access
!!$                 end interface
          
          integer function fortran_access(pathname,mode) result(status)
            character(len=*), intent(in) :: pathname
            integer, intent(in) :: mode
            
            status=-1
            call fortranaccess(pathname,mode,status)
          end function fortran_access

!!$          Function : fortran_fflush
!!$
!!$          Inputs :
!!$
!!$                 1. fd : integer : file descriptor for the pipe to be flushed
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the fflush() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to flush the pipe
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_fflush(fd) result(status)
!!$                     integer, intent(in) :: fd
!!$                   end function fortran_fflush
!!$                 end interface
 
          integer function fortran_fflush(fd) result(status)
            integer, intent(in) :: fd
            
            status=-1
            call fortranfflush(fd,status)
          end function fortran_fflush

!!$          Function : fortran_getcwd
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. buf : character string : contains the CWD upon successful exit
!!$                 2. status : integer : error flag
!!$
!!$          Purpose : Mimics the getcwd() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to get current working directory
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_getcwd(buf) result(status)
!!$                     character(len=*), intent(out) :: buf
!!$                   end function fortran_getcwd
!!$                 end interface
          
          integer function fortran_getcwd(buf) result(status)
            character(len=*), intent(out) :: buf
            
            integer :: bufsize,maxlen
            
            maxlen=len(buf)
            status=-1
            
            call fortrangetcwd(buf,status,bufsize)
            if(status.eq.0) call chopnull(buf,bufsize,maxlen)
          end function fortran_getcwd

!!$          Function : fortran_gethostname
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. buf : character string : contains the hostname upon successful exit
!!$                 2. status : integer : error flag
!!$
!!$          Purpose : Mimics the gethostname() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to get hostname
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_gethostname(hostname) result(status)
!!$                     character(len=*), intent(out) :: hostname
!!$                   end function fortran_gethostname
!!$                 end interface
          
          integer function fortran_gethostname(hostname) result(status)
            character(len=*), intent(out) :: hostname
            
            integer :: hostnamelength,maxlen
            
            maxlen=len(hostname)
            status=-1
            
            call fortrangethostname(hostname,status,hostnamelength)
            if(status.eq.0) call chopnull(hostname,hostnamelength,maxlen)
          end function fortran_gethostname

!!$          Function : fortran_getlogin
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. buf : character string : contains the login of the current user upon successful exit
!!$                 2. status : integer : error flag
!!$
!!$          Purpose : Mimics the getlogin() system call safely through the getpwuid(getuid) mechanism
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to get current user's login name
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_getlogin(loginstring) result(status)
!!$                     character(len=*), intent(out) :: loginstring
!!$                   end function fortran_getlogin
!!$                 end interface
         
          integer function fortran_getlogin(loginstring) result(status)
            character(len=*), intent(out) :: loginstring
            
            integer :: loginstringlength,maxlen
            
            maxlen=len(loginstring)
            status=-1
            
            call fortrangetlogin(loginstring,status,loginstringlength)
            if(status.eq.0) call chopnull(loginstring,loginstringlength,maxlen)
          end function fortran_getlogin

!!$          Function : fortran_getenv
!!$
!!$          Inputs :
!!$
!!$                 1. varnamestring : character string containing the name of the environment variable
!!$
!!$          Outputs :
!!$
!!$                 1. environmentstring : character string : contains the environment variable upon successful exit
!!$                 2. status : integer : error flag
!!$
!!$          Purpose : Mimics the getenv() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to get the value of the requested environment variable
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_getenv(varnamerequested,environmentstring) result(status)
!!$                     character(len=*), intent(in) :: varnamerequested
!!$                     character(len=*), intent(out) :: environmentstring
!!$                   end function fortran_getenv
!!$                 end interface

          
          integer function fortran_getenv(varnamerequested,environmentstring) result(status)
            character(len=*), intent(in) :: varnamerequested
            character(len=*), intent(out) :: environmentstring
            
            integer :: environmentstringlength,maxlen
            
            maxlen=len(environmentstring)
            status=-1
            
            call fortrangetenv(varnamerequested,environmentstring,environmentstringlength,status)
            if(status.eq.0) call chopnull(environmentstring,environmentstringlength,maxlen)
          end function fortran_getenv

!!$          Function : fortran_chdir
!!$
!!$          Inputs :
!!$
!!$                 1. pathname : character string : directory to change to
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the chdir() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to change directory
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_chdir(pathname) result(status)
!!$                     character(len=*), intent(in) :: pathname
!!$                   end function fortran_chdir
!!$                 end interface
          
          integer function fortran_chdir(pathname) result(status)
            character(len=*), intent(in) :: pathname
            
            status=-1
            
            call fortranchdir(pathname,status)
          end function fortran_chdir

!!$          Function : fortran_fputs
!!$
!!$          Inputs :
!!$
!!$                 1. buf : character string : buffer to be written to fd
!!$                 2. fd : integer : file descriptor
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the fputs() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to write the string
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_fputs(buf,fd) result(status)
!!$                     character(len=*), intent(in) :: buf
!!$                     integer, intent(in) :: fd
!!$                   end function fortran_fputs
!!$                 end interface

          
          integer function fortran_fputs(buf,fd) result(status)
            character(len=*), intent(in) :: buf
            integer, intent(in) :: fd
            
            status=-1
            
            call fortranfputs(buf,fd,status)
          end function fortran_fputs

!!$          Function : fortran_fgets
!!$
!!$          Inputs :
!!$
!!$                 1. fd : integer : file descriptor
!!$
!!$          Outputs :
!!$
!!$                 1. buf : character string : buffer to be read from fd
!!$                 2. status : integer : error flag
!!$
!!$          Purpose : Mimics the fgets() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to read the string
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_fgets(buf,fd) result(status)
!!$                     character(len=*), intent(out) :: buf
!!$                     integer, intent(in) :: fd
!!$                   end function fortran_fgets
!!$                 end interface
          
          integer function fortran_fgets(buf,fd) result(status)
            character(len=*), intent(out) :: buf
            integer, intent(in) :: fd
            
            integer :: buflength,maxlen
            
            maxlen=len(buf)
            status=-1
            
            call fortranfgets(buf,fd,buflength,status)
            if(status.eq.0) call chopnull(buf,buflength,maxlen)
          end function fortran_fgets

!!$          Function : fortran_mkdir
!!$
!!$          Inputs :
!!$
!!$                 1. dirname : character string : name of the directory to make
!!$                 2. mode : character string : string containing the usual chmod style permission bits (octal)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the mkdir() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to make directory
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_mkdir(dirname,mode) result(status)
!!$                     character(len=*), intent(in) :: dirname,mode
!!$                   end function fortran_mkdir
!!$                 end interface

          
          integer function fortran_mkdir(dirname,mode) result(status)
            character(len=*), intent(in) :: dirname,mode

            status=-1

            call fortranmkdir(dirname,mode,status)
          end function fortran_mkdir

!!$          Function : fortran_rmdir
!!$
!!$          Inputs :
!!$
!!$                 1. dirname : character string : name of the directory to remove
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer : error flag
!!$
!!$          Purpose : Mimics the rmdir() system call
!!$
!!$          Diagnostics :
!!$
!!$                 1. status = -1 : failed to remove directory
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_rmdir(dirname) result(status)
!!$                     character(len=*), intent(in) :: dirname
!!$                   end function fortran_rmdir
!!$                 end interface


          integer function fortran_rmdir(dirname) result(status)
            character(len=*), intent(in) :: dirname

            status=-1

            call fortranrmdir(dirname,status)
          end function fortran_rmdir

!!$          Function : fortran_getpid
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. pid : Process ID of the calling process
!!$
!!$          Purpose : Mimics the getpid() system call
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_getpid() result(pid)
!!$                   end function fortran_getpid
!!$                 end interface

          integer function fortran_getpid() result(pid)

            call fortrangetpid(pid)
          end function fortran_getpid

!!$          Function : fortran_getppid
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. ppid : Parent's Process ID of the calling process
!!$
!!$          Purpose : Mimics the getppid() system call
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_getppid() result(ppid)
!!$                   end function fortran_getppid
!!$                 end interface

          integer function fortran_getppid() result(ppid)

            call fortrangetpid(ppid)
          end function fortran_getppid

!!$          Function : fortran_getchar
!!$
!!$          Inputs :
!!$
!!$                 None
!!$
!!$          Outputs :
!!$
!!$                 1. readchar : character
!!$
!!$          Purpose : Reads a character from stdin
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Interface :
!!$
!!$                 interface
!!$                   integer function fortran_getchar(readchar) result(status)
!!$                     character(len=1), intent(out) :: readchar
!!$                   end function fortran_getchar
!!$                 end interface

          integer function fortran_getchar(readchar) result(status)
            character(len=1), intent(out) :: readchar
            
            call fortrangetchar(readchar,status)

          end function fortran_getchar


!!$ gnuplot_fortran95  - Fortran 95 interface to some gnuplot functions
!!$ Date - 1/27/2004 Version 0.1
!!$ Copyright(C) 2004 Madhusudan Singh
!!$
!!$    This library is free software; you can redistribute it and/or
!!$    modify it under the terms of the GNU Lesser General Public
!!$    License as published by the Free Software Foundation; either
!!$    version 2.1 of the License, or (at your option) any later version.
!!$
!!$    This library is distributed in the hope that it will be useful,
!!$    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!$    Lesser General Public License for more details.
!!$
!!$    You should have received a copy of the GNU Lesser General Public
!!$    License along with this library; if not, write to the Free Software
!!$    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!!$
!!$
!!$ This library is released under LGPL. The full text of LGPL may be found
!!$ at http://www.fsf.org/licenses/lgpl.txt
!!$
!!$ In particular, the terms of the GNU LGPL imply that you CAN use this library with
!!$ proprietary code with some restrictions (you can link, but not use the code directly). Please consult the above URL for exact details.
!!$
!!$ ADDITIONAL NOTES :
!!$
!!$ This will work for any system that :
!!$
!!$ 1. Provides popen, pclose, and other such functions.
!!$ 2. C object files append an underscore to the end of C library functions.
!!$ 3. You will need a copy of the fortranposix library ((C) Madhusudan Singh).
!!$
!!$ Please do NOT request a Windows targeted implementation as I do not use M$ Windoze or any of its variants. And probably never will.
!!$
!!$
!!$ RATIONALE :
!!$
!!$ I have seen people who prefer a simple plotting system capable of enormous
!!$ sophistication use gnuplot. I have seen people using fortran 90/95 for
!!$ serious number crunching. What I found missing to my surprise, was something that
!!$ linked the two. Granted that other libraries like pgplot, psplot, dislin, etc exist,
!!$ but there is something nice about the simplicity of gnuplot which attracts so many users.
!!$
!!$ Gnuplot is a freely available, command-driven graphical display tool for
!!$ Unix. It compiles and works quite well on a number of Unix flavours as
!!$ well as other operating systems. The following module enables sending
!!$ display requests to gnuplot through simple f90/95 calls.
!!$
!!$ This module is released under the GNU Lesser General Public License (LGPL), the
!!$ terms of which are included here for your reference.
!!$
!!$ This has been tested with GNUPlot 3.7.3.
!!$
!!$ GENERAL DIAGNOSTICS :
!!$
!!$ status = 0 <-----> SUCCESS
!!$ status <0 or NULL <-----> ERROR
!!$ status > 0 <-----> WARNING
!!$
!!$ DISCLAIMER :
!!$
!!$ Though it has been stated explicitly in the terms of the License, I would like to reiterate the following :
!!$
!!$ I disclaim ALL responsibility, express or implied, in using the library, whether any damages, losses, etc.
!!$ are deemed to be directly caused or indirectly contributed in any fashion. Please
!!$ understand that you use this library COMPLETELY at your OWN risk.
!!$
!!$ DOCUMENTATION :
!!$
!!$ This library is documented in the source itself. Please consult it before using the subroutines and functions.
!!$
!!$ CONTACT INFORMATION :
!!$
!!$ I would appreciate receiving some feedback upon your use of this code in
!!$ your programs and would request that appropriate acknowledgement be made
!!$ whenever you use this code in your programs.
!!$
!!$ If you feel that you found this library to be useful / have found a bug /
!!$ / want to discuss a feature enhancement, drop me a line at msingh@myrealbox.com.
!!$
!!$ Madhusudan Singh, Ann Arbor, early spring, 2004
!!$

      module datatypes

        implicit none

         ! These are numerical recipes style type declarations

        integer, parameter :: i4b = SELECTED_INT_KIND(9)
        integer, parameter :: dp =  KIND(1.0d0)
        integer, parameter :: lgc = KIND(.true.)
!!$        integer, parameter :: dpc = KIND((1.0d0,1.0d0))
      end module datatypes

      module gnuplot_module_data
        use datatypes
        implicit none

        private
        public :: i4b,dp,lgc,GP_CMD_SIZE,MAXLENGTH,BAD_SIGNAL,DECIMALASCIICODEFORNEWLINE,NUMPLOTSTYLES,&
             FILEPERMISSIONS,DATAFILELUNITOFFSET,TOPLEVELDIR,WORKING_DIRECTORY,SHELLPATH,VALIDPLOTSTYLE,&
             PATH_MAXNAMESZ,PATH_LENGTH,GP_TITLE_SIZE,GP_EQ_SIZE,PLOTSTYLELENGTH,gnuplot_ctrl,PI_D,plotdetails,&
             GNUPLOT_EXECUTABLE,GP_MAX_TMP_FILES,GNUPLOT_SHOWDEBUG,ACKTAG,GNUPLOT_SHOWWARNINGS

        integer(i4b), parameter :: GP_MAX_TMP_FILES = 64                             ! Maximum number of simultaenous temporary files
        integer(i4b), parameter :: GP_TMP_NAME_SIZE = 512                            ! Maximum size of a temporary file name
        integer(i4b), parameter :: GP_CMD_SIZE = 2048
        integer(i4b), parameter :: GP_TITLE_SIZE = 80
        integer(i4b), parameter :: GP_EQ_SIZE = 512
        integer(i4b), parameter :: PATH_MAXNAMESZ = 260
        integer(i4b), parameter :: PATH_LENGTH=4096
        integer(i4b), parameter :: MAXLENGTH = 200
        integer(i4b), parameter :: BAD_SIGNAL=-9999
        integer(i4b), parameter :: DECIMALASCIICODEFORNEWLINE=10 ! Consult http://www.asciitable.com
        integer(i4b), parameter :: NUMPLOTSTYLES=13
        integer(i4b), parameter :: FILEPERMISSIONS=500
        integer(i4b), parameter :: DATAFILELUNITOFFSET=99
        integer(i4b), parameter :: PLOTSTYLELENGTH=12

        real(dp), parameter :: PI_D = 3.141592653589793238462643383279502884197_dp

        character(len=1), parameter :: TOPLEVELDIR='/'
        character(len=2), parameter :: WORKING_DIRECTORY='./'
        character(len=7), parameter :: SHELLPATH='/bin/sh'
        character(len=7), parameter :: GNUPLOT_EXECUTABLE='gnuplot'
        character(len=4), parameter :: ACKTAG='1024'
        character(len=PLOTSTYLELENGTH), dimension(NUMPLOTSTYLES), parameter :: VALIDPLOTSTYLE  = (/'lines       ',&
             'points      ', 'linespoints ','impulses    ','dots        ','steps       ','errorbars   ',&
             'boxes       ','boxerrorbars','financebars ','candlesticks','vector      ','pm3d        '/)

        type :: gnuplot_ctrl
           integer(i4b) :: fd  ! file descriptor for the pipe
           integer(i4b) :: nplots  ! number of currently active plots
           integer(i4b) :: ntmp ! number of temporary files
           character(len=MAXLENGTH) :: plotstyle  ! current plotting style
           character(len=MAXLENGTH) :: linestyle ! current line style
           character(len=MAXLENGTH) :: linetype ! current line type
           character(len=MAXLENGTH) :: linewidth ! current line width
           character(len=MAXLENGTH) :: pointtype ! current pointtype
           character(len=MAXLENGTH) :: pointsize ! current pointsize
           integer(i4b) :: hardcopystatus
           character(len=MAXLENGTH) :: hardcopyformat
           character(len=MAXLENGTH) :: hardcopyfilename
           integer(i4b) :: hardcopyfilenamelength
           character(len=MAXLENGTH) :: tmp_dir
           character(len=MAXLENGTH) :: xoption
           character(len=MAXLENGTH), dimension(GP_MAX_TMP_FILES) :: DATAFILENAME
           integer(i4b), dimension(GP_MAX_TMP_FILES) :: DATAFILELUNIT
           logical :: pubflag
        end type gnuplot_ctrl

        type(gnuplot_ctrl), target :: plotdetails

        integer(i4b) :: gnuplot_pid ! Process ID for the spawned gnuplot process

        logical(lgc) :: GNUPLOT_SHOWDEBUG=.false. ! Used to echo commands passed to gnuplot_cmd
        logical(lgc) :: GNUPLOT_SHOWWARNINGS=.false.
      end module gnuplot_module_data

      module polar_cartesian
        use datatypes
        type polar
           real(dp) :: radial, theta
        end type polar

        interface assignment(=)
           module procedure cartesian_to_polar
           module procedure polar_to_cartesian
        end interface

        contains

          subroutine cartesian_to_polar(polar_out,cartesian_in)
            type(polar), intent(out) :: polar_out
            complex(dp), intent(in) :: cartesian_in

            polar_out%radial=abs(cartesian_in)
            if(cartesian_in.ne.0.0_dp) then
               polar_out%theta=atan2(aimag(cartesian_in),real(cartesian_in,dp))
            else
               polar_out%theta=0.0_dp
            end if
          end subroutine cartesian_to_polar

          subroutine polar_to_cartesian(cartesian_out,polar_in)
            complex(dp), intent(out) :: cartesian_out
            type(polar), intent(in) :: polar_in

            cartesian_out=cmplx(polar_in%radial*cos(polar_in%theta),polar_in%radial*sin(polar_in%theta),dp)
          end subroutine polar_to_cartesian

        end module polar_cartesian

      module gnuplot_module
        use gnuplot_module_data,only : gnuplot_ctrl,i4b,dp,lgc,GP_CMD_SIZE,MAXLENGTH,BAD_SIGNAL,&
             DECIMALASCIICODEFORNEWLINE,NUMPLOTSTYLES,FILEPERMISSIONS,DATAFILELUNITOFFSET,TOPLEVELDIR,WORKING_DIRECTORY,&
             SHELLPATH,VALIDPLOTSTYLE,PATH_MAXNAMESZ,GP_TITLE_SIZE,GP_EQ_SIZE,PLOTSTYLELENGTH,PATH_LENGTH,plotdetails,&
             GNUPLOT_EXECUTABLE,GP_MAX_TMP_FILES,GNUPLOT_SHOWDEBUG,ACKTAG,GNUPLOT_SHOWWARNINGS

        implicit none

        private

        public :: gnuplot_init,gnuplot_close,gnuplot_setstyle,gnuplot_setscale,gnuplot_setaxislabel,&
             gnuplot_settitle,gnuplot_reset,gnuplot_setrange,gnuplot_plot2d,gnuplot_plot3d,gnuplot_plotslope,&
             gnuplot_plotequation2d,gnuplot_hardcopy,gnuplot_resetsession,gnuplot_replot,gnuplot_setaxisformat,&
             gnuplot_setcoordinatesystem,gnuplot_setangleunits,gnuplot_plotcomplex,gnuplot_set,gnuplot_unset,&
             gnuplot_genpdf,gnuplot_plotpolarequation,gnuplot_plotequation3d,gnuplot_fitcurve,gnuplot_fitsurface,&
             gnuplot_pause,gnuplot_load,gnuplot_cmd,gnuplot_plot3dpanel,gnuplot_plot2dset

        contains

!!$          Subroutine : assert_error
!!$
!!$          Inputs :
!!$
!!$                 1. routinename : character string indicating which routine spawned the error
!!$                 2. errormessage : character string indicating the error message
!!$                 3. errornumber : integer (optional) indicating the kind of error
!!$
!!$          Outputs :
!!$
!!$                 None
!!$
!!$          Purpose : Prints error messages 
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Type : private
!!$
!!$          Interface : none

          subroutine assert_error(routinename,errormessage,errornumber)
            character(len=*), intent(in) :: routinename,errormessage
            integer(i4b), intent(in), optional :: errornumber

            write(*,'("GNUPlot-Fortran95 Interface Error Message :")')
            if(.not.present(errornumber)) then
               write(*,'(a24," : ",a100)') routinename,errormessage
            else
               write(*,'(a24," (error condition = ",i5,") : ",a70)') routinename,errornumber,errormessage
            end if

            write(*,'("GNUPlot-Fortran95 :",a24," will now exit ...")') routinename

          end subroutine assert_error

!!$          Subroutine : assert_warning
!!$
!!$          Inputs :
!!$
!!$                 1. routinename : character string indicating which routine spawned the warning
!!$                 2. warningmessage : character string indicating the warning message
!!$
!!$          Outputs :
!!$
!!$                 None
!!$
!!$          Purpose : Prints warning messages 
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Type : private
!!$
!!$          Interface : none

          subroutine assert_warning(routinename,warningmessage)
            character(len=*), intent(in) :: routinename,warningmessage

            if(GNUPLOT_SHOWWARNINGS) then
               write(*,'("GNUPlot-Fortran95 Interface Warning Message :")')
               write(*,'(a24," : ",a70)') routinename,warningmessage
               write(*,'("Program will now attempt to continue  ...")')
            end if

          end subroutine assert_warning

!!$          Function : fileexists
!!$
!!$          Inputs :
!!$
!!$                 1. filename : character string indicating the relative pathname to the file to be checked
!!$
!!$          Outputs :
!!$
!!$                 1. existence : logical indicating whether the requested file is present
!!$
!!$          Purpose : Checks for the existence or otherwise of an input file
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Type : private
!!$
!!$          Interface : none

          logical function fileexists(filename) result(existence)
            character(len=*), intent(in) :: filename

            existence=.false.
            inquire(file=filename,exist=existence)
          end function fileexists


!!$          Function : setacknowledgement
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. hcflag : integer indicating the position of the ackknowledgement label in non-publication plots
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets acknowledgement string
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = 22 : Failed to extract hostname
!!$                 2. status = 23 : Failed to extract login name
!!$                 3. status = 24 : Failed to set acknowledgement string
!!$
!!$          Type : private
!!$
!!$          Interface : none

          integer function setacknowledgement(ptr_gctrl,hcflag) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            integer(i4b), intent(in) :: hcflag

            integer(i4b), external :: fortran_gethostname,fortran_getlogin

            character(len=GP_CMD_SIZE) :: cmd
            character(len=8) :: datestring=''
            character(len=10) :: timestring=''
            character(len=5) :: zonestring=''
            character(len=MAXLENGTH) :: hostnamestring=''
            character(len=MAXLENGTH) :: loginstring=''
            character(len=13) :: ackposstring=''
            integer(i4b) :: ierror=0

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('setacknowledgement','called with a dissociated pointer',status)
               return
            end if

            call date_and_time(datestring,timestring,zonestring)

            ierror=fortran_gethostname(hostnamestring)

            if(ierror.ne.0) then
               status=22
               call assert_warning('setacknowledgement','failed to extract host name')
               return
            end if

            ierror=fortran_getlogin(loginstring)

            if(ierror.ne.0) then
               status=23
               call assert_warning('setacknowledgement','failed to extract login name')
               return
            end if

            if(hcflag.eq.0) then
               ackposstring='0.1,1.01,-1.0'
            else if (hcflag.eq.1) then
               ackposstring='1.01,0.1,-1.0'
            else if (hcflag.eq.2) then
               ackposstring='1.21,1.1,-0.8'
            end if

            cmd='set label '//trim(ACKTAG)//' "Created using GNUPLOT-FORTRAN95 interface (Madhusudan Singh) (C) 2004 : by '//&
                 trim(loginstring)//'@'//trim(hostnamestring)//'" at graph '//trim(ackposstring)//' rotate font "Helvetica,7"'
            ierror=gnuplot_cmd(ptr_gctrl,cmd)
            if(ierror.ne.0) then
               status=24
               call assert_warning('gnuplot_init','failed to set acknowledgement string')
               return
            end if

            status=0

          end function setacknowledgement

!!$          Function : registernewdatafile
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 1. filename : character string indicating the relative pathname to the file to be registered
!!$
!!$          Outputs :
!!$
!!$                 1. ok : logical value indicating if the requested file was successfully registered
!!$
!!$          Purpose : Registers a new temporary data file with the interface
!!$
!!$          Diagnostics :
!!$
!!$                 None
!!$
!!$          Type : private
!!$
!!$          Interface : none


          logical function registernewdatafile(ptr_gctrl,filename) result(ok)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: filename

            character(len=MAXLENGTH) :: filesuffix
            integer(i4b) :: suffixlength

            ok=.false.

            ptr_gctrl%ntmp=ptr_gctrl%ntmp+1

            suffixlength=int(log10(real(ptr_gctrl%ntmp)))+1
            write(filesuffix,'(i10)') ptr_gctrl%ntmp

            if(ptr_gctrl%ntmp.le.GP_MAX_TMP_FILES) then
               ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp)=trim(filename)//"."//trim(adjustl(filesuffix))
               ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp)=DATAFILELUNITOFFSET+ptr_gctrl%ntmp
               ok=.true.
            else
               call assert_error('registernewdatafile','too many temporary files')
            end if

          end function registernewdatafile

!!$          Subroutine : gnuplot_get_program_path
!!$
!!$          Inputs :
!!$
!!$                 1. pname : character string (name of the program - in this case, gnuplot)
!!$
!!$          Outputs :
!!$
!!$                 1. buf : character string (full path to the executable)
!!$                 2. status : integer (error flag)
!!$
!!$          Purpose : Retrieves the full path to the requested executable
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = -1 : Cannot get current working directory
!!$                 2. status = -2 : Cannot get path string
!!$                 3. status = -3 : Access denied (bad parameter passed)
!!$
!!$          Type : private
!!$
!!$          Interface : none


          subroutine gnuplot_get_program_path(pname,status)

            character(len=*), intent(in) :: pname
            integer(i4b), intent(out) :: status

            integer(i4b), external :: fortran_getenv,fortran_access,fortran_getcwd

            integer(i4b) :: len_path,ierror=0,enoent=0,enotdir=0,eacces=0,i=1,k=1,lenpathstring
            character(len=PATH_LENGTH) :: pathstring
            character(len=MAXLENGTH) :: cwd
            character(len=PATH_MAXNAMESZ) :: buf

            status=0

            status=fortran_getenv('PATH',pathstring)
            lenpathstring=len_trim(pathstring)

            if((pathstring.eq.'').or.(status.ne.0)) then
               status=-2
               call assert_error('gnuplot_get_program_path','path string not found',status)
               return
            end if

            ierror=fortran_getcwd(cwd)

            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_get_program_path','failed to extract current working directory',status)
               return
            end if

            buf=trim(cwd)//"/"//pname
            len_path=len_trim(buf)

            ierror=fortran_access(trim(buf),1)
            if(ierror.eq.0) return ! Path found in current working directory

            buf=''
            k=1
            checkpath : do while (k.le.lenpathstring)
               i=index(pathstring(k:lenpathstring),':')
               buf=pathstring(k:k+i-2)//'/'//trim(pname)
               len_path=len_trim(buf)
               ierror=fortran_access(trim(buf),1)
               if(ierror.eq.0) exit checkpath
               k=k+i
            end do checkpath

            if(ierror.eq.0) return

            if(ierror.eq.-1) then
               status=-3
               call assert_error('gnuplot_get_program_path','Bad parameter passed to fortran_access',status)
            end if

            call assert_error('gnuplot_get_program_path','you do not have execute permissions for gnuplot !')


          end subroutine gnuplot_get_program_path

!!$          Function : gnuplot_setdefaults
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets various plotting paramters to their default values
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Failure
!!$
!!$          Type : private
!!$
!!$          Interface : none

          integer function gnuplot_setdefaults(ptr_gctrl) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setdefaults','called with a dissociated pointer',status)
               return
            end if

            ptr_gctrl%linestyle='1'
            ptr_gctrl%linetype='1'
            ptr_gctrl%linewidth='1'
            ptr_gctrl%pointtype='1'
            ptr_gctrl%pointsize='1'
            ptr_gctrl%ntmp=0
            ptr_gctrl%nplots=0
            ptr_gctrl%plotstyle='linespoints'
            ptr_gctrl%hardcopyformat='X11'
            status=0
          end function gnuplot_setdefaults


!!$          Function : gnuplot_init
!!$
!!$          Inputs :
!!$
!!$                 1. xoption : character string that passes any X11 options to the GNUPlot session, such as -persist,-background,-clear,-raise,-noraise,-tvtwm,-gray,or -mono.
!!$
!!$          Outputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctrl.
!!$
!!$          Purpose : Initiates a gnuplot session
!!$
!!$          Diagnostics :
!!$
!!$                 0. ptr_gctrl = associated : Success
!!$                 1. ptr_gctrl = NULL : Failure
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$             interface
!!$                function gnuplot_init(xoption) result(ptr_gctrl)
!!$                  character(len=*), intent(in) :: xoption
!!$                end function gnuplot_init
!!$             end interface

          function gnuplot_init(xoption) result(ptr_gctrl)
            character(len=*), intent(in) :: xoption

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            integer(i4b), external :: fortran_getenv,fortran_popen,fortran_gethostname,fortran_getlogin,fortran_mkdir

            character(len=MAXLENGTH) :: displaystring
            character(len=GP_CMD_SIZE) :: cmd
            integer(i4b) :: istat=0,writefd=0,ierror=0
            character(len=8) :: datestring=''
            character(len=10) :: timestring=''
            character(len=5) :: zonestring=''
            character(len=MAXLENGTH) :: hostnamestring=''
            character(len=MAXLENGTH) :: loginstring=''

            nullify(ptr_gctrl)
            if(.not.associated(ptr_gctrl,plotdetails)) then
               allocate(ptr_gctrl,STAT=istat)
               if(istat.ne.0) then
                  call assert_error('gnuplot_init','could not allocate pointer',istat)
               end if
               ptr_gctrl=>plotdetails
            end if

            istat=fortran_getenv('DISPLAY',displaystring)

            if((len_trim(displaystring).eq.0).or.(istat.ne.0)) then
               call assert_error('gnuplot_init','error identifying display',istat)
               nullify(ptr_gctrl)
               return
            end if

            write(*,'("gnuplot_init : current display identified as ",a6,"  {",i1,"}")') &
                 trim(displaystring),len_trim(displaystring)

            call gnuplot_get_program_path(GNUPLOT_EXECUTABLE,ierror)
            if(ierror.ne.0) then
               call assert_error('gnuplot_init','cannot find gnuplot in your $PATH',ierror)
               nullify(ptr_gctrl)
               return
            end if

            if(len_trim(xoption).ne.0) write(*,'("No checking of ",a12," as a valid X option will be done ...")') xoption
            ptr_gctrl%xoption=xoption

            writefd=fortran_popen(GNUPLOT_EXECUTABLE//" "//trim(xoption),'w')

            if(writefd.eq.-1) then
               call assert_error('gnuplot_init','failed to initiate IPC',writefd)
               nullify(ptr_gctrl)
               return ! ptr_gctrl is dissociated at this point
            end if

            ptr_gctrl%fd=writefd

            ierror=gnuplot_setdefaults(ptr_gctrl)

            if(ierror.ne.0) then
               call assert_error('gnuplot_init','failed to set plot defaults',ierror)
               nullify(ptr_gctrl)
               return
            end if

            call date_and_time(datestring,timestring,zonestring)

            ierror=fortran_gethostname(hostnamestring)

            if(ierror.ne.0) then
               call assert_error('gnuplot_init','failed to extract host name',ierror)
               nullify(ptr_gctrl)
               return
            end if

            ierror=fortran_getlogin(loginstring)

            if(ierror.ne.0) then
               call assert_error('gnuplot_init','failed to extract login name',ierror)
               nullify(ptr_gctrl)
               return
            end if

            ptr_gctrl%tmp_dir='.gnuplot-f95-'//trim(datestring)//trim(zonestring)//trim(timestring)//&
                 trim(hostnamestring)

            ierror=fortran_mkdir(trim(ptr_gctrl%tmp_dir),'700')

            if(ierror.ne.0) then
               call assert_error('gnuplot_init','failed to make temporary directory',ierror)
               nullify(ptr_gctrl)
               return
            end if

            ierror=setacknowledgement(ptr_gctrl,0)
            ptr_gctrl%hardcopystatus=1 ! Necessary because gnuplot 4.0 requires a set output
          end function gnuplot_init

!!$          Function : gnuplot_close
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Concludes a gnuplot session
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -8888 : PDF conversion requested and failed.
!!$                 3. status = -7777 : Cannot close communication pipe to gnuplot
!!$                 4. status = -1 : Could not remove temporary directory
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$             interface
!!$                integer function gnuplot_close(ptr_gctrl) result(status)
!!$                  type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$                end function gnuplot_close
!!$             end interface


          integer function gnuplot_close(ptr_gctrl) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            integer(i4b), external :: fortran_pclose,fortran_rmdir

            integer(i4b) :: i=0,ierror=0,istat=0
            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_close','called with a dissociated pointer',status)
               return
            end if

            cmd='quit'
            status=gnuplot_cmd(ptr_gctrl,cmd)

            ierror=fortran_pclose(ptr_gctrl%fd)
            if(ierror.ne.0) then
               status=-7777
               call assert_error('gnuplot_close','failed to close IPC with GNUPlot',status)
               return
            end if

!!$ This is called after issuing the quit command, should ensure that the plotting is complete before datafile is deleted, the second predicate ? just being paranoid :)
            ifplottingdonedeletedatafile : if((status.eq.0).and.(ptr_gctrl%ntmp.ge.1)) then
               deletealltmpfiles : do i=1,ptr_gctrl%ntmp
                  open(unit=ptr_gctrl%DATAFILELUNIT(i),file=ptr_gctrl%DATAFILENAME(i),status='old',iostat=istat)
                  if(istat.eq.0) close(ptr_gctrl%DATAFILELUNIT(i),status='delete',iostat=ierror)
!!$ Shows the superiority of FORTRAN over C - do not need to worry about deleting temporary files explicitly !
                  if(ierror.ne.0) then
                     status=22
                     call assert_warning('gnuplot_close','failed to delete temporary file '//trim(ptr_gctrl%DATAFILENAME(i)))
                  end if
               end do deletealltmpfiles
            end if ifplottingdonedeletedatafile

            ierror=fortran_rmdir(trim(WORKING_DIRECTORY//ptr_gctrl%tmp_dir))
            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_close','failed to remove temporary directory',status)
               return
            end if

            nullify(ptr_gctrl)

          end function gnuplot_close

!!$          Function : gnuplot_cmd
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. cmd : character string (command to be executed in the gnuplot session)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Executes a single line of commands in the gnuplot session via a pipe
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -8888 : Failed to write to pipe
!!$                 3. status = -7777 : Command too long for pipe buffer
!!$                 4. status = -6666 : Failed to flush the pipe
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$              interface
!!$                 integer function gnuplot_cmd(ptr_gctrl,cmd) result(status)
!!$                    type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$                    character(len=*), intent(in) :: cmd
!!$                 end function gnuplot_cmd
!!$              end interface

          integer function  gnuplot_cmd(ptr_gctrl,cmd) result(status)
!!$ I will confess I do not understand the purpose of using a va_list in the C implementation of this subroutine (by N. Devillard). Its a string either way. Why worry about types ??

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=GP_CMD_SIZE), intent(in) :: cmd

            integer(i4b), external :: fortran_pipebufsize,fortran_fputs,fortran_fflush

            integer(i4b) :: ierror=0,lencmd=0,max_pipe_buf=0
            character(len=GP_CMD_SIZE) :: pipedcmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_cmd','called with a dissociated pointer',status)
               return
            end if

            lencmd=len_trim(cmd)

            pipedcmd=cmd(1:lencmd)//achar(DECIMALASCIICODEFORNEWLINE)
            lencmd=len_trim(pipedcmd)

            max_pipe_buf=fortran_pipebufsize()

            if(lencmd.ge.max_pipe_buf) then
               status=-7777
               call assert_error('gnuplot_cmd','too long a command passed to pipe buffer',status)
               return
            end if

            if(GNUPLOT_SHOWDEBUG) write(*,'(a)') trim(pipedcmd)

            ierror=fortran_fputs(trim(pipedcmd),ptr_gctrl%fd)

            if(ierror.ne.0) then
               status=-8888
               call assert_error('gnuplot_cmd','failed to write to pipe',status)
               return
            end if

            ierror=fortran_fflush(ptr_gctrl%fd)

            if(ierror.ne.0) then
               status=-6666
               call assert_error('gnuplot_cmd','failed to flush the pipe',status)
               return
            end if

            status=0

          end function gnuplot_cmd

!!$          Function : gnuplot_setstyle
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. plotstyle : character string (plot style - see VALIDPLOTSTYLE)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets the plotting style of a gnuplot session
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = 22 : Invalid plotting style requested, replaced with points
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_setstyle(ptr_gctrl,plotstyle) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: plotstyle
!!$             end function gnuplot_setstyle
!!$          end interface

          integer function gnuplot_setstyle(ptr_gctrl,plotstyle) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: plotstyle

            integer(i4b) :: i=0
            logical :: foundvalid=.false.

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setstyle','called with an dissociated pointer',status)
               return
            end if

            i=1
            foundvalid=.false.
            checkvalid : do while(.not.foundvalid)
               if(plotstyle.eq.trim(VALIDPLOTSTYLE(i))) then
                  ptr_gctrl%plotstyle=plotstyle
                  foundvalid=.true.
                  status=0
                  exit checkvalid
               end if
               i=i+1
               if(i.gt.NUMPLOTSTYLES) exit checkvalid
            end do checkvalid

            if(.not.foundvalid) then
               status=22
               !call assert_warning('gnuplot_setstyle','unknown requested style, '//trim(plotstyle)//', using points')
               ptr_gctrl%plotstyle='points'
            end if

            status=0
          end function gnuplot_setstyle

!!$          Function : gnuplot_setscale
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. axis : character string indicating the axis, such as x,x2,y,y2,z,z2,v,t,u, etc.
!!$                 3. scale : 3-character string (choices - 'LOG' or 'NLG')
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets scaling along different axes
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set requested scaling for the specified axis
!!$
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_setscale(ptr_gctrl,axis,scale) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: axis
!!$               character(len=3), intent(in) :: scale
!!$             end function gnuplot_setscale
!!$          end interface

          integer function gnuplot_setscale(ptr_gctrl,axis,scale) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: axis
            character(len=3), intent(in) :: scale

            character(len=GP_CMD_SIZE) :: cmd
            integer(i4b) :: ierror=0

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setscale','called with an dissociated pointer',status)
               return
            end if

            cmd='set '
            scaling : select case(scale)
            case('LOG')
               cmd=trim(cmd)//' logscale '//trim(axis)
            case('NLG')
               cmd=trim(cmd)//' nologscale '//trim(axis)
            end select scaling

            ierror=gnuplot_cmd(ptr_gctrl,cmd)
            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_setscale','failed to set requested scaling in '//trim(axis),status)
               return
            end if

            status=0

          end function gnuplot_setscale

!!$          Function : gnuplot_setaxislabel
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. axis : character string indicating the axis
!!$                 3. axislabel : character string indicating the label
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets axis label
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set requested label for the axis
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_setaxislabel(ptr_gctrl,axis,axislabel) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: axis,axislabel
!!$             end function gnuplot_setaxislabels
!!$          end interface

          integer function gnuplot_setaxislabel(ptr_gctrl,axis,axislabel) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: axis,axislabel

            integer(i4b) :: ierror=0
            character(len=GP_CMD_SIZE) :: axislabelcmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setaxislabel','called with a dissociated pointer',status)
               return
            end if

            axislabelcmd='set '//trim(axis)//'label "'//trim(axislabel)//'"'
            ierror=gnuplot_cmd(ptr_gctrl,axislabelcmd)

            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_setaxislabel','failed to set label for axis '//trim(axis),status)
               return
            end if
            status=0
          end function gnuplot_setaxislabel


!!$          Function : gnuplot_settitle
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. title : character string (title for the plot)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets the plot title
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set title
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_settitle(ptr_gctrl,title) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: title
!!$             end function gnuplot_settitle
!!$          end interface


          integer function gnuplot_settitle(ptr_gctrl,title) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: title

            integer(i4b) :: ierror=0
            character(len=GP_CMD_SIZE) :: titlecmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_settitle','called with a dissociated pointer',status)
               return
            end if

            titlecmd='set title "'//trim(title)//'"'
            ierror=gnuplot_cmd(ptr_gctrl,titlecmd)

            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_settitle','failed invocation of gnuplot_cmd',status)
               return
            end if
            status=0
          end function gnuplot_settitle

!!$          Function : gnuplot_resetsession
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Resets the session
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = 22 : Failed to delete temporary file
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_resetsession(ptr_gctrl) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$             end function gnuplot_resetsession
!!$          end interface

          integer function gnuplot_resetsession(ptr_gctrl) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            integer(i4b) :: ierror=0,istat=0,i

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_resetsession','called with a dissociated pointer',status)
               return
            end if

            ptr_gctrl%nplots=0
            status=0
          end function gnuplot_resetsession

!!$          Function : gnuplot_set
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. setstring : character string (can be anything that has not been covered by other gnuplot_set type of commands)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets a GNUPlot variable
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set with requested string
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_set(ptr_gctrl,setstring) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: setstring
!!$             end function gnuplot_set
!!$          end interface

          integer function gnuplot_set(ptr_gctrl,setstring) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: setstring

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_set','called with a dissociated pointer',status)
               return
            end if

            cmd='set '//trim(setstring)

            status=gnuplot_cmd(ptr_gctrl,cmd)

            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_set','cannot set requested string',status)
               return
            end if

          end function gnuplot_set

!!$          Function : gnuplot_unset
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. setstring : character string (can be anything that has not been covered by other gnuplot_set type of commands)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Unsets a GNUPlot variable
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot unset with requested string
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_unset(ptr_gctrl,unsetstring) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: unsetstring
!!$             end function gnuplot_unset
!!$          end interface

          integer function gnuplot_unset(ptr_gctrl,unsetstring) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: unsetstring

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_set','called with a dissociated pointer',status)
               return
            end if

            cmd='unset '//trim(unsetstring)

            status=gnuplot_cmd(ptr_gctrl,cmd)

            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_set','cannot unset requested string',status)
               return
            end if

          end function gnuplot_unset

!!$          Function : gnuplot_reset
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Resets the parameters for the current session
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Failed to reset parameters
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_reset(ptr_gctrl) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$             end function gnuplot_reset
!!$          end interface

          integer function gnuplot_reset(ptr_gctrl) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            character(len=GP_CMD_SIZE) :: cmd
            integer(i4b) :: ierror=0

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_reset','called with a dissociated pointer',status)
               return
            end if

            cmd='reset'
            ierror=gnuplot_cmd(ptr_gctrl,cmd)
            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_reset','failed to reset parameters for current session',status)
               return
            end if

            ierror=setacknowledgement(ptr_gctrl,0)

            status=0
          end function gnuplot_reset

!!$          Function : gnuplot_load
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. loadfile : character string : name of the load file
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Runs a gnuplot load file
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot find the specified load file
!!$                 3. status = -2 : You do not have read permissions for the load file
!!$                 4. status = -3 : Cannot run the specified load file
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_load(ptr_gctrl,loadfile) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: loadfile
!!$             end function gnuplot_load
!!$          end interface

          integer function gnuplot_load(ptr_gctrl,loadfile) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: loadfile

            integer(i4b), external :: fortran_access

            character(len=GP_CMD_SIZE) :: cmd
            integer(i4b) :: ierror=0

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_load','called with a dissociated pointer',status)
               return
            end if

            ierror=fortran_access(trim(loadfile),0)
            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_load','cannot find the specified load file '//trim(loadfile),status)
               return
            end if

            ierror=fortran_access(trim(loadfile),4)
            if(status.ne.0) then
               status=-2
               call assert_error('gnuplot_load','you do not have read permissions for the specified load&
                    &file '//trim(loadfile),status)
               return
            end if


            cmd='load "'//trim(loadfile)//'"'

            status=gnuplot_cmd(ptr_gctrl,cmd)

            if(status.ne.0) then
               status=-3
               call assert_error('gnuplot_load','cannot run the specified load file',status)
               return
            end if

          end function gnuplot_load


!!$          Function : gnuplot_replot
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Send a replot command to the session
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Failed to replot
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_replot(ptr_gctrl) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$             end function gnuplot_replot
!!$          end interface

          integer function gnuplot_replot(ptr_gctrl) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_reset','called with a dissociated pointer',status)
               return
            end if

            cmd='replot'
            status=gnuplot_cmd(ptr_gctrl,cmd)
            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_reset','failed to reset parameters for current session',status)
               return
            end if

          end function gnuplot_replot

!!$          Function : gnuplot_setaxisformat
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. axis : character string indicating the axis
!!$                 3. axisformat : character string indicating the format
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Set the format for an axis
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Failed to set format for the specified axis
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_setaxisformat(ptr_gctrl) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: axis,axisformat
!!$             end function gnuplot_setaxisformat
!!$          end interface

          integer function gnuplot_setaxisformat(ptr_gctrl,axis,axisformat) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            character(len=*), intent(in) :: axis,axisformat

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_reset','called with a dissociated pointer',status)
               return
            end if

            cmd='set format '//trim(axis)//' "'//trim(axisformat)//'"'
            status=gnuplot_cmd(ptr_gctrl,cmd)
            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_setaxisformat','failed to set format for axis'//trim(axis),status)
               return
            end if

          end function gnuplot_setaxisformat

!!$          Function : gnuplot_setcoordinatesystem
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. coordsystem : character string indicating the kind of coordinate system (polar or nopolar)
!!$
!!$          Outputs :
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets the coordinate system
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set requested coordinate system
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_setcoordsystem(ptr_gctrl,coordsystem) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: coordsystem
!!$             end function gnuplot_setcoordsystem
!!$          end interface

          integer function gnuplot_setcoordinatesystem(ptr_gctrl,coordsystem) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: coordsystem

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setcoordinatesystem','called with a dissociated pointer',status)
               return
            end if

            cmd='set '//trim(coordsystem)

            status=gnuplot_cmd(ptr_gctrl,cmd)

            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_setcoordinatesystem','cannot set requested coordinate system',status)
               return
            end if

          end function gnuplot_setcoordinatesystem

!!$          Function : gnuplot_setangleunits
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. angleunits : character string indocating the units of angle measurement (radians or degrees) - this may effect some transcendental functions - consult GNUPlot help
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets the angle units
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set requested angle units
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$
!!$             integer function gnuplot_setangleunits(ptr_gctrl,angles) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: angleunits
!!$             end function gnuplot_setangleunits
!!$          end interface

          integer function gnuplot_setangleunits(ptr_gctrl,angleunits) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: angleunits

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setangleunits','called with a dissociated pointer',status)
               return
            end if

            cmd='set '//trim(angleunits)

            status=gnuplot_cmd(ptr_gctrl,cmd)

            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_setangleunits','cannot set requested angle units',status)
               return
            end if

          end function gnuplot_setangleunits

!!$          Function : gnuplot_setrange
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. axis : character string indicating the axis
!!$                 3. lo : real (lower limit of the range)
!!$                 4. hi : real (upper limit of the range)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Sets axes ranges
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot set requested range for the axis
!!$                 3. status = -8888 : Invalid range specification
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_setrange(ptr_gctrl,axis,lo,hi) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: axis
!!$               real(dp), intent(in) :: lo,hi
!!$             end function gnuplot_setrange
!!$          end interface

          integer function gnuplot_setrange(ptr_gctrl,axis,lo,hi) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: axis
            real(dp), intent(in) :: lo,hi

            integer(i4b) :: ierror=0
            character(len=GP_CMD_SIZE) :: rangecmd
            character(len=MAXLENGTH) :: tmp1,tmp2

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_setrange','called with a dissociated pointer',status)
               return
            end if

            if(lo.ge.hi) then
               status=-8888
               call assert_error('gnuplot_setrange','called with invalid '//trim(axis)//'range',status)
               return
            end if

            write(tmp1,'(es12.5)') lo
            write(tmp2,'(es12.5)') hi

            rangecmd='set '//trim(axis)//'range ['//trim(tmp1)//':'//trim(tmp2)//']'
            ierror=gnuplot_cmd(ptr_gctrl,rangecmd)
            if(ierror.ne.0) then
               status=-1
               call assert_error('gnuplot_setrange','failed invocation of gnuplot_cmd',status)
               return
            end if

            status=0

          end function gnuplot_setrange


!!$          Function : gnuplot_pause
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. t : integer time for which to pause
!!$                 3. displaystring : string to be displayed (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Pauses the session
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot pause
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_pause(ptr_gctrl,t,displaystring) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: t
!!$               character(len=*), intent(in), optional :: displaystring
!!$             end function gnuplot_pause
!!$          end interface

          integer function gnuplot_pause(ptr_gctrl,t,displaystring) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: t
            character(len=*), intent(in), optional :: displaystring

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_pause','called with a dissociated pointer',status)
               return
            end if

            cmd='pause '//trim(t)

            if(present(displaystring)) cmd=trim(cmd)//' '//displaystring

            status=gnuplot_cmd(ptr_gctrl,cmd)

            if(status.ne.0) then
               status=-1
               call assert_error('gnuplot_pause','cannot set requested coordinate system',status)
               return
            end if

          end function gnuplot_pause

!!$          Function : gnuplot_plot2d
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. n : integer (number of points to be plotted)
!!$                 3. x : real, dimension(n) (abcissae)
!!$                 4. y1 : real, dimension(n) (first ordinate)
!!$                 5. y1title : character string (title for y1) (optional)
!!$                 6. y1axis : character string (axis for y1) (optional)
!!$                 7. y2 : real, dimension(n) (second ordinate) (optional)
!!$                 8. y2title : character string (title for y2) (optional)
!!$                 9. y2axis : character string (axis for y2) (optional)
!!$                 10. y3 : real, dimension(n) (third ordinate) (optional)
!!$                 11. y3title : character string (title for y3) (optional)
!!$                 12. y3axis : character string (axis for y3) (optional)
!!$                 13. y4 : real, dimension(n) (fourth ordinate) (optional)
!!$                 14. y4title : character string (title for y4) (optional)
!!$                 15. y4axis : character string (axis for y4) (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots x,y1,y2,y3,y4 plot (max of 4 plots possible)
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -3 : Called with less than one point
!!$                 3. status = -1 : Impossible error
!!$                 4. status = -2 : Cannot create temporary file
!!$                 5. status = 22 : Failed to delete temporary file
!!$                 6. status = -4 : Failed to write to temporary file
!!$                 7. status = -5 : Failed to open temporary file in temporary directory
!!$                 8. status = -6 : Failed to execute the plot / replot command
!!$                 9. status = -7 : Invalid y axis specification supplied
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$            integer function gnuplot_plot2d(ptr_gctrl,n,x,y1,y1title,y1axis,y2,y2title,y2axis,&
!!$               y3,y3title,y3axis,y4,y4title,y4axis) result(status)
!!$              type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$              real(dp), dimension(:), intent(in) :: x,y1
!!$              real(dp), dimension(:), intent(in), optional :: y2,y3,y4
!!$              integer(i4b), intent(in) :: n
!!$              character(len=*), intent(in), optional :: y1title,y2title,y3title,y4title
!!$              character(len=*), intent(in), optional :: y1axis,y2axis,y3axis,y4axis
!!$            end function gnuplot_plot2d
!!$          end interface

          integer function gnuplot_plot2d(ptr_gctrl,n,x,y1,y1title,y1axis,y2,y2title,y2axis,&
               y3,y3title,y3axis,y4,y4title,y4axis) result(status)
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            real(dp), dimension(:), intent(in) :: x,y1
            real(dp), dimension(:), intent(in), optional :: y2,y3,y4
            integer(i4b), intent(in) :: n
            character(len=*), intent(in), optional :: y1title,y2title,y3title,y4title
            character(len=*), intent(in), optional :: y1axis,y2axis,y3axis,y4axis

            integer(i4b) :: i=0,temp_fd=0,ierror=0,presenceindex=0
            character(len=MAXLENGTH) :: tmpdirpath=''
            character(len=GP_CMD_SIZE) :: cmd
            logical(lgc) :: datafileexists,getnewdatafile

            status=BAD_SIGNAL
            datafileexists=.false.

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plot2d','called with a dissociated pointer',status)
               return
            end if

            if((size(x).lt.1).or.(size(y1).lt.1).or.(n.lt.1)) then
               status=-3
               call assert_error('gnuplot_plot2d','called with less than 1 point',status)
               return
            end if

            if(present(y2)) then
               if(size(y2).lt.1) then
                  status=-3
                  call assert_error('gnuplot_plot2d','second variable has less than 1 point',status)
               end if
            end if

            if(present(y3)) then
               if(size(y3).lt.1) then
                  status=-3
                  call assert_error('gnuplot_plot2d','third variable has less than 1 point',status)
               end if
            end if

            if(present(y4)) then
               if(size(y4).lt.1) then
                  status=-3
                  call assert_error('gnuplot_plot2d','fourth variable has less than 1 point',status)
               end if
            end if


            tmpdirpath=WORKING_DIRECTORY//ptr_gctrl%tmp_dir
            getnewdatafile=registernewdatafile(ptr_gctrl,trim(tmpdirpath)//'/xy.dat')

            if(.not.getnewdatafile) then
               status=-2
               call assert_error('gnuplot_plot2d','cannot create a temporary file',status)
               return
            end if

!!$ A most interesting effect of interaction of parent and child processes here. File I/O (parent) is slow, so the child could concievably get to the plot command
!!$ before the writing of the file is even complete. This is guarded against by the loop checkifdatafilewrittenyet

            presenceindex=0

            if(present(y2)) then
               presenceindex=presenceindex+1
            end if

            if(present(y3)) then
               presenceindex=presenceindex+1
            end if

            if(present(y4)) then
               presenceindex=presenceindex+1
            end if

            if(present(y1axis)) then
               if(len_trim(y1axis).ne.2) then
                  status=-7
                  call assert_error('gnuplot_plot2d','invalid axis specification for y1',status)
                  return
               end if
            end if

            if(present(y2axis)) then
               if(len_trim(y2axis).ne.2) then
                  status=-7
                  call assert_error('gnuplot_plot2d','invalid axis specification for y2',status)
                  return
               end if
            end if

            if(present(y3axis)) then
               if(len_trim(y3axis).ne.2) then
                  status=-7
                  call assert_error('gnuplot_plot2d','invalid axis specification for y3',status)
                  return
               end if
            end if

            if(present(y4axis)) then
               if(len_trim(y4axis).ne.2) then
                  status=-7
                  call assert_error('gnuplot_plot2d','invalid axis specification for y4',status)
                  return
               end if
            end if

10          format(5(es12.5,1x))

            open(unit=ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),file=trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp)),iostat=ierror)

            if(ierror.ne.0) then
               status=-5
               call assert_error('gnuplot_plot2d','failed to open temporary file in temporary directory',status)
               return
            end if

            writedatatotempfile : do i = 1,n
               presencecases : select case (presenceindex)
               case(0) ! Only x and y1 are present
                  write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),10,iostat=ierror) x(i),y1(i)
               case(1) ! x, y1, and y2 are present
                  write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),10,iostat=ierror) x(i),y1(i),y2(i)
               case(2) ! x, y1, y2, and y3 are present
                  write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),10,iostat=ierror) x(i),y1(i),y2(i),y3(i)
               case(3) ! x, y1, y2, y3, and y4 are present
                  write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),10,iostat=ierror) x(i),y1(i),y2(i),y3(i),y4(i)
               case default
                  status=-1
                  call assert_error('gnuplot_plot2d','impossible value of presenceindex encountered',status)
                  close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                  if(ierror.ne.0) then
                     status=22
                     call assert_warning('gnuplot_plot2d','failed to delete temporary file')
                  end if
                  return
               end select presencecases
               if(ierror.ne.0) then
                  status=-4
                  call assert_error('gnuplot_plot2d','failed to write to temporary file, possibly lack of file space',status)
                  close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                  if(ierror.ne.0) then
                     status=22
                     call assert_warning('gnuplot_plot2d','failed to delete temporary file')
                  end if
                  return
               end if
            end do writedatatotempfile

            close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp))


            checkifdatafilewrittenyet : do while (.not.datafileexists)
               datafileexists=fileexists(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))
               ifdatafileexists : if(datafileexists) then
                  if(ptr_gctrl%nplots.gt.0) then
                     cmd='replot '
                  else
                     cmd='plot '
                  end if

                  cmd=trim(cmd)//' "'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 1:2 ' ! Using 'u' - reduces the size of the buffer to be piped
                  if(present(y1axis)) cmd=trim(cmd)//' axes x1'//trim(y1axis)
                  if(present(y1title)) cmd=trim(cmd)//' t "'//trim(y1title)//'"' ! Using the dimunitive for "title"
                  cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle) ! Using the acronymn for "with"

                  if(present(y2)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 1:3 '
                  if(present(y2axis)) cmd=trim(cmd)//' axes x1'//trim(y2axis)
                  if(present(y2title)) cmd=trim(cmd)//' t "'//trim(y2title)//'"'
                  if(present(y2)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)

                  if(present(y3)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 1:4 '
                  if(present(y3axis)) cmd=trim(cmd)//' axes x1'//trim(y3axis)
                  if(present(y3title)) cmd=trim(cmd)//' t "'//trim(y3title)//'"'
                  if(present(y3)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)

                  if(present(y4)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 1:5 '
                  if(present(y4axis)) cmd=trim(cmd)//' axes x1'//trim(y4axis)
                  if(present(y4title)) cmd=trim(cmd)//' t "'//trim(y4title)//'"'
                  if(present(y4)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)

                  status=gnuplot_cmd(ptr_gctrl,cmd)
                  if(status.ne.0) then
                     status=-6
                     call assert_error('gnuplot_plot2d','could not execute plot / replot command',status)
                  end if

                  ptr_gctrl%nplots=ptr_gctrl%nplots+1
               end if ifdatafileexists
            end do checkifdatafilewrittenyet

          end function gnuplot_plot2d

!!$          Function : gnuplot_fitcurve
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. n : integer (number of points to be plotted)
!!$                 3. np : integer (number of adjustable parameters)
!!$                 4. x : real, dimension(n) (abcissae)
!!$                 5. y : real, dimension(n) (ordinate)
!!$                 6. ffit : character string : function to be fit (must be the RHS of the expression f(x)=...)
!!$                 7. vialist : character(len=1), dimension(np) : the via list (like 'a,b,c' (constants used in ffit)
!!$                 8. paraminit : real, dimension(np) : array containing initial values for the parameters (optional)
!!$                 9. fitlimit : real : FIT_LIMIT , default is usually 1e-6 (depends upon your gnuplot installation) (optional)
!!$                 10. fitlog : character string : filename of the fitting log (optional)
!!$                 11. plotresult : integer whether to plot the result (=1) or not (=0) with the data (optional)
!!$                 12. xlo : real xrange low (optional)
!!$                 13. xhi : real xrange high (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Fits f(x) to supplied y vs x and plots the result if desired
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -3 : Called with less than one point
!!$                 3. status = -1 : Impossible error
!!$                 4. status = -2 : Cannot create temporary file
!!$                 5. status = 22 : Failed to delete temporary file
!!$                 6. status = -4 : Failed to write to temporary file
!!$                 7. status = -5 : No function specified
!!$                 8. status = -7 : Invalid function specification
!!$                  9. status = -8 : Invalid plot directive
!!$                 10. status = 23 : Specified fit log file will overwrite an existing file
!!$                 11. status = -9 : Failed to open temporary file in temporary directory
!!$                 12. status = -10 : Could not execute the fit
!!$                 13. status = -11 : Could not plot the result
!!$                 14. status = -12 : Wrong initial paraminit array size
!!$                 15. status = -13 : Wrong vialist size
!!$                 16. status = -14 : Invalid input xrange
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_fitcurve(ptr_gctrl,n,np,x,y,ffit,vialist,paraminit,xlo,xhi,fitlimit,fitlog,plotresult) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               real(dp), dimension(:), intent(in) :: x,y
!!$               character(len=*), intent(in) :: ffit
!!$               character(len=1), dimension(:) :: vialist
!!$               real(dp), intent(in), optional :: fitlimit
!!$               integer(i4b), intent(in) :: n,np
!!$               integer(i4b), intent(in), optional :: plotresult
!!$               real(dp), dimension(:), intent(in), optional :: paraminit
!!$               character(len=*), intent(in), optional :: fitlog
!!$             end function gnuplot_fitcurve
!!$          end interface

          integer function gnuplot_fitcurve(ptr_gctrl,n,np,x,y,ffit,vialist,paraminit,xlo,xhi,fitlimit,fitlog,plotresult) &
               result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            real(dp), dimension(:), intent(in) :: x,y
            character(len=*), intent(in) :: ffit
            character(len=1), dimension(:),intent(in) :: vialist
            real(dp), intent(in), optional :: fitlimit,xlo,xhi
            integer(i4b), intent(in) :: n,np
            integer(i4b), intent(in), optional :: plotresult
            character(len=*), intent(in), optional :: fitlog
            real(dp), dimension(:), intent(in), optional :: paraminit

            integer(i4b) :: i=0,temp_fd=0,ierror=0
            character(len=MAXLENGTH) :: tmpdirpath=''
            character(len=GP_CMD_SIZE) :: cmd
            character(len=MAXLENGTH) :: junkstring,vialiststring

            logical(lgc) :: datafileexists,getnewdatafile,fexists

            status=BAD_SIGNAL
            datafileexists=.false.

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_fitcurve','called with a dissociated pointer',status)
               return
            end if

            if(((size(x).lt.1).or.(size(y).lt.1).or.(n.lt.1))) then
               status=-3
               call assert_error('gnuplot_fitcurve','called with less than 1 point',status)
               return
            end if
            if(len_trim(ffit).eq.0) then
               status=-5
               call assert_error('gnuplot_fitcurve','no function specified',status)
               return
            end if

            if(index(trim(ffit),'=').ne.0) then
               status=-7
               call assert_error('gnuplot_fitcurve','LHS and = supplied, only RHS needed',status)
               return
            end if

            if(present(plotresult)) then
               if((plotresult.ne.1).and.(plotresult.ne.0)) then
                  status=-8
                  call assert_error('gnuplot_fitcurve','Invalid plot directive',status)
                  return
               end if
            end if

            if(present(fitlog)) then
               fexists=.false.
               inquire(file=trim(fitlog),exist=fexists)
               if(fexists) then
                  status=23
                  call assert_warning('gnuplot_fitcurve','fit log file already exists and will be overwritten')
               end if
            end if

            if(present(paraminit)) then
               if(size(paraminit).ne.np) then
                  status=-12
                  call assert_error('gnuplot_fitcurve','input initial value array is of incorrect size',status)
                  return
               end if
            end if

            if(size(vialist).ne.np) then
               status=-13
               call assert_error('gnuplot_fitcurve','input list of adjustable parameters is of the wrong size',status)
               return
            end if

            if((present(xlo).and..not.present(xhi)).or.(present(xhi).and..not.present(xlo)).or.(xlo.ge.xhi)) then
               status=-14
               call assert_error('gnuplot_fitcurve','invalid x range specification',status)
               return
            end if

11          format(2(es12.5,1x))
            tmpdirpath=WORKING_DIRECTORY//ptr_gctrl%tmp_dir
            getnewdatafile=registernewdatafile(ptr_gctrl,trim(tmpdirpath)//'/xyfit.dat')

            if(.not.getnewdatafile) then
               status=-2
               call assert_error('gnuplot_fitcurve','cannot create a temporary file',status)
               return
            end if

            open(unit=ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),file=trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp)),iostat=ierror)
            if(ierror.ne.0) then
               status=-9
               call assert_error('gnuplot_fitcurve','failed to open temporary file in temporary directory',status)
               return
            end if
            writedatatotempfile : do i = 1,n
               write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),11,iostat=ierror) x(i),y(i)
               if(ierror.ne.0) then
                  status=-4
                  call assert_error('gnuplot_fitcurve','failed to write to temporary file, possibly lack of file space',status)
                  close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                  if(ierror.ne.0) then
                     status=22
                     call assert_warning('gnuplot_fitcurve','failed to delete temporary file')
                  end if
                  return
               end if
            end do writedatatotempfile

            close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp))


            checkifdatafilewrittenyet : do while (.not.datafileexists)
               datafileexists=fileexists(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))
               ifdatafileexists : if(datafileexists) then
                  cmd='f(x) = '//trim(ffit)
                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
                  if(ierror.ne.0) then
                     status=-10
                     call assert_error('gnuplot_fitcurve','fitting interrupted',status)
                     return
                  end if

                  if(present(fitlimit)) then
                     write(junkstring,'(es12.5)') fitlimit
                     cmd='FIT_LIMIT = '//trim(junkstring)
                     ierror=gnuplot_cmd(ptr_gctrl,cmd)
                     if(ierror.ne.0) then
                        status=-10
                        call assert_error('gnuplot_fitcurve','failed to set FIT_LIMIT',status)
                        return
                     end if
                  end if

                  vialiststring=''
                  do i = 1,np
                     vialiststring=trim(vialiststring)//vialist(i)
                     if(i.lt.np) vialiststring=trim(vialiststring)//','
                  end do

                  if(present(paraminit)) then
                     setparaminits : do i = 1,np
                        write(junkstring,'(es12.5)') paraminit(i)
                        cmd=vialist(i)//'='//trim(junkstring)
                        ierror=gnuplot_cmd(ptr_gctrl,cmd)
                        if(ierror.ne.0) then
                           status=-10
                           call assert_error('gnuplot_fitcurve','failed to set initial value for '//&
                                vialist(i)//'='//trim(junkstring),status)
                           return
                        end if
                     end do setparaminits
                  end if

                  cmd='fit '
                  if(present(xlo)) then
                     cmd=trim(cmd)//' ['
                     write(junkstring,'(es12.5)') xlo
                     cmd=trim(cmd)//trim(junkstring)//':'
                     write(junkstring,'(es12.5)') xhi
                     cmd=trim(cmd)//trim(junkstring)//'] '
                  end if

                  cmd=trim(cmd)//' f(x) "'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" via '//trim(vialiststring)
                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
                  if(ierror.ne.0) then
                     status=-10
                     call assert_error('gnuplot_fitcurve','failed to fit curve',status)
                     return
                  end if
                  if(present(fitlog)) then
                     cmd='! mv fit.log '//trim(fitlog)
                     ierror=gnuplot_cmd(ptr_gctrl,cmd)
                     if(ierror.ne.0) then
                        status=-10
                        call assert_error('gnuplot_fitcurve','could not relocate the file fit.log to '//trim(fitlog),status)
                     end if
                  end if

                  checkifplotrequested : if(present(plotresult)) then
                     if(plotresult.eq.1) then
                        cmd='plot '
                        if(present(xlo)) then
                           cmd=trim(cmd)//' ['
                           write(junkstring,'(es12.5)') xlo
                           cmd=trim(cmd)//trim(junkstring)//':'
                           write(junkstring,'(es12.5)') xhi
                           cmd=trim(cmd)//trim(junkstring)//'] '
                        end if
                        cmd=trim(cmd)//' f(x) t "f(x)='//trim(ffit)//'", "'//&
                             trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 1:2 '
                        cmd=trim(cmd)//' t "User data" w '//trim(ptr_gctrl%plotstyle)
                        ierror=gnuplot_cmd(ptr_gctrl,cmd)
                        if(ierror.ne.0) then
                           status=-11
                           call assert_error('gnuplot_fitcurve','could not execute plot the fitted curve',status)
                        end if
                        ptr_gctrl%nplots=ptr_gctrl%nplots+1
                     end if
                  end if checkifplotrequested
               end if ifdatafileexists
            end do checkifdatafilewrittenyet
          end function gnuplot_fitcurve

!!$          Function : gnuplot_plot3d
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. n : integer (number of points along x axis to be plotted)
!!$                 3. m : integer (number of points along y axis to be plotted)
!!$                 4. x : real, dimension(n) (abcissae)
!!$                 5. y : real, dimension(m) (abcissae)
!!$                 6. z1 : real, dimension(n X m) (first surface)
!!$                 7. z1title : character string (title for z1) (optional)
!!$                 8. z2 : real, dimension(n X m) (second surface) (optional)
!!$                 9. z2title : character string (title for z2) (optional)
!!$                 10. z3 : real, dimension(n X m) (third surface) (optional)
!!$                 11. z3title : character string (title for z3) (optional)
!!$                 12. z4 : real, dimension(n X m) (fourth surface) (optional)
!!$                 13. z4title : character string (title for z4) (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots x,y,z1,z2,z3,z4 plot (max of 4 surface plots possible)
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -3 : Called with less than one point
!!$                 3. status = -1 : Impossible error
!!$                 4. status = -2 : Cannot create temporary file
!!$                 5. status = 22 : Failed to delete temporary file
!!$                 6. status = -4 : Failed to write to temporary file
!!$                 7. status = 24 : Failed to unset previously set acknowledgement string
!!$                 8. status = 23 : Failed to set acknowledgement string
!!$                 9. status = -5 : Failed to open temporary file in temporary directory
!!$                10. status = -6 : Failed to execute the splot / replot command
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plot3d(ptr_gctrl,n,m,x,y,z1,z1title,z2,z2title,z3,z3title,z4,z4title) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               real(dp), dimension(:), intent(in) :: x,y
!!$               real(dp), dimension(:,:), intent(in) :: z1
!!$               real(dp), dimension(:,:), intent(in), optional :: z2,z3,z4
!!$               integer(i4b), intent(in) :: n.m
!!$               character(len=*), intent(in), optional :: z1title,z2title,z3title,z4title
!!$             end function gnuplot_plot3d
!!$          end interface

          integer function gnuplot_plot3d(ptr_gctrl,n,m,x,y,z1,z1title,z2,z2title,z3,z3title,z4,z4title) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            real(dp), dimension(:), intent(in) :: x,y
            real(dp), dimension(:,:), intent(in) :: z1
            real(dp), dimension(:,:), intent(in), optional :: z2,z3,z4
            integer(i4b), intent(in) :: n,m
            character(len=*), intent(in), optional :: z1title,z2title,z3title,z4title

            integer(i4b) :: i=0,j=0,temp_fd=0,ierror=0,presenceindex=0
            character(len=MAXLENGTH) :: tmpdirpath=''
            character(len=GP_CMD_SIZE) :: cmd
            logical(lgc) :: datafileexists,getnewdatafile

            status=BAD_SIGNAL
            datafileexists=.false.

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plot3d','called with a dissociated pointer',status)
               return
            end if

            if((size(x).lt.1).or.(size(y).lt.0).or.(size(z1,1).lt.1).or.(size(z1,2).lt.1).or.(n.lt.1).or.(m.lt.1)) then
               status=-3
               call assert_error('gnuplot_plot3d','called with less than 1 point',status)
               return
            end if

            if(present(z2)) then
               if((size(z2,1).lt.1).or.(size(z2,2).lt.1)) then
                  status=-3
                  call assert_error('gnuplot_plot3d','second variable has less than 1 point',status)
                  return
               end if
            end if

            if(present(z3)) then
               if((size(z3,1).lt.1).or.(size(z3,2).lt.1)) then
                  status=-3
                  call assert_error('gnuplot_plot3d','third variable has less than 1 point',status)
                  return
               end if
            end if

            if(present(z4)) then
               if((size(z4,1).lt.1).or.(size(z4,2).lt.1)) then
                  status=-3
                  call assert_error('gnuplot_plot3d','fourth variable has less than 1 point',status)
                  return
               end if
            end if

            if(.not.ptr_gctrl%pubflag) then
               cmd='set nolabel '//trim(ACKTAG)
               ierror=gnuplot_cmd(ptr_gctrl,cmd)
               if(ierror.ne.0) then
                  status=24
                  call assert_warning('gnuplot_plot3d','failed to unset prior acknowledgement string')
               end if
               ierror=setacknowledgement(ptr_gctrl,2)
               if(ierror.ne.0) then
                  status=23
                  call assert_warning('gnuplot_plot3d','failed to set acknowledgement string')
               end if
            end if

            tmpdirpath=WORKING_DIRECTORY//ptr_gctrl%tmp_dir
            getnewdatafile=registernewdatafile(ptr_gctrl,trim(tmpdirpath)//'/xyz.dat')

            if(.not.getnewdatafile) then
               status=-2
               call assert_error('gnuplot_plot3d','cannot create a temporary file',status)
               return
            end if

            presenceindex=0

            if(present(z2)) then
               presenceindex=presenceindex+1
            end if

            if(present(z3)) then
               presenceindex=presenceindex+1
            end if

            if(present(z4)) then
               presenceindex=presenceindex+1
            end if

30          format(6(es12.5,1x))

            open(unit=ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),file=trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp)),&
                 iostat=ierror)
            if(ierror.ne.0) then
               status=-5
               call assert_error('gnuplot_plot3d','failed to open temporary file in temporary directory',status)
               return
            end if
            writedatatotempfilex : do i = 1,n
               writedatatotempfiley : do j = 1,m
                  presencecases : select case (presenceindex)
                  case(0) ! Only x, y and z1 are present
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),30,iostat=ierror) x(i),y(j),z1(i,j)
                  case(1) ! x, y, z1, and z2 are present
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),30,iostat=ierror) x(i),y(j),z1(i,j),z2(i,j)
                  case(2) ! x, y, z1, z2, and z3 are present
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),30,iostat=ierror) x(i),y(j),z1(i,j),z2(i,j),z3(i,j)
                  case(3) ! x, y, z1, z2, z3, and z4 are present
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),30,iostat=ierror) x(i),y(j),z1(i,j),z2(i,j),z3(i,j),&
                          z4(i,j)
                  case default
                     status=-1
                     call assert_error('gnuplot_plot3d','impossible value of presenceindex encountered',status)
                     close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                     if(ierror.ne.0) then
                        status=22
                        call assert_warning('gnuplot_plot3d','failed to delete temporary file')
                     end if
                     return
                  end select presencecases
                  if(ierror.ne.0) then
                     status=-4
                     call assert_error('gnuplot_plot3d','failed to write to temporary file, possibly lack of file space',&
                          status)
                     close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                     if(ierror.ne.0) then
                        status=22
                        call assert_warning('gnuplot_plot3d','failed to delete temporary file')
                     end if
                     return
                  end if
               end do writedatatotempfiley
               write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),'(a)') ''
            end do writedatatotempfilex

            close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp))

            checkifdatafilewrittenyet : do while (.not.datafileexists)
               datafileexists=fileexists(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))
               ifdatafileexists : if(datafileexists) then
                  if(ptr_gctrl%nplots.gt.0) then
                     cmd='replot '
                  else
                     cmd='splot '
                  end if
                  cmd=trim(cmd)//' "'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" using 1:2:3 '
                  if(present(z1title)) cmd=trim(cmd)//' t "'//trim(z1title)//'"'
                  cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)

                  if(present(z2)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" using 1:2:4 '
                  if(present(z2title)) cmd=trim(cmd)//' t "'//trim(z2title)//'"'
                  if(present(z2)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                  if(present(z3)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" using 1:2:5 '
                  if(present(z3title)) cmd=trim(cmd)//' t "'//trim(z3title)//'"'
                  if(present(z3)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                  if(present(z4)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" using 1:2:6 '
                  if(present(z4title)) cmd=trim(cmd)//' t "'//trim(z4title)//'"'
                  if(present(z4)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                  status=gnuplot_cmd(ptr_gctrl,cmd)
                  if(status.ne.0) then
                     status=-6
                     call assert_error('gnuplot_plot3d','could not execute splot / replot command',status)
                  end if
                  ptr_gctrl%nplots=ptr_gctrl%nplots+1
               end if ifdatafileexists
            end do checkifdatafilewrittenyet

          end function gnuplot_plot3d

!!$          Function : gnuplot_fitsurface
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. n : integer (number of points to be plotted along x)
!!$                 3. m : integer (number of points to be plotted along y)
!!$                 4. np : integer (number of adjustable parameters)
!!$                 5. x : real, dimension(n) (abcissae)
!!$                 6. y : real, dimension(m) (abcissae)
!!$                 7. z : real, dimension(n X m) (ordinate)
!!$                 8. gfit : character string : function to be fit (must be the RHS of the expression g(x,y)=...)
!!$                 9. vialist : character(len=1), dimension(np) : the via list (like 'a,b,c' (constants used in ffit)
!!$                 10. paraminit : real, dimension(np) : array containing initial values for the parameters (optional)
!!$                 11. fitlimit : real : FIT_LIMIT , default is usually 1e-6 (depends upon your gnuplot installation) (optional)
!!$                 12. fitlog : character string : filename of the fitting log (optional)
!!$                 13. plotresult : integer whether to plot the result (=1) or not (=0) with the data (optional)
!!$                 14. xlo : real : xrange low (optional)
!!$                 15. xhi : real : xrange high (optional)
!!$                 16. ylo : real : yrange low (optional)
!!$                 17. yhi : real : yrange high (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Fits g(x,y) to supplied z vs x,y and plots the result if desired
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -3 : Called with less than one point
!!$                 3. status = -1 : Impossible error
!!$                 4. status = -2 : Cannot create temporary file
!!$                 6. status = 22 : Failed to delete temporary file
!!$                 7. status = -4 : Failed to write to temporary file
!!$                 8. status = -5 : No function specified
!!$                 9. status = -7 : Invalid function specification
!!$                 10. status = -8 : Invalid plot directive
!!$                 11. status = 23 : Specified fit log file will overwrite an existing file
!!$                 12. status = -9 : Failed to open temporary file temporary directory
!!$                 13. status = -10 : Could not execute the fit
!!$                 14. status = -11 : Could not plot the result
!!$                 15. status = -12 : Wrong initial paraminit array size
!!$                 16. status = -13 : Wrong vialist size
!!$                 17. status = -14 : Invalid input x/y range
!!$                 18. status = 24 : Failed to unset previously set acknowledgement string
!!$                 19. status = 25 : Failed to set acknowledgement string
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_fitsurface(ptr_gctrl,n,m,np,x,y,z,gfit,vialist,paraminit,xlo,xhi,ylo,yhi,fitlimit,fitlog,plotresult) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               real(dp), dimension(:), intent(in) :: x,y
!!$               real(dp), dimension(:,:), intent(in) :: z
!!$               character(len=*), intent(in) :: gfit
!!$               character(len=1), dimension(:) :: vialist
!!$               real(dp), intent(in), optional :: fitlimit
!!$               integer(i4b), intent(in) :: n,m,np
!!$               integer(i4b), intent(in), optional :: plotresult
!!$               real(dp), dimension(:), intent(in), optional :: paraminit
!!$               character(len=*), intent(in), optional :: fitlog
!!$             end function gnuplot_fitsurface
!!$          end interface

          integer function gnuplot_fitsurface(ptr_gctrl,n,m,np,x,y,z,gfit,vialist,paraminit,xlo,xhi,ylo,yhi,&
               fitlimit,fitlog,plotresult) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            real(dp), dimension(:), intent(in) :: x,y
            real(dp), dimension(:,:), intent(in) :: z
            character(len=*), intent(in) :: gfit
            character(len=1), dimension(:),intent(in) :: vialist
            real(dp), intent(in), optional :: fitlimit,xlo,xhi,ylo,yhi
            integer(i4b), intent(in) :: n,m,np
            integer(i4b), intent(in), optional :: plotresult
            character(len=*), intent(in), optional :: fitlog
            real(dp), dimension(:), intent(in), optional :: paraminit

            integer(i4b) :: i=0,j=0,temp_fd=0,ierror=0
            character(len=MAXLENGTH) :: tmpdirpath=''
            character(len=GP_CMD_SIZE) :: cmd
            character(len=MAXLENGTH) :: junkstring,vialiststring

            logical(lgc) :: datafileexists,getnewdatafile,fexists

            status=BAD_SIGNAL
            datafileexists=.false.

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_fitsurface','called with a dissociated pointer',status)
               return
            end if

            if(((size(x).lt.1).or.(size(y).lt.1).or.(n.lt.1).or.(m.lt.1).or.(size(z,1).lt.1).or.(size(z,2).lt.1))) then
               status=-3
               call assert_error('gnuplot_fitsurface','called with less than 1 point',status)
               return
            end if
            if(len_trim(gfit).eq.0) then
               status=-5
               call assert_error('gnuplot_fitsurface','no function specified',status)
               return
            end if

            if(index(trim(gfit),'=').ne.0) then
               status=-7
               call assert_error('gnuplot_fitsurface','LHS and = supplied, only RHS needed',status)
               return
            end if

            if(present(plotresult)) then
               if((plotresult.ne.1).and.(plotresult.ne.0)) then
                  status=-8
                  call assert_error('gnuplot_fitsurface','Invalid plot directive',status)
                  return
               end if
            end if

            if(present(fitlog)) then
               fexists=.false.
               inquire(file=trim(fitlog),exist=fexists)
               if(fexists) then
                  status=23
                  call assert_warning('gnuplot_fitsurface','fit log file already exists and will be overwritten')
               end if
            end if

            if(present(paraminit)) then
               if(size(paraminit).ne.np) then
                  status=-12
                  call assert_error('gnuplot_fitsurface','input initial value array is of incorrect size',status)
                  return
               end if
            end if

            if(size(vialist).ne.np) then
               status=-13
               call assert_error('gnuplot_fitsurface','input list of adjustable parameters is of the wrong size',status)
               return
            end if

            if((present(xlo).and..not.present(xhi)).or.(present(xhi).and..not.present(xlo)).or.(xlo.ge.xhi)) then
               status=-14
               call assert_error('gnuplot_fitsurface','invalid x range specification',status)
               return
            end if

            if((present(ylo).and..not.present(yhi)).or.(present(yhi).and..not.present(ylo)).or.(ylo.ge.xhi)) then
               status=-14
               call assert_error('gnuplot_fitsurface','invalid y range specification',status)
               return
            end if

31          format(3(es12.5,1x))
            tmpdirpath=WORKING_DIRECTORY//ptr_gctrl%tmp_dir
            getnewdatafile=registernewdatafile(ptr_gctrl,trim(tmpdirpath)//'/xyzfit.dat')

            if(.not.getnewdatafile) then
               status=-2
               call assert_error('gnuplot_fitsurface','cannot create a temporary file',status)
               return
            end if

            open(unit=ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),file=trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp)),&
                 iostat=ierror)
            if(ierror.ne.0) then
               status=-9
               call assert_error('gnuplot_fitsurface','failed to open temporary file in temporary directory',status)
               return
            end if
            writedatatotempfilex : do i = 1,n
               writedatatotempfiley : do j = 1,m
                  write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),31,iostat=ierror) x(i),y(j),z(i,j)
                  if(ierror.ne.0) then
                     status=-4
                     call assert_error('gnuplot_fitsurface','failed to write to temporary file, possibly lack of file space',&
                          status)
                     close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                     if(ierror.ne.0) then
                        status=22
                        call assert_warning('gnuplot_fitsurface','failed to delete temporary file')
                     end if
                     return
                  end if
               end do writedatatotempfiley
            end do writedatatotempfilex

            close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp))


            checkifdatafilewrittenyet : do while (.not.datafileexists)
               datafileexists=fileexists(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))
               ifdatafileexists : if(datafileexists) then
                  cmd='g(x,y) = '//trim(gfit)
                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
                  if(ierror.ne.0) then
                     status=-10
                     call assert_error('gnuplot_fitsurface','fitting interrupted',status)
                     return
                  end if

                  if(present(fitlimit)) then
                     write(junkstring,'(es12.5)') fitlimit
                     cmd='FIT_LIMIT = '//trim(junkstring)
                     ierror=gnuplot_cmd(ptr_gctrl,cmd)
                     if(ierror.ne.0) then
                        status=-10
                        call assert_error('gnuplot_fitsurface','failed to set FIT_LIMIT',status)
                        return
                     end if
                  end if

                  vialiststring=''
                  do i = 1,np
                     vialiststring=trim(vialiststring)//vialist(i)
                     if(i.lt.np) vialiststring=trim(vialiststring)//','
                  end do

                  if(present(paraminit)) then
                     setparaminits : do i = 1,np
                        write(junkstring,'(es12.5)') paraminit(i)
                        cmd=vialist(i)//'='//trim(junkstring)
                        ierror=gnuplot_cmd(ptr_gctrl,cmd)
                        if(ierror.ne.0) then
                           status=-10
                           call assert_error('gnuplot_fitsurface','failed to set initial value for '//vialist(i)&
                                //'='//trim(junkstring),status)
                           return
                        end if
                     end do setparaminits
                  end if

                  cmd='fit '
                  if(present(xlo)) then
                     cmd=trim(cmd)//' ['
                     write(junkstring,'(es12.5)') xlo
                     cmd=trim(cmd)//trim(junkstring)//':'
                     write(junkstring,'(es12.5)') xhi
                     cmd=trim(cmd)//trim(junkstring)//'] '
                  end if

                  if(present(ylo)) then
                     cmd=trim(cmd)//' ['
                     write(junkstring,'(es12.5)') ylo
                     cmd=trim(cmd)//trim(junkstring)//':'
                     write(junkstring,'(es12.5)') yhi
                     cmd=trim(cmd)//trim(junkstring)//'] '
                  end if

                  cmd=trim(cmd)//' g(x,y) "'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" via '//trim(vialiststring)
                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
                  if(ierror.ne.0) then
                     status=-10
                     call assert_error('gnuplot_fitsurface','failed to fit surface',status)
                     return
                  end if
                  if(present(fitlog)) then
                     cmd='! mv fit.log '//trim(fitlog)
                     ierror=gnuplot_cmd(ptr_gctrl,cmd)
                     if(ierror.ne.0) then
                        status=-10
                        call assert_error('gnuplot_fitsurface','could not relocate the file fit.log to '//trim(fitlog),status)
                     end if
                  end if

                  checkifplotrequested : if(present(plotresult).and.(plotresult.eq.1)) then
                     if(.not.ptr_gctrl%pubflag) then
                        cmd='set nolabel '//trim(ACKTAG)
                        ierror=gnuplot_cmd(ptr_gctrl,cmd)
                        if(ierror.ne.0) then
                           status=24
                           call assert_warning('gnuplot_fitsurface','failed to unset prior acknowledgement string')
                        end if
                        ierror=setacknowledgement(ptr_gctrl,2)
                        if(ierror.ne.0) then
                           status=25
                           call assert_warning('gnuplot_fitsurface','failed to set acknowledgement string')
                        end if
                     end if
                     cmd='splot '
                     if(present(xlo)) then
                        cmd=trim(cmd)//' ['
                        write(junkstring,'(es12.5)') xlo
                        cmd=trim(cmd)//trim(junkstring)//':'
                        write(junkstring,'(es12.5)') xhi
                        cmd=trim(cmd)//trim(junkstring)//'] '
                     end if

                     if(present(ylo)) then
                        cmd=trim(cmd)//' ['
                        write(junkstring,'(es12.5)') ylo
                        cmd=trim(cmd)//trim(junkstring)//':'
                        write(junkstring,'(es12.5)') yhi
                        cmd=trim(cmd)//trim(junkstring)//'] '
                     end if

                     cmd=trim(cmd)//' g(x,y) t "g(x,y)='//trim(gfit)//'", "'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))&
                          //'" u 1:2:3 '
                     cmd=trim(cmd)//' t "User data" w '//trim(ptr_gctrl%plotstyle)
                     ierror=gnuplot_cmd(ptr_gctrl,cmd)
                     if(ierror.ne.0) then
                        status=-11
                        call assert_error('gnuplot_fitsurface','could not execute plot the fitted surface',status)
                     end if

                     ptr_gctrl%nplots=ptr_gctrl%nplots+1
                  end if checkifplotrequested
               end if ifdatafileexists
            end do checkifdatafilewrittenyet
          end function gnuplot_fitsurface

!!$          Function : gnuplot_plotslope
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. a : real (in a*x+b)
!!$                 3. b : real (in a*x+b)
!!$                 4. xlo : real (lower limit of xrange) (optional) ignored if a previous plot exists
!!$                 5. xhi : real (upper limit of xrange) (optional) ignored if a previous plot exists
!!$                 6. slopetitle : character string (title for slope)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots the specified slope line
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Called with invalid range (if present)
!!$                 3. status = -2 : Failed to plot slope
!!$                 4. status = -4 : Missing one of the two xrange specifications (if present)
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plotslope(ptr_gctrl,a,b,slopetitle,xlo,xhi) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               real(dp), intent(in) :: a,b
!!$               real(dp), intent(in), optional :: xlo,xhi
!!$               character(len=*), intent(in) :: slopetitle
!!$             end function gnuplot_plotslope
!!$          end interface
          integer function gnuplot_plotslope(ptr_gctrl,a,b,slopetitle,xlo,xhi) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            real(dp), intent(in) :: a,b
            character(len=*), intent(in) :: slopetitle
            real(dp), intent(in), optional :: xlo,xhi

            character(len=GP_CMD_SIZE) :: cmd
            character(len=MAXLENGTH) :: tmp1,tmp2

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plotslope','called with a dissociated pointer',status)
               return
            end if
            if(ptr_gctrl%nplots.gt.0) then
               cmd='replot '
            else
               cmd='plot '
            end if

            ifxrangepresent : if(present(xlo).and.(ptr_gctrl%nplots.eq.0)) then ! Cannot set ranges with replot
               if(.not.present(xhi)) then
                  status=-4
                  call assert_error('gnuplot_plotslope','missing one of the two x arguments',status)
                  return
               end if

               if(xlo.ge.xhi) then
                  status=-1
                  call assert_error('gnuplot_plotslope','called with an invalid x range',status)
                  return
               end if

               write(tmp1,'(es12.5)') xlo
               write(tmp2,'(es12.5)') xhi

               cmd=trim(cmd)//' [x = '//trim(tmp1)//':'//trim(tmp2)//'] '
            end if ifxrangepresent

            write(tmp1,'(es12.5)') a
            write(tmp2,'(es12.5)') b

            cmd=trim(cmd)//' '//trim(tmp1)//'*x + '//trim(tmp2)//' '
            cmd=trim(cmd)//' t "'//trim(slopetitle)//'"'
            cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
            status=gnuplot_cmd(ptr_gctrl,cmd)
            if(status.ne.0) then
               status=-2
               call assert_error('gnuplot_plotslope','failed to plot slope',status)
            end if

            ptr_gctrl%nplots=ptr_gctrl%nplots+1
            status=0
          end function gnuplot_plotslope

!!$          Function : gnuplot_plotequation2d
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. equationstring : character string (string defining the function to be plotted)
!!$                 3. equationtitle : character string (title for equation)
!!$                 4. xlo : real (lower limit of xrange) (optional) ignored if a previous plot exists
!!$                 5. xhi : real (upper limit of xrange) (optional) ignored if a previous plot exists
!!$
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots a specified equation
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Called with invalid range (if present)
!!$                 3. status = -2 : Failed to plot equation
!!$                 4. status = -4 : Missing one of the two xrange specifications (if present)
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plotequation2d(ptr_gctrl,equationstring,equationtitle,xlo,xhi) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*) :: equationstring
!!$               real(dp), intent(in), optional :: xlo,xhi
!!$               character(len=*), intent(in) :: equationtitle
!!$             end function gnuplot_plotequation2d
!!$          end interface


          integer function gnuplot_plotequation2d(ptr_gctrl,equationstring,equationtitle,xlo,xhi) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*) :: equationstring
            character(len=*), intent(in) :: equationtitle
            real(dp), intent(in), optional :: xlo,xhi

            character(len=GP_CMD_SIZE) :: cmd
            character(len=MAXLENGTH) :: tmp1,tmp2

            integer(i4b) :: ierror=0

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plotequation2d','called with a dissociated pointer',status)
               return
            end if
            if(ptr_gctrl%nplots.gt.0) then
               cmd='replot '
            else
               cmd='plot '
            end if
            ifxrangepresent : if(present(xlo).and.(ptr_gctrl%nplots.eq.0)) then ! Cannot set ranges with replot
               if(.not.present(xhi)) then
                  status=-4
                  call assert_error('gnuplot_plotequation2d','missing one of the two x arguments for equation plotting',status)
                  return
               end if

               if(xlo.ge.xhi) then
                  status=-1
                  call assert_error('gnuplot_plotequation2d','called with an invalid x range for equation plotting',status)
                  return
               end if
               write(tmp1,'(es12.5)') xlo
               write(tmp2,'(es12.5)') xhi
               cmd=trim(cmd)//' [x='//trim(tmp1)//':'//trim(tmp2)//'] '
            end if ifxrangepresent

            cmd=trim(cmd)//' '//trim(equationstring)//' '
            cmd=trim(cmd)//' t "'//trim(equationtitle)//'"'
            cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
            ierror=gnuplot_cmd(ptr_gctrl,cmd)

            if(ierror.ne.0) then
               status=-2
               call assert_error('gnuplot_plotequation2d','failed to plot equation',status)
            end if

            ptr_gctrl%nplots=ptr_gctrl%nplots+1
            status=0
          end function gnuplot_plotequation2d

!!$          Function : gnuplot_plotequation3d
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. equationstring : character string (string defining the function to be plotted)
!!$                 3. equationtitle : character string (title for equation)
!!$                 4. xlo : real (lower limit of xrange) (optional) ignored if a previous plot exists
!!$                 5. xhi : real (upper limit of xrange) (optional) ignored if a previous plot exists
!!$                 6. ylo : real (lower limit of yrange) (optional) ignored if a previous plot exists
!!$                 7. yhi : real (upper limit of yrange) (optional) ignored if a previous plot exists
!!$
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots a specified equation
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Called with invalid range (if present)
!!$                 3. status = -2 : Failed to plot equation
!!$                 4. status = -4 : Missing one of the two x/y range specifications (if present)
!!$                 5. status = 24 : Failed to unset previously set acknowledgement string
!!$                 6. status = 23 : Failed to set acknowledgement string
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plotequation3d(ptr_gctrl,equationstring,equationtitle,xlo,xhi) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*) :: equationstring
!!$               real(dp), intent(in), optional :: xlo,xhi,ylo,yhi
!!$               character(len=*), intent(in) :: equationtitle
!!$             end function gnuplot_plotequation32
!!$          end interface

          integer function gnuplot_plotequation3d(ptr_gctrl,equationstring,equationtitle,xlo,xhi,ylo,yhi) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*) :: equationstring
            character(len=*), intent(in) :: equationtitle
            real(dp), intent(in), optional :: xlo,xhi,ylo,yhi

            character(len=GP_CMD_SIZE) :: cmd
            character(len=MAXLENGTH) :: tmp1,tmp2
            integer(i4b) :: ierror=0

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plotequation3d','called with a dissociated pointer',status)
               return
            end if
            if(.not.ptr_gctrl%pubflag) then
               cmd='set nolabel '//trim(ACKTAG)
               ierror=gnuplot_cmd(ptr_gctrl,cmd)
               if(ierror.ne.0) then
                  status=24
                  call assert_warning('gnuplot_plotequation3d','failed to unset prior acknowledgement string')
               end if
               ierror=setacknowledgement(ptr_gctrl,2)
               if(ierror.ne.0) then
                  status=23
                  call assert_warning('gnuplot_plotequation3d','failed to set acknowledgement string')
               end if
            end if

            if(ptr_gctrl%nplots.gt.0) then
               cmd='replot '
            else
               cmd='splot '
            end if
            ifxyrangepresent : if(present(xlo).and.(ptr_gctrl%nplots.eq.0)) then ! Cannot set ranges with replot
               if(.not.present(xhi).or..not.present(ylo).or..not.present(yhi)) then
                  status=-4
                  call assert_error('gnuplot_plotequation3d','missing one of the two x/y arguments for equation plotting',status)
                  return
               end if

               if((xlo.ge.xhi).or.(ylo.ge.yhi)) then
                  status=-1
                  call assert_error('gnuplot_plotequation3d','called with an invalid x/y range for equation plotting',status)
                  return
               end if
               write(tmp1,'(es12.5)') xlo
               write(tmp2,'(es12.5)') xhi
               cmd=trim(cmd)//' [x='//trim(tmp1)//':'//trim(tmp2)//'] '

               write(tmp1,'(es12.5)') ylo
               write(tmp2,'(es12.5)') yhi
               cmd=trim(cmd)//' [y='//trim(tmp1)//':'//trim(tmp2)//'] '
            end if ifxyrangepresent

            cmd=trim(cmd)//' '//trim(equationstring)//' '
            cmd=trim(cmd)//' t "'//trim(equationtitle)//'"'
            cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
            ierror=gnuplot_cmd(ptr_gctrl,cmd)
            if(ierror.ne.0) then
               status=-2
               call assert_error('gnuplot_plotequation3d','failed to plot equation',status)
            end if

            ptr_gctrl%nplots=ptr_gctrl%nplots+1
            status=0
          end function gnuplot_plotequation3d



!!$          Function : gnuplot_plotpolarequation
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. equationstring : character string (string defining the function to be plotted)
!!$                 3. equationtitle : character string (title for equation)
!!$
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots a specified equation r=f(t) in polar coordinates - invoke gnuplot_setcoordinatesystem first to set the polar coordinate system
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 3. status = -2 : Failed to plot equation
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plotpolarequation(ptr_gctrl,equationstring,equationtitle) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*) :: equationstring
!!$               character(len=*), intent(in) :: equationtitle
!!$             end function gnuplot_plotpolarequation
!!$          end interface


          integer function gnuplot_plotpolarequation(ptr_gctrl,equationstring,equationtitle) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*) :: equationstring
            character(len=*), intent(in) :: equationtitle

            character(len=GP_CMD_SIZE) :: cmd

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plotpolarequation','called with a dissociated pointer',status)
               return
            end if
            if(ptr_gctrl%nplots.gt.0) then
               cmd='replot '
            else
               cmd='plot '
            end if

            cmd=trim(cmd)//' '//trim(equationstring)//' '
            cmd=trim(cmd)//' t "'//trim(equationtitle)//'"'
            cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
            status=gnuplot_cmd(ptr_gctrl,cmd)
            if(status.ne.0) then
               status=-2
               call assert_error('gnuplot_plotpolarequation','failed to plot equation',status)
            end if

            ptr_gctrl%nplots=ptr_gctrl%nplots+1
            status=0
          end function gnuplot_plotpolarequation

!!$          Function : gnuplot_plotcomplex
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. n : integer (number of points to be plotted)
!!$                 3. c1 : complex, dimension(n)
!!$                 4. c1title : character string (title for c1) (optional)
!!$                 5. c2 : complex, dimension(n) (optional)
!!$                 6. c2title : character string (title for c2) (optional)
!!$                 7. c3 : complex, dimension(n) (optional)
!!$                 8. c3title : character string (title for c3) (optional)
!!$                 9. c4 : complex, dimension(n) (optional)
!!$                 10. c4title : character string (title for c4) (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Polar complex plots c1,c2,c3,c4 plot (max of 4 plots possible) in the z plane
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -3 : Called with less than one point
!!$                 3. status = -1 : Impossible error
!!$                 4. status = -2 : Cannot create temporary file
!!$                 5. status = -4 : Cannot allocate an internal array
!!$                 6. status = -5 : Cannot deallocate an internal array
!!$                 7. status = 22 : Failed to delete temporary file
!!$                 8. status = -6 : Failed to write to temporary file
!!$                 9. status = -7 : Failed to open temporary file in temporary directory
!!$                10. status = -8 : Failed to run plot / replot command
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plotcomplex(ptr_gctrl,n,c1,c1title,c2,c2title,c3,c3title,c4,c4title) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               complex(dp), dimension(:), intent(in) :: c1
!!$               complex(dp), dimension(:), intent(in), optional :: c2,c3,c4
!!$               integer(i4b), intent(in) :: n
!!$               character(len=*), intent(in), optional :: c1title,c2title,c3title,c4title
!!$             end function gnuplot_plotcomplex
!!$          end interface

          integer function gnuplot_plotcomplex(ptr_gctrl,n,c1,c1title,c2,c2title,c3,c3title,c4,c4title) result(status)
            use polar_cartesian

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            complex(dp), dimension(:), intent(in) :: c1
            complex(dp), dimension(:), intent(in), optional :: c2,c3,c4
            integer(i4b), intent(in) :: n
            character(len=*), intent(in), optional :: c1title,c2title,c3title,c4title

            integer(i4b) :: i=0,temp_fd=0,ierror=0,presenceindex=0,istat=0
            character(len=MAXLENGTH) :: tmpdirpath=''
            character(len=GP_CMD_SIZE) :: cmd
            logical(lgc) :: datafileexists,getnewdatafile
            type(polar), dimension(:), allocatable :: p1,p2,p3,p4

            status=BAD_SIGNAL
            datafileexists=.false.

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_plotcomplex','called with a dissociated pointer',status)
               return
            end if

            if((size(c1).lt.1).or.(n.lt.1)) then
               status=-3
               call assert_error('gnuplot_plotcomplex','called with less than 1 point',status)
               return
            end if

            if(present(c2)) then
               if(size(c2).lt.1) then
                  status=-3
                  call assert_error('gnuplot_plotcomplex','second variable has less than 1 point',status)
                  return
               end if
            end if

            if(present(c3)) then
               if(size(c3).lt.1) then
                  status=-3
                  call assert_error('gnuplot_plotcomplex','third variable has less than 1 point',status)
                  return
               end if
            end if

            if(present(c4)) then
               if(size(c4).lt.1) then
                  status=-3
                  call assert_error('gnuplot_plotcomplex','fourth variable has less than 1 point',status)
                  return
               end if
            end if

            tmpdirpath=WORKING_DIRECTORY//ptr_gctrl%tmp_dir
            getnewdatafile=registernewdatafile(ptr_gctrl,trim(tmpdirpath)//'/rt.dat')

            if(.not.getnewdatafile) then
               status=-2
               call assert_error('gnuplot_plotcomplex','cannot create a temporary file',status)
               return
            end if

            presenceindex=0

            if(present(c2)) then
               presenceindex=presenceindex+1
            end if

            if(present(c3)) then
               presenceindex=presenceindex+1
            end if

            if(present(c4)) then
               presenceindex=presenceindex+1
            end if

20          format(8(es12.5,1x))

            open(unit=ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),file=trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp)),iostat=ierror)
            if(ierror.ne.0) then
               status=-5
               call assert_error('gnuplot_plotcomplex','failed to open temporary file in temporary directory',status)
               return
            end if
            allocatearrays : select case (presenceindex)
               case(0)
                  allocate(p1(n),stat=istat)
                  if(istat.ne.0) then
                     status=-4
                     call assert_error('gnuplot_plotcomplex','cannot allocate array for polar type',status)
                     return
                  end if
               case(1)
                  allocate(p1(n),p2(n),stat=istat)
                  if(istat.ne.0) then
                     status=-4
                     call assert_error('gnuplot_plotcomplex','cannot allocate array for polar type',status)
                     return
                  end if
               case(2)
                  allocate(p1(n),p2(n),p3(n),stat=istat)
                  if(istat.ne.0) then
                     status=-4
                     call assert_error('gnuplot_plotcomplex','cannot allocate array for polar type',status)
                     return
                  end if
               case(3)
                  allocate(p1(n),p2(n),p3(n),p4(n),stat=istat)
                  if(istat.ne.0) then
                     status=-4
                     call assert_error('gnuplot_plotcomplex','cannot allocate array for polar type',status)
                     return
                  end if
               end select allocatearrays

               writedatatotempfile : do i = 1,n
                  presencecases : select case (presenceindex)
                  case(0) ! Only c1 is present
                     p1(i)=c1(i)
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),20,iostat=ierror) p1(i)
                  case(1) ! c1, and c2 are present
                     p1(i)=c1(i)
                     p2(i)=c2(i)
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),20,iostat=ierror) p1(i),p2(i)
                  case(2) ! c1, c2, and c3 are present
                     p1(i)=c1(i)
                     p2(i)=c2(i)
                     p3(i)=c3(i)
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),20,iostat=ierror) p1(i),p2(i),p3(i)
                  case(3) ! c1, c2, c3, and c4 are present
                     p1(i)=c1(i)
                     p2(i)=c2(i)
                     p3(i)=c3(i)
                     p4(i)=c4(i)
                     write(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),20,iostat=ierror) p1(i),p2(i),p3(i),p4(i)
                  case default
                     status=-1
                     call assert_error('gnuplot_plotcomplex','impossible value of presenceindex encountered',&
                          status)
                     close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                     if(ierror.ne.0) then
                        status=22
                        call assert_warning('gnuplot_plotcomplex','failed to delete temporary file')
                     end if
                     return
                  end select presencecases
                  if(ierror.ne.0) then
                     status=-6
                     call assert_error('gnuplot_plotcomplex','failed to write to temporary file,&
                          &possibly lack of file space',status)
                     close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp),status='delete',iostat=ierror)
                     if(ierror.ne.0) then
                        status=22
                        call assert_warning('gnuplot_plotcomplex','failed to delete temporary file')
                     end if
                     return
                  end if
               end do writedatatotempfile
               close(ptr_gctrl%DATAFILELUNIT(ptr_gctrl%ntmp))

               deallocatearrays : select case (presenceindex)
               case(0)
                  deallocate(p1,stat=istat)
                  if(istat.ne.0) then
                     status=-5
                     call assert_error('gnuplot_plotcomplex','cannot deallocate array for polar type',status)
                     return
                  end if
               case(1)
                  deallocate(p1,p2,stat=istat)
                  if(istat.ne.0) then
                     status=-5
                     call assert_error('gnuplot_plotcomplex','cannot deallocate array for polar type',status)
                     return
                  end if
               case(2)
                  deallocate(p1,p2,p3,stat=istat)
                  if(istat.ne.0) then
                     status=-5
                     call assert_error('gnuplot_plotcomplex','cannot deallocate array for polar type',status)
                     return
                  end if
               case(3)
                  deallocate(p1,p2,p3,p4,stat=istat)
                  if(istat.ne.0) then
                     status=-5
                     call assert_error('gnuplot_plotcomplex','cannot deallocate array for polar type',status)
                     return
                  end if
               end select deallocatearrays
               checkifdatafilewrittenyet : do while (.not.datafileexists)
                  datafileexists=fileexists(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))
                  ifdatafileexists : if(datafileexists) then
                     if(ptr_gctrl%nplots.gt.0) then
                        cmd='replot '
                     else
                        cmd='plot '
                     end if
                     cmd=trim(cmd)//' "'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 2:1 '
                     if(present(c1title)) cmd=trim(cmd)//' t "'//trim(c1title)//'"'
                     cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                     if(present(c2)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 4:3 '
                     if(present(c2title)) cmd=trim(cmd)//' t "'//trim(c2title)//'"'
                     if(present(c2)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                     if(present(c3)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 6:5 '
                     if(present(c3title)) cmd=trim(cmd)//' t "'//trim(c3title)//'"'
                     if(present(c3)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                     if(present(c4)) cmd=trim(cmd)//',"'//trim(ptr_gctrl%DATAFILENAME(ptr_gctrl%ntmp))//'" u 8:7 '
                     if(present(c4title)) cmd=trim(cmd)//' t "'//trim(c4title)//'"'
                     if(present(c4)) cmd=trim(cmd)//' w '//trim(ptr_gctrl%plotstyle)
                     status=gnuplot_cmd(ptr_gctrl,cmd)
                     if(status.ne.0) then
                        status=-8
                        call assert_error('gnuplot_plotcomplex','could not execute plot / replot command',status)
                     end if
                     ptr_gctrl%nplots=ptr_gctrl%nplots+1
                  end if ifdatafileexists
               end do checkifdatafilewrittenyet
             end function gnuplot_plotcomplex

!!$          Function : gnuplot_hardcopy
!!$
!!$          Inputs :
!!$
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$                 2. plotfileformat : character string (format of the requested file)
!!$                 3. plotfilename : character string (file name for the plot file)
!!$                 4. extraarguments : character string (any user specified commands that come after "set terminal <terminal-type>", 
!!$                                     such as the enhanced mode for PostScript) (optional)
!!$                 5. pubflag : 3 character string (optional) (if set to 'PUB', the acknowledgement and date string will be suppressed)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Produces a hard(soft? :) )copy of the plot in a user specified format. Valid file formats are :
!!$                   a. 'PS' : PostScript
!!$                   b. 'EPS' : Encapsulated PostScript
!!$                   c. 'PDF' : Portable Document Format
!!$                   d. 'PNG' : Portable Network Graphics
!!$                   e. 'LATEX' : LaTeX format
!!$                   f. 'PSLATEX'
!!$                   g. 'EPSLATEX'
!!$                   h. 'PSTRICKS'
!!$                   i. 'HPGL' : HP Graphics Language
!!$                   j. 'HPDJ'
!!$                   k. 'FIG'
!!$                   l. 'DXF' : AutoCad format
!!$                   m. 'EEPIC'
!!$
!!$                  in addition to x11, which of course does not generate any "hardcopy" :)
!!$
!!$                    This function must be called immediately after gnuplot_init and before any plotting or set commands. 
!!$                    In case an unknown format is chosen, the function reverts to PostScript.
!!$                    If you choose the PDF format, you must have the epstopdf command available in your $PATH variable. The actual
!!$                    conversion from intermediary ps to pdf takes place in gnuplot_genpdf. 
!!$                    So gnuplot_genpdf must be called after this to ensure that the PDF file is generated
!!$                    The function does no file extension checking, so make sure that you supply the correct output filename.
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = 22 : Unknown format requested
!!$                 3. status = -2 : Failed to set requested (valid) terminal type
!!$                 4. status = -3 : Failed to set output file
!!$                 5. status = 23 : Failed to set acknowledgement string
!!$                 6. status = 24 : Failed to unset prior acknowledgement label
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_hardcopy(ptr_gctrl,plotfileformat,plotfilename,extraarguments,pubflag) result(status)
!!$               type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$               character(len=*), intent(in) :: plotfileformat,plotfilename
!!$               character(len=*), intent(in), optional :: extraarguments
!!$               character(len=3), intent(in), optional :: pubflag
!!$             end function gnuplot_hardcopy
!!$          end interface

          integer function gnuplot_hardcopy(ptr_gctrl,plotfileformat,plotfilename,extraarguments,pubflag) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=*), intent(in) :: plotfileformat,plotfilename
            character(len=*), intent(in), optional :: extraarguments
            character(len=3), intent(in), optional :: pubflag

            character(len=GP_CMD_SIZE) :: cmd
            integer(i4b) :: ierror=0,pubint=1

            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_hardcopy','called with a dissociated pointer',status)
               return
            end if

            ptr_gctrl%hardcopyformat=trim(plotfileformat)
            ptr_gctrl%hardcopyfilename=trim(plotfilename)
            ptr_gctrl%hardcopyfilenamelength=len_trim(plotfilename)

            cmd='set terminal '

            selectformat : select case(trim(plotfileformat))
               case('X11')
                  cmd=trim(cmd)//' x11'
               case('PS')
                  cmd=trim(cmd)//' postscript '
               case('PDF')
                  cmd=trim(cmd)//' postscript '
                  ptr_gctrl%hardcopyfilename=plotfilename(1:ptr_gctrl%hardcopyfilenamelength-3)//'ps'
                  ptr_gctrl%hardcopyfilenamelength=len_trim(ptr_gctrl%hardcopyfilename)
               case('EPS')
                  cmd=trim(cmd)//' postscript eps '
               case('PNG')
                  cmd=trim(cmd)//' png '
               case('LATEX')
                  cmd=trim(cmd)//' latex '
               case('PSLATEX')
                  cmd=trim(cmd)//' pslatex '
               case('EPSLATEX')
                  cmd=trim(cmd)//' epslatex '
               case('PSTRICKS')
                  cmd=trim(cmd)//' pstricks '
               case('HPGL')
                  cmd=trim(cmd)//' hpgl '
               case('HPDJ')
                  cmd=trim(cmd)//' hpdj '
               case('FIG')
                  cmd=trim(cmd)//' fig '
               case('DXF')
                  cmd=trim(cmd)//' dxf'
               case('EEPIC')
                  cmd=trim(cmd)//' eepic '
               case default
                  cmd=trim(cmd)//' postscript '
                  status=22
                  call assert_warning('gnuplot_hardcopy','unknown format requested, reverting to PostScript')
                  ptr_gctrl%hardcopyformat='PS'
               end select selectformat

               if(present(extraarguments)) then
                  cmd=trim(cmd)//' '//trim(extraarguments)
               end if

               ierror=gnuplot_cmd(ptr_gctrl,cmd)
               if(ierror.ne.0) then
                  status=-2
                  call assert_error('gnuplot_hardcopy','failed to set terminal type',status)
                  return
               end if

               ptr_gctrl%pubflag=.true.

               pubint=1

               if(.not.present(pubflag)) then
                  pubint=0
               else
                  cmd='set nolabel '//trim(ACKTAG)
                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
                  if(ierror.ne.0) then
                     status=24
                     call assert_warning('gnuplot_hardcopy','failed to unset prior acknowledgement string')
                  end if
                  if(pubflag.ne.'PUB') pubint=0
               end if


               ackdisplay : if(pubint.eq.0) then
                  ptr_gctrl%pubflag=.false.
!!$                  cmd='set nolabel '//trim(ACKTAG)
!!$                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
!!$                  if(ierror.ne.0) then
!!$                     status=24
!!$                     call assert_warning('gnuplot_hardcopy','failed to unset prior acknowledgement string')
!!$                  end if
                  ierror=setacknowledgement(ptr_gctrl,1)
                  if(ierror.ne.0) then
                     status=23
                     call assert_warning('gnuplot_hardcopy','failed to set acknowledgement string')
                  end if
               end if ackdisplay
               status=0

               if(trim(plotfileformat).ne.'X11') then
                  cmd='set output "'//trim(ptr_gctrl%hardcopyfilename)//'"'
                  ierror=gnuplot_cmd(ptr_gctrl,cmd)
                  if(ierror.ne.0) then
                     status=-3
                     call assert_error('gnuplot_hardcopy','failed to set output file',status)
                     return
                  end if

                  ptr_gctrl%hardcopystatus=status
               else
                  ptr_gctrl%hardcopystatus=1
               end if
             end function gnuplot_hardcopy

!!$          Function : gnuplot_genpdf
!!$
!!$          Inputs :
!!$                 1. ptr_gctrl : pointer to structure of type gnuplot_ctr
!!$
!!$          Outputs :
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Produces a pdf - must be called after all the plot commands you want included have been executed
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Cannot generate PDF from the PS file.
!!$                 3. status = 22 : Called when the o/p format was not PDF
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$            integer function gnuplot_genpdf(ptr_gctrl) result(status)
!!$              type(gnuplot_ctrl), pointer :: ptr_gctrl
!!$            end function gnuplot_genpdf
!!$          end interface

          integer function gnuplot_genpdf(ptr_gctrl) result(status)

            type(gnuplot_ctrl), pointer :: ptr_gctrl

            character(len=GP_CMD_SIZE) :: cmd
            integer(i4b) :: ierror=0,hardcopylength=0
            character(len=MAXLENGTH) :: pdffile


            status=BAD_SIGNAL

            if(.not.associated(ptr_gctrl)) then
               call assert_error('gnuplot_genpdf','called with a dissociated pointer',status)
               return
            end if

            hardcopylength=len_trim(ptr_gctrl%hardcopyformat)

            if(ptr_gctrl%hardcopyformat(1:hardcopylength).ne.'PDF') then
               status=-2
               call assert_warning('gnuplot_genpdf','requested output format was not PDF')
               return
            end if

            pdfconvert : if((hardcopylength.ne.0).and.(ptr_gctrl%hardcopystatus.eq.0)) then
               pdffile=ptr_gctrl%hardcopyfilename(1:ptr_gctrl%hardcopyfilenamelength-2)//'pdf'
               cmd='! epstopdf -d --nogs -hires '//trim(ptr_gctrl%hardcopyfilename)//' > '//trim(WORKING_DIRECTORY//&
                    ptr_gctrl%tmp_dir)//'/testcorr.ps; epstopdf --outfile '//trim(WORKING_DIRECTORY)//trim(pdffile)//' '//&
                    trim(WORKING_DIRECTORY//ptr_gctrl%tmp_dir)//'/testcorr.ps; rm '//trim(WORKING_DIRECTORY//ptr_gctrl%tmp_dir)//&
                    '/*.ps; '//'rm '//trim(ptr_gctrl%hardcopyfilename)
               ierror=gnuplot_cmd(ptr_gctrl,cmd)
               if(ierror.ne.0) then
                  status=-8888
                  call assert_error('gnuplot_genpdf','failed to convert PS file to PDF file',status)
                  return
               end if
            end if pdfconvert

            status=0
          end function gnuplot_genpdf

!!$          Function : gnuplot_plot2dset
!!$
!!$          Inputs :
!!$
!!$                 1. pf : character string (plot file name)
!!$                 2. pformat : character string (format of hardcopy)
!!$                 3. commontitle : character string (title used for the whole plot)
!!$                 4. xarray : real, dimension(n) (x array)
!!$                 5. xlabel : character string (x label)
!!$                 6. ylabel : character string (label for the y1 axis)
!!$                 7. y2label : character string (label for the y2 axis) (optional)
!!$                 8. yarray :  real, dimension(n) (y array)
!!$                 9. yaxis : character string (len=2) (axis used - y1 or y2) (optional)
!!$                 10. ykey : character string (ykey)
!!$                 11. yarray1 : real, dimension(n) (yarray1) (optional)
!!$                 12. yaxis1 : character string (len=2) (axis used - y1 or y2) (optional)
!!$                 13. ykey1 : character string (ykey1) (optional)
!!$                 14. yarray2 : real, dimension(n) (yarray2) (optional)
!!$                 15. yaxis2 : character string (len=2) (axis used - y1 or y2) (optional)
!!$                 16. ykey2 : character string (ykey2) (optional)
!!$                 17. yarray3 : real, dimension(n) (yarray3) (optional)
!!$                 18. yaxis3 : character string (len=2) (axis used - y1 or y2) (optional)
!!$                 19. ykey3 : character string (ykey3) (optional)
!!$                 20. logy1 : character string (len=3) (log scale or not -- for y1 axis) (optional)
!!$                 21. logy2 : character string (len=3) (log scale or not -- for y2 axis) (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots x,y,y1,y2,y3 plot (max of 4 line plots possible)
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Array size mismatch
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plot2dset(pf,pformat,commontitle,xarray,xlabel,ylabel,y2label,yarray,yaxis,ykey,&
!!$                  yarray1,yaxis1,ykey1,yarray2,yaxis2,ykey2,yarray3,yaxis3,ykey3,logy1,logy2) result(status)
!!$               character(len=*), intent(in) :: pf,pformat,commontitle
!!$               real(dp), dimension(:), intent(in) :: xarray,yarray
!!$               character(len=*), intent(in) :: xlabel,ylabel,ykey
!!$               real(dp), dimension(:), intent(in), optional :: yarray1,yarray2,yarray3
!!$               character(len=*), intent(in), optional :: y2label,ykey1,ykey2,ykey3
!!$               character(len=*), intent(in), optional :: yaxis,yaxis1,yaxis2,yaxis3
!!$               character(len=*), intent(in), optional :: logy1,logy2
!!$             end function gnuplot_plot2dset
!!$          end interface

          integer function gnuplot_plot2dset(pf,pformat,commontitle,xarray,xlabel,ylabel,y2label,yarray,yaxis,ykey,&
               yarray1,yaxis1,ykey1,yarray2,yaxis2,ykey2,yarray3,yaxis3,ykey3,logy1,logy2) result(status)
            character(len=*), intent(in) :: pf,pformat,commontitle
            real(dp), dimension(:), intent(in) :: xarray,yarray
            character(len=*), intent(in) :: xlabel,ylabel,ykey
            real(dp), dimension(:), intent(in), optional :: yarray1,yarray2,yarray3
            character(len=*), intent(in), optional :: y2label,ykey1,ykey2,ykey3
            character(len=*), intent(in), optional :: yaxis,yaxis1,yaxis2,yaxis3
            character(len=*), intent(in), optional :: logy1,logy2

            integer(i4b) :: numplots,ix,iy,i
            type(gnuplot_ctrl), pointer :: ptr_gctrl
            character(len=2) :: yaxisused,yaxis1used,yaxis2used,yaxis3used

            status=BAD_SIGNAL

            yaxisused='y1'
            yaxis1used='y1'
            yaxis2used='y1'
            yaxis3used='y1'

            numplots=1

            if(present(yarray1)) numplots=numplots+1
            if(present(yarray2)) numplots=numplots+1
            if(present(yarray3)) numplots=numplots+1

            ix=size(xarray,1)
            iy=size(yarray,1)

            if(ix.ne.iy) then
               status=-1
               call assert_error('gnuplot_plot2dset','array size mismatch',status)
!!$               write(*,'("array size mismatch : ",i3," vs ",i3,".")') &
!!$                    ix,iy
               return
            end if

            GNUPLOT_SHOWWARNINGS=.false.
            GNUPLOT_SHOWDEBUG=.false.

            i=0
            ptr_gctrl=>gnuplot_init('-persist')

            if(.not.associated(ptr_gctrl)) then
               status=BAD_SIGNAL
               call assert_error('gnuplot_plot2dset','called with a dissociated pointer',status)
               return
            end if

            status=gnuplot_hardcopy(ptr_gctrl,trim(pformat),trim(pf),'enhanced color solid','PUB')

            status=gnuplot_setstyle(ptr_gctrl,'lines')
            status=gnuplot_settitle(ptr_gctrl,trim(commontitle))
            status=gnuplot_setaxislabel(ptr_gctrl,'x',trim(xlabel))
            status=gnuplot_setaxislabel(ptr_gctrl,'y',trim(ylabel))
            if(present(logy1)) then
               if(logy1.eq.'y') status=gnuplot_setscale(ptr_gctrl,'y','LOG')
            end if
            if(present(y2label)) then
               if(present(logy2)) then
                  if(logy2.eq.'y') status=gnuplot_setscale(ptr_gctrl,'y2','LOG')
               end if
               status=gnuplot_set(ptr_gctrl,'y2tics')
               status=gnuplot_setaxislabel(ptr_gctrl,'y2',trim(y2label))
               if(present(yaxis)) yaxisused=yaxis
               if(present(yaxis1)) yaxis1used=yaxis1
               if(present(yaxis2)) yaxis2used=yaxis2
               if(present(yaxis3)) yaxis3used=yaxis3
            end if
            status=gnuplot_set(ptr_gctrl,'size 1.0,1.0')
            selectplot2dcmd : select case(numplots)
            case(1)
               status=gnuplot_plot2d(ptr_gctrl,ix,xarray(1:ix),yarray(1:ix),y1axis=yaxisused,y1title=trim(ykey))
            case(2)
               status=gnuplot_plot2d(ptr_gctrl,ix,xarray(1:ix),yarray(1:ix),y1axis=yaxisused,y1title=trim(ykey),y2=yarray1(1:ix),&
                    y2axis=yaxis1used,y2title=trim(ykey1))
            case(3)
               status=gnuplot_plot2d(ptr_gctrl,ix,xarray(1:ix),yarray(1:ix),y1axis=yaxisused,y1title=trim(ykey),y2=yarray1(1:ix),&
                    y2axis=yaxis1used,y2title=trim(ykey1),y3=yarray2(1:ix),y3axis=yaxis2used,y3title=trim(ykey2))
            case(4)
               status=gnuplot_plot2d(ptr_gctrl,ix,xarray(1:ix),yarray(1:ix),y1axis=yaxisused,y1title=trim(ykey),y2=yarray1(1:ix),&
                    y2axis=yaxis1used,y2title=trim(ykey1),y3=yarray2(1:ix),y3axis=yaxis2used,y3title=trim(ykey2),y4=yarray3(1:ix),&
                    y4axis=yaxis3used,y4title=trim(ykey3))
            end select selectplot2dcmd

            status=gnuplot_set(ptr_gctrl,'output')
            if(trim(pformat).eq.'PDF') status=gnuplot_genpdf(ptr_gctrl)
            status=gnuplot_close(ptr_gctrl)
          end function gnuplot_plot2dset

!!$          Function : gnuplot_plot3dpanel
!!$
!!$          Inputs :
!!$
!!$                 1. pf : character string (plot file name)
!!$                 2. pformat : character string (format of hardcopy)
!!$                 3. commontitle : character string (title used for the whole plot)
!!$                 4. xarray : real, dimension(n) (x array)
!!$                 5. xlabel : character string (x label)
!!$                 6. yarray : real, dimension(m) (y array)
!!$                 7. ylabel : character string (y label)
!!$                 8. zarray : real, dimension (n X m) (z array)
!!$                 9. zlabel : character string (z label)
!!$                 10. zarray1 : real, dimension (n X m) (z1 array) (optional)
!!$                 11. zlabel1 : character string (z1 label) (optional)
!!$                 12. zarray2 : real, dimension (n X m) (z2 array) (optional)
!!$                 13. zlabel2 : character string (z2 label) (optional)
!!$                 14. zarray3 : real, dimension (n X m) (z3 array) (optional)
!!$                 15. zlabel3 : character string (z3 label) (optional)
!!$                 16. logz : character string (len=3) (log scale or not) (optional)
!!$
!!$          Outputs :
!!$
!!$                 1. status : integer (error flag)
!!$
!!$          Purpose : Plots x,y,z,z1,z2,z3 plot (max of 4 surface plots possible)
!!$
!!$          Diagnostics :
!!$
!!$                 0. status = 0 : Success
!!$                 1. status = BAD_SIGNAL : Called with a dissociated pointer
!!$                 2. status = -1 : Array size mismatch
!!$
!!$          Type : public
!!$
!!$          Interface :
!!$
!!$          interface
!!$             integer function gnuplot_plot3dpanel(pf,pformat,commontitle,xarray,xlabel,yarray,ylabel,zarray,zlabel,&
!!$                  zarray1,zlabel1,zarray2,zlabel2,zarray3,zlabel3,logz) result(status)
!!$               character(len=*), intent(in) :: pf,pformat,commontitle
!!$               real(dp), dimension(:), intent(in) :: xarray,yarray
!!$               character(len=*), intent(in) :: xlabel,ylabel,zlabel
!!$               real(dp), dimension(:,:), intent(in) :: zarray
!!$               real(dp), dimension(:,:), intent(in), optional :: zarray1,zarray2,zarray3
!!$               character(len=*), intent(in), optional :: zlabel1,zlabel2,zlabel3
!!$               character(len=*), intent(in), optional :: logz
!!$             end function gnuplot_plot3dpanel
!!$          end interface

          integer function gnuplot_plot3dpanel(pf,pformat,commontitle,xarray,xlabel,yarray,ylabel,zarray,zlabel,&
               zarray1,zlabel1,zarray2,zlabel2,zarray3,zlabel3,logz) result(status)
            character(len=*), intent(in) :: pf,pformat,commontitle
            real(dp), dimension(:), intent(in) :: xarray,yarray
            character(len=*), intent(in) :: xlabel,ylabel,zlabel
            real(dp), dimension(:,:), intent(in) :: zarray
            real(dp), dimension(:,:), intent(in), optional :: zarray1,zarray2,zarray3
            character(len=*), intent(in), optional :: zlabel1,zlabel2,zlabel3
            character(len=*), intent(in), optional :: logz


            integer(i4b) :: numplots,ix,iy,i
            type(gnuplot_ctrl), pointer :: ptr_gctrl

            status=BAD_SIGNAL
            numplots=1

            if(present(zarray1)) numplots=numplots+1
            if(present(zarray2)) numplots=numplots+1
            if(present(zarray3)) numplots=numplots+1

            ix=size(xarray,1)
            iy=size(yarray,1)

            if((size(zarray,1).ne.ix).or.(size(zarray,2).ne.iy)) then
               status=-1
               call assert_error('gnuplot_plot3dpanel','z array size mismatch',status)
!!$               write(*,'("z array size mismatch : ",i3," X ",i3,/,"Instead of ",i3," X ",i3,".")') &
!!$                    size(zarray,1),size(zarray,2),ix,iy
               return
            end if

            GNUPLOT_SHOWWARNINGS=.false.
            GNUPLOT_SHOWDEBUG=.false.

            i=0
            ptr_gctrl=>gnuplot_init('-persist')

            if(.not.associated(ptr_gctrl)) then
               status=BAD_SIGNAL
               call assert_error('gnuplot_plot3dpanel','called with a dissociated pointer',status)
               return
            end if

            status=gnuplot_hardcopy(ptr_gctrl,trim(pformat),trim(pf),'enhanced color solid','PUB')

            status=gnuplot_setstyle(ptr_gctrl,'pm3d')
            status=gnuplot_set(ptr_gctrl,'view map')
            status=gnuplot_unset(ptr_gctrl,'surface')
            status=gnuplot_set(ptr_gctrl,'style data pm3d')
            status=gnuplot_set(ptr_gctrl,'ticslevel 0')
            status=gnuplot_set(ptr_gctrl,'pm3d at b')
            if(present(logz)) then
               if(logz.eq.'log') status=gnuplot_setscale(ptr_gctrl,'zcb','LOG')
            end if
            status=gnuplot_set(ptr_gctrl,'palette')
            status=gnuplot_set(ptr_gctrl,'nokey')
            status=gnuplot_setaxislabel(ptr_gctrl,'x',trim(xlabel))
            status=gnuplot_setaxislabel(ptr_gctrl,'y',trim(ylabel))
            setplot3dsize : select case(numplots)
            case(1)
               status=gnuplot_set(ptr_gctrl,'size 1.0,1.0')
            case(2)
               status=gnuplot_set(ptr_gctrl,'size 1.0,2.0')
            case(3)
               status=gnuplot_set(ptr_gctrl,'size 1.0,3.0')
            case(4)
               status=gnuplot_set(ptr_gctrl,'size 2.0,2.0')
            end select setplot3dsize

            status=gnuplot_set(ptr_gctrl,"label """//trim(commontitle)//""" at screen 0.02,0.1,0.0 rotate font ""Helvetica,10""")

            if(numplots.gt.1) status=gnuplot_set(ptr_gctrl,'multiplot')
            status=gnuplot_set(ptr_gctrl,'origin 0.0,0.0')
            if(numplots.gt.1) status=gnuplot_set(ptr_gctrl,'size 0.9,0.9')
            status=gnuplot_settitle(ptr_gctrl,trim(zlabel))
            status=gnuplot_set(ptr_gctrl,'zrange [*:*]')
            status=gnuplot_plot3d(ptr_gctrl,ix,iy,xarray(1:ix),yarray(1:iy),zarray(1:ix,1:iy))
            i=i+1
            if(numplots.gt.1) then
               status=gnuplot_resetsession(ptr_gctrl)
               status=gnuplot_set(ptr_gctrl,'origin 0.0,1.0')
               status=gnuplot_settitle(ptr_gctrl,trim(zlabel1))
               status=gnuplot_set(ptr_gctrl,'zrange [*:*]')
               status=gnuplot_plot3d(ptr_gctrl,ix,iy,xarray(1:ix),yarray(1:iy),zarray1(1:ix,1:iy))
               i=i+1
            end if

            if(numplots.ge.3) then
               if(numplots.eq.3) status=gnuplot_set(ptr_gctrl,'origin 0.0,2.0')
               if(numplots.eq.4) status=gnuplot_set(ptr_gctrl,'origin 1.0,0.0')
               status=gnuplot_settitle(ptr_gctrl,trim(zlabel2))
               status=gnuplot_set(ptr_gctrl,'zrange [*:*]')
               status=gnuplot_plot3d(ptr_gctrl,ix,iy,xarray(1:ix),yarray(1:iy),zarray2(1:ix,1:iy))
               i=i+1
            end if

            if(numplots.eq.4) then
               status=gnuplot_set(ptr_gctrl,'origin 1.0,1.0')
               status=gnuplot_settitle(ptr_gctrl,trim(zlabel3))
               status=gnuplot_set(ptr_gctrl,'zrange [*:*]')
               status=gnuplot_plot3d(ptr_gctrl,ix,iy,xarray(1:ix),yarray(1:iy),zarray3(1:ix,1:iy))
               i=i+1
            end if

            if(numplots.gt.1) status=gnuplot_unset(ptr_gctrl,'multiplot')
            status=gnuplot_set(ptr_gctrl,'output')
            if(trim(pformat).eq.'PDF') status=gnuplot_genpdf(ptr_gctrl)
            status=gnuplot_close(ptr_gctrl)
          end function gnuplot_plot3dpanel

        end module gnuplot_module

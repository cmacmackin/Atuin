!  Copyright 2014 Chris MacMackin
!
!  This file is part of Atuin.
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  


!======================================================================!
!                     B E G I N    P R O G R A M :                     !
!                             M K H I S T                              !
!======================================================================!
!                                                                      !
!   AUTHOR:         Christopher MacMackin                              !
!   WRITTEN:        March, 2014                                        !
!   MODIFICATIONS:  None                                               !
!                                                                      !
!   PURPOSE:        Reads in dump files from Atuin simulations and     !
!                   converts them to a movie using Gnuplot and         !
!                   ImageMagick. Makes use of third party libraries to !
!                   interface with Gnuplot.                            !
!                                                                      !
!   EXTERNALS:      None                                               !
!                                                                      !
!----------------------------------------------------------------------!
PROGRAM mkhist
    USE datatypes
    USE gnuplot_module_data
    USE gnuplot_module
    USE physvals
    USE particle
    USE pdstrbtn
    IMPLICIT NONE
    
    CHARACTER(LEN=1)                    ::  remove
    CHARACTER(LEN=32)                   ::  arange,                    &
                                            dumpname,                  &
                                            listfile,                  &
                                            namefile,                  &
                                            plotname,                  &
                                            time,                      &
                                            xsun,                      &
                                            ysun
    CHARACTER(LEN=GP_CMD_SIZE)          ::  command
    INTEGER                             ::  argstat,                   &
                                            iter,                      &
                                            iter2,                     &
                                            iter3,                     &
                                            res,                       &
                                            retval
    REAL(8)                             ::  cpubeg = 0.d0,             &
                                            cpuend = 0.d0,             &
                                            cputot = 0.d0,             &
                                            xmin,                      &
                                            xmax,                      &
                                            ymin,                      &
                                            ymax
    REAL(8), ALLOCATABLE, DIMENSION(:)  ::  xvals,                     &
                                            yvals
    REAL(8), ALLOCATABLE, DIMENSION(:,:)::  zvals
    TYPE(gnuplot_ctrl), POINTER         ::  instance
!----------------------------------------------------------------------!

    command = 'splot "tmpmatrix.dat" matrix#'

    ! Start CPU timer
    CALL CPU_TIME (cpubeg  )

    ! Get name-list file from command-line arguments
    CALL GET_COMMAND_ARGUMENT(1,namefile,STATUS=argstat)
    IF ( argstat > 0 ) THEN
        WRITE( 6, 666)
        STOP
    ELSE IF ( argstat == -1 ) THEN
        WRITE( 6, 667)
        STOP
    END IF
    
    ! Load list file
    OPEN(10, FILE=namefile, STATUS='OLD', IOSTAT=argstat)
    IF ( argstat /= 0 ) THEN
        WRITE( 6, 668) namefile
        STOP
    END IF
    
    IF ( COMMAND_ARGUMENT_COUNT() < 6 ) WRITE( 6,1000)
    
    CALL GET_COMMAND_ARGUMENT(3,arange)
    READ(arange,    *) xmin
    CALL GET_COMMAND_ARGUMENT(4,arange)
    READ(arange,    *) xmax
    CALL GET_COMMAND_ARGUMENT(5,arange)
    READ(arange,    *) ymin
    CALL GET_COMMAND_ARGUMENT(6,arange)
    READ(arange,    *) ymax
        
    ! Iterate through the dump-files
    READ (10,    *, IOSTAT=argstat) dumpname
    DO WHILE ( argstat == 0 )
        CALL readin(filename=dumpname)
        WRITE(time,'(F12.3)') simtime/secsyear
        WRITE(xsun,   *) partlist(1)%xpos/auval
        WRITE(ysun,   *) partlist(1)%ypos/auval
        plotname = dumpname(1:INDEX(dumpname,'.',.TRUE.))//'png'
        
        ! Initialize gnuplot
        instance => gnuplot_init('')
        retval = gnuplot_hardcopy(instance,'PNG',plotname,             &
                 'size 500,500','PUB')
        retval = gnuplot_set(instance,'pm3d map')
        retval = gnuplot_set(instance,'pm3d interpolate 0,0')
        !retval = gnuplot_set(instance,'palette defined (0 "#000000", ' &
        !         //'1 "#0057FF")')
        retval = gnuplot_setaxisformat(instance,'x','')
        retval = gnuplot_setaxisformat(instance,'y','')
        retval = gnuplot_set(instance,'palette rgb 7,5,15')
        retval = gnuplot_set(instance,'size ratio -1')
        retval = gnuplot_set(instance,'lmargin screen 0.05')
        retval = gnuplot_set(instance,'bmargin screen 0.05')
        retval = gnuplot_set(instance,'rmargin screen 0.95')
        retval = gnuplot_set(instance,'tmargin screen 0.9')
        retval = gnuplot_unset(instance,'key')
        retval = gnuplot_unset(instance,'tics')
        retval = gnuplot_unset(instance,'colorbox')
        retval = gnuplot_settitle(instance,'Time = '//                 &
                 TRIM(ADJUSTL(time))//' yrs')
        
        ! Get positions of particles
        res = MIN(INT(SQRT(numnnull/2.d0)),100)
        ALLOCATE(xvals(res),yvals(res),zvals(res,res))
        FORALL(iter=1:res+1) xvals(iter) = xmin + (DBLE(iter-1) + 0.5) &
                                           *(xmax - xmin)/DBLE(res-1)
        FORALL(iter=1:res+1) yvals(iter) = ymin + (DBLE(iter-1) + 0.5) &
                                           *(ymax - ymin)/DBLE(res-1)
        zvals = 0
        DO iter = 2,ttlprtcl
            DO iter2 = 1,(res-1)
                DO iter3 = 1,(res-1)
                IF ( ( partlist(iter)%xpos/auval >= xvals(iter2) ).AND.&
                   ( partlist(iter)%xpos/auval < xvals(iter2 + 1) )    &
                   .AND. ( partlist(iter)%ypos/auval >= yvals(iter3) ) &
                   .AND. ( partlist(iter)%ypos/auval < yvals(iter3     &
                   + 1) ) ) THEN
                    zvals(iter2,iter3) = zvals(iter2,iter3) +          &
                                         partlist(iter)%mass
                END IF 
                END DO
            END DO
        END DO
        
        OPEN(12,FILE='tmpmatrix.dat',STATUS='UNKNOWN')
        WRITE(arange,    *) res
        DO iter = 1,res
            WRITE(12,   *) zvals(iter,:)
        END DO
        CLOSE(12)
        
        ! Produce plot
        retval = gnuplot_cmd(instance,command)
        retval = gnuplot_close(instance)
        
        ! Prepare for next iteration
        CALL SYSTEM('rm tmpmatrix.dat')
        DEALLOCATE(xvals,yvals,partlist,zvals)
        NULLIFY(instance)
        READ (10,   *, IOSTAT=argstat) dumpname
    END DO
    CLOSE(10)
    
    ! Create video from images
    CALL SYSTEM('convert -delay 10 *.png atuinvid.gif')
           
    ! Check if user wants to delete images
    IF ( COMMAND_ARGUMENT_COUNT() >= 2 ) THEN
        CALL GET_COMMAND_ARGUMENT(2,remove)
        IF ( remove /= 's' ) CALL SYSTEM('rm *.png')
    ELSE
        DO
            WRITE( 6,2000)
            READ ( 5,   *) remove
            SELECT CASE (remove)
                CASE ('y')
                    CALL SYSTEM('rm *.png')
                    EXIT
                CASE ('Y')
                    CALL SYSTEM('rm *.png')
                    EXIT
                CASE ('n')
                    EXIT
                CASE ('N')
                    EXIT
                CASE DEFAULT
                    WRITE( 6,2001)
            END SELECT
        END DO
    END IF
    
    ! End CPU time counter, and report time used.
    CALL CPU_TIME (cpuend  )
    cputot = cpuend - cpubeg
    WRITE( 6, 3000) cputot

!----------------------------------------------------------------------!
!                        Write format statements                       !
!----------------------------------------------------------------------!
     666 FORMAT('MKVIDEO : This program converts Atuin dump files ',   &
                'into a video of the simulation.',/,                   &
                'MKVIDEO : Please provide as arguments name of file ', &
                'with list of dump files ,',/,                         &
                'MKVIDEO: a specifier for whether to save image ',     &
                'files, and the x- and y-ranges.')
     667 FORMAT('MKVIDEO : List file name must be no more than 32 ',   &
                'characters.')
     668 FORMAT('MKVIDEO : Provided list file ',A32,' does ',/,        &
                'MKVIDEO : not exist.')
    1000 FORMAT('MKVIDEO : Wrong number of arguments.')
    2000 FORMAT('MKVIDEO : Delete individual images? (y/n)')
    2001 FORMAT('MKVIDEO : Respond with "y" or "n".')
    3000 FORMAT('MKVIDEO : ',/                                         &
               ,'MKVIDEO : Total cpu usage for this run is ',1pg12.5   &
               ,' seconds.')
!----------------------------------------------------------------------!

    STOP
END PROGRAM mkhist    
!======================================================================!
!                       E N D    P R O G R A M :                       !
!                             M K H I S T                              !
!======================================================================!

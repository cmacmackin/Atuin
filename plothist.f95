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
    IMPLICIT NONE
    
    CHARACTER(LEN=1)                    ::  remove
    CHARACTER(LEN=32)                   ::  matrix,                    &
                                            namefile,                  &
                                            plotname,                  &
                                            time
    CHARACTER(LEN=GP_CMD_SIZE)          ::  command
    INTEGER                             ::  argstat,                   &
                                            retval
    REAL(8)                             ::  cpubeg = 0.d0,             &
                                            cpuend = 0.d0,             &
                                            cputot = 0.d0
    TYPE(gnuplot_ctrl), POINTER         ::  instance
!----------------------------------------------------------------------!

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
    
    ! Iterate through the dump-files
    READ (10, '(A32,A12)', IOSTAT=argstat) matrix, time
    DO WHILE ( argstat == 0 )
        plotname = matrix(1:INDEX(matrix,'.',.TRUE.))//'png'
        
        ! Initialize gnuplot
        instance => gnuplot_init('')
        retval = gnuplot_hardcopy(instance,'PNG',plotname,             &
                 'size 500,500','PUB')
        retval = gnuplot_set(instance,'pm3d map')
        retval = gnuplot_set(instance,'pm3d interpolate 0,0')
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
        
        ! Produce plot
        command = 'splot "'//TRIM(matrix)//'" matrix#'
        retval = gnuplot_cmd(instance,command)
        retval = gnuplot_close(instance)
        
        ! Prepare for next iteration
        NULLIFY(instance)
        READ (10, '(A32,A12)', IOSTAT=argstat) matrix, time
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

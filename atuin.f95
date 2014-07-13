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
!                              A T U I N                               !
!======================================================================!
!                                                                      !
!   AUTHOR:         Christopher MacMackin                              !
!   WRITTEN:        February, 2014                                     !
!   MODIFICATIONS:  None                                               !
!                                                                      !
!   PURPOSE:        Solves the n-body problem of the formation of      !
!                   rocky planets, using a treecode algorithm. This    !
!                   program is named after the turtle Great A'Tuin in  !
!                   Terry Pratchett's *Discworld* novels, which        !
!                   carries the world on its back.                     !
!                                                                      !
!   EXTERNALS:      None                                               !
!                                                                      !
!----------------------------------------------------------------------!
PROGRAM atuin
    USE forces
    USE leap2d
    USE particle
    USE pdstrbtn
    USE physvals
    USE tree2d
    IMPLICIT NONE

    ! Declare variables
    CHARACTER(LEN=9)    ::  iterchar = ''
    CHARACTER(LEN=32)   ::  dumpname,                                  &
                            infile,                                    &
                            namefile,                                  &
                            outfile
    INTEGER             ::  argstat,                                   &
                            dumpfreq,                                  &
                            iter,                                      &
                            numpart,                                   &
                            startnum
    LOGICAL             ::  new
    REAL(8)             ::  cpubeg,                                    &
                            cpuend,                                    &
                            cputot,                                    &
                            density,                                   &
                            discmass,                                  &
                            massdev,                                   &
                            outrad,                                    &
                            inrad,                                     &
                            stardnst,                                  &
                            starmass,                                  &
                            ttltime,                                   &
                            veldev

    ! Specify name-lists
    NAMELIST /simspecs/ density , discmass, infile  , inrad   ,        &
                        massdev , new     , numpart , outrad  ,        &
                        stardnst, starmass, ttltime , veldev
    NAMELIST /outspecs/ dumpfreq, dumpname, outfile , startnum
!----------------------------------------------------------------------!

    ! Start CPU timer
    CALL CPU_TIME(cpubeg)

    ! Get name-list file from command-line arguments
    CALL GET_COMMAND_ARGUMENT(1,namefile,STATUS=argstat)
    IF ( argstat > 0 ) THEN
        WRITE( 6, 666)
        STOP
    ELSE IF ( argstat == -1 ) THEN
        WRITE( 6, 667)
        STOP
    END IF

    ! Load name-lists
    OPEN(10, FILE=namefile, STATUS='OLD', RECL=80, DELIM='APOSTROPHE', &
         IOSTAT=argstat)
    IF ( argstat /= 0 ) THEN
        WRITE( 6, 668) namefile
        STOP
    END IF
    READ (10,NML=simspecs,IOSTAT=argstat)
    READ (10,NML=outspecs,IOSTAT=argstat)
    CLOSE(10)
    IF ( argstat /= 0 ) THEN
        WRITE( 6, 669) namefile
        STOP
    END IF

    ! Set up particle list
    WRITE( 6,1000)
    IF ( new .EQV. .TRUE. ) THEN
        WRITE( 6,1010) numpart
        CALL mkdisc(density,discmass,inrad,massdev,numpart,outrad,     &
                    stardnst,starmass,veldev)
        simtime = 0.d0
        iter = 0
    ELSE
        WRITE( 6,1020) infile
        CALL readin(infile)
        iter = startnum
    END IF
    
    ! Produce the first dump of particle distribution
    WRITE(iterchar,'(I9.9)') iter
    CALL writeout(TRIM(dumpname)//'-'//iterchar//'.out')
    
    ! Produce header and first entry of file tracking simulation
    OPEN(12, FILE=outfile, STATUS='UNKNOWN')
    WRITE(12,1030) 
    WRITE(12,1035) iter, simtime, deltat, xcom, ycom, energy, angmom,  &
                   xlinmom, ylinmom, numnnull

    ! Integrate forward until have reached end of simulation time or 
    ! have only one particle remaining
    WRITE( 6,1040)
    DO WHILE ( ( simtime < ttltime ) .AND. ( numnnull > 1 ) )
        CALL updtall
        CALL allaccel(gravity)
        CALL snglleap
        iter = iter + 1
        simtime = simtime + deltat
        IF ( MOD(iter,MAX(1,dumpfreq/1000)) == 0 ) WRITE(12,1035) iter,&
            simtime, deltat, xcom, ycom, energy, angmom, xlinmom,      &
            ylinmom, numnnull
        
        ! Dump particle distribution (if necessary)
        IF ( MOD(iter,dumpfreq) == 0 ) THEN
            WRITE(iterchar,'(I9.9)') iter
            CALL writeout(TRIM(dumpname)//'-'//iterchar//'.out')
            WRITE( 6,1045) iter, simtime/secsyear
        END IF
    END DO
    close(9)
    ! Prepare to terminate program
    DEALLOCATE(partlist)
    WRITE( 6,1050) simtime/secsyear, iter, nchngdt
    
    ! End CPU time counter, and report time used.
    CALL CPU_TIME(cpuend)
    cputot = cpuend - cpubeg
    WRITE( 6,3000) cputot

!----------------------------------------------------------------------!
!                        Write format statements                       !
!----------------------------------------------------------------------!
     666 FORMAT('ATUIN   : This is an n-body solver. Please provide ', &
                'as an argument the ',/,                               &
                'ATUIN   : name of a file containing a name-list with',&
                'parameters for the ',/,                               &
                'ATUIN   : simulation.')
     667 FORMAT('ATUIN   : Name-list file name must be no more than ', &
                '32 characters.')
     668 FORMAT('ATUIN   : Provided name-list file ',A32,' does ',/,   &
                'ATUIN   : not exist.')
     669 FORMAT('ATUIN   : Provided name-list file ',A32,' not ',/,    &
                'ATUIN   : properly formatted.')
    1000 FORMAT('ATUIN   : ==========================================',&
                '================',/,                                  &
                'ATUIN   :  ___  ______ __ __ __ __  __ ',/,           &
                'ATUIN   : // \\ | || | || || || ||\ || ',/,           &
                'ATUIN   : ||=||   ||   || || || ||\\|| ',/,           &
                'ATUIN   : || ||   ||   \\_// || || \|| ',/,           &
                'ATUIN   :                              ',/,           &
                'ATUIN   : __  __  ____    ___   ____   _  _  ',/,     &
                'ATUIN   : ||\ ||  || ))  // \\  || \\  \\//  ',/,     &
                'ATUIN   : ||\\||  ||=)  ((   )) ||  ))  )/   ',/,     &
                'ATUIN   : || \||  ||_))  \\_//  ||_//  //    ',/,     &
                'ATUIN   :                                    ',/,     &
                'ATUIN   :  __  __ ___  ___ __ __ __     ___  ______ ',&
                '  ___   ____ ',/,                                     &
                'ATUIN   : (( \ || ||\\//|| || || ||    // \\ | || | ',&
                ' // \\  || \\',/,                                     &
                'ATUIN   :  \\  || || \/ || || || ||    ||=||   ||   ',&
                '((   )) ||_//',/,                                     &
                'ATUIN   : \_)) || ||    || \\_// ||__| || ||   ||   ',&
                ' \\_//  || \\',/,                                     &
                "ATUIN   :                      _,.---.---.---.--..",  &
                "_",/,                                                 &
                "ATUIN   :                  _.-' `--.`---.`---'-. _,`",&
                "--.._",/,                                             &
                "ATUIN   :                 /`--._ .'.     `.     `,`-",&
                ".`-._\",/,                                            &
                "ATUIN   :                ||   \  `.`---.__`__..-`. ,",&
                "'`-._/",/,                                            &
                "ATUIN   :           _  ,`\ `-._\   \    `.    `_.-`-",&
                "._,``-.",/,                                           &
                "ATUIN   :        ,`   `-_ \/ `-.`--.\    _\_.-'\__.-",&
                "`-.`-._`.",/,                                         &
                "ATUIN   :       (_.o> ,--. `._/'--.-`,--`  \_.-'    ",&
                "   \`-._ \",/,                                        &
                "ATUIN   :        `---'    `._ `---._/__,----`       ",&
                "    `-. `-\",/,                                       &
                "ATUIN   :                  /_, ,  _..-'             ",&
                "       `-._\",/,                                      &
                "ATUIN   :                  \_, \/ ._(",/,             &
                "ATUIN   :                   \_, \/ ._\",/,            &
                "ATUIN   :                    `._,\/ ._\",/,           &
                "ATUIN   :                      `._// ./`-._",/,       &
                "ATUIN   :                        `-._-_-_.-'",/,      &
                'ATUIN   : ==========================================',&
                '================')
    1010 FORMAT('ATUIN   : Initializing protoplanetary disc of ',I8,   &
                ' particles.')
    1020 FORMAT('ATUIN   : Reading particle distribution from file ',/,&
                'ATUIN   : ',A32,'.')
    1030 FORMAT('#ATUIN: step, time, delta t, CofM x, CofM y, ',       &
                'energy, ang momentum, lin momenutm x, lin momentum ', &
                'y, non-null')
    1035 FORMAT(I12,1P8G24.15,I12)
    1040 FORMAT('ATUIN   : Beginning integration...')
    1045 FORMAT('ATUIN   : On step ',I9,' at time ',1PG12.5,' years.')
    1050 FORMAT('ATUIN   : ',/,                                        &
                'ATUIN   : Simulation ran for ',1PG12.5,' years.',/,   &
                'ATUIN   : Required ',I9,' time steps, with ',I5,      &
                ' adjustments.' )
    3000 FORMAT('ATUIN   : ',/,                                        &
                'ATUIN   : Total cpu usage for this run is ',1PG12.5,  &
                ' seconds.')
!----------------------------------------------------------------------!

    STOP
END PROGRAM atuin    
!======================================================================!
!                       E N D    P R O G R A M :                       !
!                              A T U I N                               !
!======================================================================!

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
!                      B E G I N    M O D U L E :                      !
!                           P D S T R B T N                            !
!======================================================================!
!                                                                      !
!   PURPOSE:    Holds and processes the particles used in an n-body    !
!               simulation. Contains the information regarding the     !
!               state of the simulated system.                         !
!   CONTAINS:                                                          !
!   EXTERNALS:  particle (module), physvals (module), randgen (module) !
!                                                                      !
!----------------------------------------------------------------------!
MODULE pdstrbtn
    USE particle
    USE physvals
    USE randgen
    IMPLICIT NONE

    ! Variable declarations
    INTEGER                                         ::  ionum = 2,     &
                                                        nchngdt = -1,  &
                                                        numnnull,      &
                                                        ttlprtcl
    REAL(8)                                         ::  angmom,        &
                                                        deltat =       &
                                                        HUGE(deltat),  &
                                                        energy,        &
                                                        simtime,       &
                                                        totmass,       &
                                                        xcom,          &
                                                        xlargest,      &
                                                        xlinmom,       &
                                                        xsmallst,      &
                                                        ycom,          &
                                                        ylargest,      &
                                                        ylinmom,       &
                                                        ysmallst
    TYPE(prtcl), DIMENSION(:), ALLOCATABLE, TARGET  ::  partlist
    
!----------------------------------------------------------------------!
CONTAINS

    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                           M K 2 B O D Y                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        March, 2014                                    !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Initializes the particle distribution to be a  !
    !                   system of two similar massed objects (0.8 and  !
    !                   0.6 solar masses) orbiting each other with an  !
    !                   apocentre of 1AU and eccentricity of 0.6.      !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE mk2body
        IMPLICIT NONE

        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) THEN
            DEALLOCATE(partlist)
            WRITE(6,1000)
        END IF
                
        ALLOCATE(partlist(2))
        ttlprtcl = 2
        
        ! Initialize first particle
        partlist(1)%mass = 8.d-1*m_sol
        partlist(1)%density = d_sol
        partlist(1)%xpos = -3.d0/2.8d1*auval
        partlist(1)%yvel = -3.d0/7.d0*SQRT(8.d0*gconst*1.4d0*m_sol/    &
                            (1.25d0*auval))*COS(pi*163859.15724971273/ &
                            (0.289318781*secsyear))
        partlist(1)%xvel = -3.d0/7.d0*SQRT(8.d0*gconst*1.4d0*m_sol/    &
                            (1.25d0*auval))*SIN(pi*163859.15724971273/ &
                            (0.289318781*secsyear))
        
        ! Initialize second particle
        partlist(2)%mass = 6.d-1*m_sol
        partlist(2)%density = d_sol
        partlist(2)%xpos = 4.d0/2.8d1*auval
        partlist(2)%yvel = 4.d0/7.d0*SQRT(8.d0*gconst*1.4d0*m_sol/     &
                            (1.25d0*auval))*COS(pi*163859.15724971273/ &
                            (0.289318781*secsyear))
        partlist(2)%xvel = 4.d0/7.d0*SQRT(8.d0*gconst*1.4d0*m_sol/     &
                            (1.25d0*auval))*SIN(pi*163859.15724971273/ &
                            (0.289318781*secsyear))
        
        numnnull = 2
        
        xsmallst = MIN(partlist(1)%xpos,partlist(2)%xpos)
        ysmallst = MIN(partlist(1)%ypos,partlist(2)%ypos)
        xlargest = MAX(partlist(1)%xpos,partlist(2)%xpos) + 1.d0
        ylargest = MAX(partlist(1)%ypos,partlist(2)%ypos) + 1.d0
        simtime = 0.d0

    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MK2BODY : WARNING: List of particles already',    &
                    ' present and was deallocated',/,                  &
                    'MK2BODY : prior to creating new one.')
    !------------------------------------------------------------------!
        
        RETURN
    END SUBROUTINE mk2body
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                           M K 2 B O D Y                          !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                           M K 4 B O D Y                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        March, 2014                                    !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Initializes the particle distribution to be a  !
    !                   system of four particles of solar mass. The    !
    !                   particles are in pairs, separated by 2AU,      !
    !                   which are in a circular orbit about each       !
    !                   other. The particles in each pair are          !
    !                   separated by 0.1AU and are also in circular    !
    !                   orbits about each other.                       !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE mk4body
        IMPLICIT NONE

        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) THEN
            DEALLOCATE(partlist)
            WRITE(6,1000)
        END IF
                
        ALLOCATE(partlist(4))
        ttlprtcl = 4
        
        ! Initialize first particle
        partlist(1)%mass = m_sol
        partlist(1)%density = d_sol
        partlist(1)%xpos = 1.d0*auval
        partlist(1)%ypos = -5.d-2*auval
        partlist(1)%yvel = 5.d-1*SQRT(2.d0*gconst*m_sol/auval)
        partlist(1)%xvel = 5.d-1*SQRT(2.d1*gconst*m_sol/auval)
        
        ! Initialize second particle
        partlist(2)%mass = m_sol
        partlist(2)%density = d_sol
        partlist(2)%xpos = 1.d0*auval
        partlist(2)%ypos = 5.d-2*auval
        partlist(2)%yvel = 5.d-1*SQRT(2.d0*gconst*m_sol/auval)
        partlist(2)%xvel = -5.d-1*SQRT(2.d1*gconst*m_sol/auval)
        
        ! Initialize third particle
        partlist(3)%mass = m_sol
        partlist(3)%density = d_sol
        partlist(3)%xpos = -1.d0*auval
        partlist(3)%ypos = -5.d-2*auval
        partlist(3)%yvel = -5.d-1*SQRT(2.d0*gconst*m_sol/auval)
        partlist(3)%xvel = 5.d-1*SQRT(2.d1*gconst*m_sol/auval)
        
        ! Initialize fourth particle
        partlist(4)%mass = m_sol
        partlist(4)%density = d_sol
        partlist(4)%xpos = -1.d0*auval
        partlist(4)%ypos = 5.d-2*auval
        partlist(4)%yvel = -5.d-1*SQRT(2.d0*gconst*m_sol/auval)
        partlist(4)%xvel = -5.d-1*SQRT(2.d1*gconst*m_sol/auval)
        
        numnnull = 4
        
        xsmallst = -1.5d1*auval
        ysmallst = -1.5d1*auval
        xlargest = 1.5d1*auval
        ylargest = 1.5d1*auval
        simtime = 0.d0

    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MK2BODY : WARNING: List of particles already',    &
                    ' present and was deallocated',/,                  &
                    'MK2BODY : prior to creating new one.')
    !------------------------------------------------------------------!
        
        RETURN
    END SUBROUTINE mk4body
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                           M K 4 B O D Y                          !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                            M K D I S C                           !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Initializes the particle distribution to be a  !
    !                   protoplanetary disc with the desired           !
    !                   properties.                                    !
    !                                                                  !
    !   ARGUMENTS:     *discdnst, the density of disc particles.       !
    !                  *discmass, the total mass (not including star)  !
    !                             of the protoplanetary disc. Note     !
    !                             that the final disc mass will not    !
    !                             exactly match this value, but should !
    !                             be quite close (especially for large !
    !                             numpart).                            !
    !                  *inner   , the inner radius of the disc.        !
    !                  *massdev , standard deviation in the mass of    !
    !                             the disc particles, as fraction of   !
    !                             average particle mass.               !
    !                  *numpart , number of disc particles.            !
    !                  *outer   , outer radius of the disc.            !
    !                  *stardnst, density of the central star.         !
    !                  *starmass, mass of the star which the           !
    !                             protoplanetary disc orbits.          !
    !                  *veldev  , standard deviation in particle       !
    !                             velocity, as fraction of velocity    !
    !                             that would result in circular orbit. !
    !                                                                  !
    !                   Arguments marked with an asterisk (*) are ones !
    !                   for which the user must specify a value when   !
    !                   the subroutine is called. Those without an     !
    !                   asterisk (*) are used only to return a value   !
    !                   to the caller.                                 !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE mkdisc  ( discdnst, discmass, inner   , massdev ,       &
                         numpart , outer   , stardnst, starmass,       &
                         veldev   )
        IMPLICIT NONE
        
        ! Input and output variables:
        INTEGER, INTENT(IN) ::  numpart
        REAL(8), INTENT(IN) ::  discdnst,                              &
                                discmass,                              &
                                inner,                                 &
                                massdev,                               &
                                outer,                                 &
                                stardnst,                              &
                                starmass,                              &
                                veldev
        
        ! Other variables:
        INTEGER ::  iter
        REAL(8) ::  angpos,                                            &
                    avemass,                                           &
                    orbvel,                                            &
                    radpos
    !------------------------------------------------------------------!

        ! Initialize and allocate variables
        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) THEN
            DEALLOCATE(partlist)
            WRITE(6,1000)
        END IF
        ALLOCATE(partlist(numpart + 1))
        numnnull = numpart + 1
        ttlprtcl = numnnull
        avemass = discmass/DBLE(numpart)
        xsmallst = -1.1d0*outer
        ysmallst = -1.1d0*outer
        xlargest = 1.1d0*inner
        ylargest = 1.1d0*inner
        simtime = 0.d0
        
        ! Initialize star
        partlist(1)%ptype = 2
        partlist(1)%mass = starmass
        partlist(1)%density = stardnst
        
        ! Initialize random number generator. Use numpart as initial
        ! seed (as it is an integer that is unlikely to be zero).
        radpos = unirand(numpart)
        
        ! Initialize particles
        DO iter = 2,(numpart + 1)
            ! Randomly select properties
            radpos = ((SQRT(outer) - SQRT(inner))*unirand() +          &
                     SQRT(inner))**2.d0
            angpos = unirand()*2.d0*pi
            orbvel = SQRT(gconst*starmass/radpos)*(1.d0 + veldev*      &
                     normrand())
            
            ! Assign these properties
            partlist(iter)%mass = ABS(avemass*(1.d0 + massdev*         &
                                  normrand()))
            partlist(iter)%density = discdnst
            partlist(iter)%xpos = radpos*COS(angpos)
            partlist(iter)%ypos = radpos*SIN(angpos)
            partlist(iter)%xvel = -orbvel*SIN(angpos)
            partlist(iter)%yvel = orbvel*COS(angpos)

            ! Add opposite momentum to sun, so that total momentum 
            ! remains zero
            partlist(1)%xvel = partlist(1)%xvel - partlist(iter)%xvel* &
                               partlist(iter)%mass/starmass
            partlist(1)%yvel = partlist(1)%yvel - partlist(iter)%yvel* &
                               partlist(iter)%mass/starmass
        END DO

    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MKDISC  : WARNING: List of particles already',    &
                    ' present and was deallocated',/,                  &
                    'MKDISC  : prior to creating new one.')
    !------------------------------------------------------------------!
                
        RETURN
    END SUBROUTINE mkdisc
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                            M K D I S C                           !
    !==================================================================!


    
    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                          M K E A R T H 1                         !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Initializes the particle distribution to be a  !
    !                   planet the size of the earth orbiting a star   !
    !                   with the properties of the sun. Meant for an   !
    !                   explicit integrator.                           !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE mkearth1
        IMPLICIT NONE
        
        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) THEN
            DEALLOCATE(partlist)
            WRITE(6,1000)
        END IF
                
        ALLOCATE(partlist(2))
        ttlprtcl = 2
        
        ! Initialize Sun
        partlist(1)%ptype = 2
        partlist(1)%mass = m_sol
        partlist(1)%lumin = l_sol
        partlist(1)%density = d_sol
        
        ! Initialize Earth
        partlist(2)%mass = m_earth
        partlist(2)%density = d_earth
        partlist(2)%ypos = -auval
        partlist(2)%xvel = 2.d0*pi*auval/secsyear
        
        numnnull = 2
        
        xsmallst = MIN(partlist(1)%xpos,partlist(2)%xpos)
        ysmallst = MIN(partlist(1)%ypos,partlist(2)%ypos)
        xlargest = MAX(partlist(1)%xpos,partlist(2)%xpos) + 1.d0
        ylargest = MAX(partlist(1)%ypos,partlist(2)%ypos) + 1.d0
        simtime = 0.d0
        
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MKEARTH1: WARNING: List of particles already',    &
                    ' present and was deallocated',/,                  &
                    'MKEARTH1: prior to creating new one.')
    !------------------------------------------------------------------!
        
        RETURN
    END SUBROUTINE mkearth1
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                          M K E A R T H 1                         !
    !==================================================================!

    

    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                          M K E A R T H 2                         !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Initializes the particle distribution to be a  !
    !                   planet the size of the earth orbiting a star   !
    !                   with the properties of the sun. Meant for a    !
    !                   leapfrog integrator.                           !
    !                                                                  !
    !   ARGUMENTS:     *deltat  , the size of the time-step being      !
    !                             used for the simulation.             !
    !                                                                  !
    !                   Arguments marked with an asterisk (*) are ones !
    !                   for which the user must specify a value when   !
    !                   the subroutine is called. Those without an     !
    !                   asterisk (*) are used only to return a value   !
    !                   to the caller.                                 !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE mkearth2
        IMPLICIT NONE
        
        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) THEN
            DEALLOCATE(partlist)
            WRITE(6,1000)
        END IF
                
        ALLOCATE(partlist(2))
        ttlprtcl = 2
        
        ! Initialize Sun
        partlist(1)%ptype = 2
        partlist(1)%mass = m_sol
        partlist(1)%lumin = l_sol
        partlist(1)%density = d_sol
        partlist(1)%xvel = -2.d0*pi*auval/secsyear*m_earth/(m_earth +  &
                           m_sol)
        
        ! Initialize Earth
        partlist(2)%mass = m_earth
        partlist(2)%density = d_earth
        partlist(2)%ypos = -auval
        partlist(2)%xvel = 2.d0*pi*auval/secsyear
        
        numnnull = 2
        
        xsmallst = MIN(partlist(1)%xpos,partlist(2)%xpos)
        ysmallst = MIN(partlist(1)%ypos,partlist(2)%ypos)
        xlargest = MAX(partlist(1)%xpos,partlist(2)%xpos) + 1.d0
        ylargest = MAX(partlist(1)%ypos,partlist(2)%ypos) + 1.d0
        simtime = 0.d0

    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MKEARTH2: WARNING: List of particles already',    &
                    ' present and was deallocated',/,                  &
                    'MKEARTH2: prior to creating new one.')
    !------------------------------------------------------------------!
        
        RETURN
    END SUBROUTINE mkearth2
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                          M K E A R T H 2                         !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                           M K E L L I P                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Initializes the particle distribution to be a  !
    !                   planet the size of the earth orbiting  in an   !
    !                   elliptical orbit around a star with the        !
    !                   properties of the sun. Meant for a leapfrog    !
    !                   integrator.                                    !
    !                                                                  !
    !   ARGUMENTS:     *deltat  , the size of the time-step being      !
    !                             used for the simulation.             !
    !                                                                  !
    !                   Arguments marked with an asterisk (*) are ones !
    !                   for which the user must specify a value when   !
    !                   the subroutine is called. Those without an     !
    !                   asterisk (*) are used only to return a value   !
    !                   to the caller.                                 !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE mkellip
        IMPLICIT NONE

        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) THEN
            DEALLOCATE(partlist)
            WRITE(6,1000)
        END IF
                
        ALLOCATE(partlist(2))
        ttlprtcl = 2
        
        ! Initialize Sun
        partlist(1)%ptype = 2
        partlist(1)%mass = m_sol
        partlist(1)%lumin = l_sol
        partlist(1)%density = d_sol
        partlist(1)%xvel = -1.2d0*2.d0*pi*auval/secsyear * m_earth/m_sol
        
        ! Initialize Earth
        partlist(2)%mass = m_earth
        partlist(2)%density = d_earth
        partlist(2)%ypos = -auval
        partlist(2)%xvel = 1.2d0*2.d0*pi*auval/secsyear
        !partlist(2)%yvel = -1.2d0*2.d0*pi*auval/secsyear

        numnnull = 2
        
        xsmallst = MIN(partlist(1)%xpos,partlist(2)%xpos)
        ysmallst = MIN(partlist(1)%ypos,partlist(2)%ypos)
        xlargest = MAX(partlist(1)%xpos,partlist(2)%xpos) + 1.d0
        ylargest = MAX(partlist(1)%ypos,partlist(2)%ypos) + 1.d0
        simtime = 0.d0

    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MKELLIP : WARNING: List of particles already',    &
                    ' present and was deallocated',/,                  &
                    'MKELLIP : prior to creating new one.')
    !------------------------------------------------------------------!
        
        RETURN
    END SUBROUTINE mkellip
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                           M K E L L I P                          !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                           R E A D I N                            !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Reads in a list of particles from the file     !
    !                   specified. These are placed in the provided    !
    !                   array, or the partlist array if none is        !
    !                   provided in the arguments.                     !
    !                                                                  !
    !   ARGUMENTS:     *filename, the name of the file from which to   !
    !                             read.                                !
    !                   filetime, the time in the simulation at which  !
    !                             the file being read in was           !
    !                             originally written out at.           !
    !                             (OPTIONAL)                           !
    !                   list    , an allocatable array into which the  !
    !                             particles from the file will be      !
    !                             read. (OPTIONAL)                     !
    !                                                                  !
    !                   Arguments marked with an asterisk (*) are ones !
    !                   for which the user must specify a value when   !
    !                   the subroutine is called. Those without an     !
    !                   asterisk (*) are used only to return a value   !
    !                   to the caller.                                 !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE readin   ( filename )
        IMPLICIT NONE
        
        ! Input and output variables:
        CHARACTER(LEN=*), INTENT(IN)    ::  filename
        
        ! Other variables:
        INTEGER ::  iter
    !------------------------------------------------------------------!
        
        OPEN(ionum, FILE=filename, STATUS='UNKNOWN', FORM='UNFORMATTED')
        
        ! Read number of particles from file
        READ (ionum) ttlprtcl
        READ (ionum) simtime
                
        IF ( ALLOCATED(partlist) .EQV. .TRUE. ) DEALLOCATE(partlist)
        ALLOCATE(partlist(ttlprtcl)) 
        numnnull = 0
        xlargest = -HUGE(xlargest)
        xsmallst = HUGE(xlargest)
        ylargest = -HUGE(ylargest)
        ysmallst = HUGE(ylargest)
        
        ! Read in all particles from file
        DO iter = 1, ttlprtcl
            READ (ionum) partlist(iter)%ptype , partlist(iter)%lumin , &
                         partlist(iter)%mass  , partlist(iter)%xaccel, &
                         partlist(iter)%xpos  , partlist(iter)%xvel  , &
                         partlist(iter)%yaccel, partlist(iter)%ypos  , &
                         partlist(iter)%yvel
            IF ( partlist(iter)%ptype > 0 ) THEN
                numnnull = numnnull + 1
                xlargest = MAX(xlargest,partlist(iter)%xpos)
                xsmallst = MIN(xlargest,partlist(iter)%xpos)
                ylargest = MAX(ylargest,partlist(iter)%ypos)
                ysmallst = MIN(ylargest,partlist(iter)%ypos)
            END IF
        END DO
                
        ! Close file
        CLOSE(ionum)
        
        RETURN
    END SUBROUTINE readin
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                           R E A D I N                            !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         W R I T E O U T                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Writes the provided list of particles (or, if  !
    !                   none provided, the list of all particles) to a !
    !                   binary file. These can be read in again later. !
    !                                                                  !
    !   ARGUMENTS:     *filename, the name of the file to which to     !
    !                             write.                               !
    !                  *list    , the list of particles to write to a  !
    !                             file. (OPTIONAL)                     !
    !                                                                  !
    !                   Arguments marked with an asterisk (*) are ones !
    !                   for which the user must specify a value when   !
    !                   the subroutine is called. Those without an     !
    !                   asterisk (*) are used only to return a value   !
    !                   to the caller.                                 !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE writeout ( filename, list     )
        IMPLICIT NONE
        
        ! Input and output variables:
        CHARACTER(LEN=*), INTENT(IN)                    ::  filename
        TYPE(prtcl), INTENT(IN), DIMENSION(:), OPTIONAL ::  list

        ! Other variables:
        INTEGER                                 ::  iter,              &
                                                    partnum
        TYPE(prtcl), ALLOCATABLE, DIMENSION(:)  ::  printlst
    !------------------------------------------------------------------!
        
        IF ( PRESENT(list) .EQV. .TRUE. ) THEN
            ALLOCATE(printlst(SIZE(list)))
            printlst = list
            partnum = SIZE(list)
        ELSE
            ALLOCATE(printlst(ttlprtcl))
            printlst = partlist
            partnum = ttlprtcl
        END IF
        
        OPEN(ionum, FILE=filename, STATUS='UNKNOWN', FORM='UNFORMATTED')
        
        ! Write number of particles in file
        WRITE(ionum) SIZE(printlst)
        WRITE(ionum) simtime
        
        ! Loop through particle list and print
        DO iter = 1, partnum
            WRITE(ionum) partlist(iter)%ptype , partlist(iter)%lumin , &
                         partlist(iter)%mass  , partlist(iter)%xaccel, &
                         partlist(iter)%xpos  , partlist(iter)%xvel  , &
                         partlist(iter)%yaccel, partlist(iter)%ypos  , &
                         partlist(iter)%yvel
        END DO
        
        ! Write end of file identifier
        CLOSE(ionum)
        
        RETURN
    END SUBROUTINE writeout
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                         W R I T E O U T                          !
    !==================================================================!

END MODULE pdstrbtn
!======================================================================!
!                        E N D    M O D U L E :                        !
!                           P D S T R B T N                            !
!======================================================================!

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
!                             L E A P 2 D                              !
!======================================================================!
!                                                                      !
!   PURPOSE:    Uses a leapfrog algorithm to calculate the motion of a !
!               Newtonian n-body system. Designed for a 2D system (ie: !
!               all particles constrained to a single plain). Contains !
!               only the subroutines and data needed to take           !
!               simulation forward one step, none of the information   !
!               about the state of the simulated system.               !
!   CONTAINS:   snglleap (subroutine)                                  !
!   EXTERNALS:  particle (module), pdstrbtn (module)                   !
!                                                                      !
!----------------------------------------------------------------------!
MODULE leap2d
    USE particle
    USE pdstrbtn
    IMPLICIT NONE

    ! Variable declarations
    
!----------------------------------------------------------------------!
CONTAINS

    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         S N G L L E A P                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        A simple, single leapfrog algorithm, which can !
    !                   be used when none of the forces acting on a    !
    !                   particle are functions of the particle's       !
    !                   velocity.                                      !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE snglleap
        IMPLICIT NONE
        
        ! Variables declarations:
        INTEGER                             ::  iter
    !------------------------------------------------------------------!
        
        ! Initialize values:
        xsmallst = HUGE(xsmallst)
        xlargest = -HUGE(xlargest)
        ysmallst = HUGE(ysmallst)
        ylargest = -HUGE(ylargest)
        
        ! Iterate through particles, and integrate forward one time-step
        DO iter = 1, ttlprtcl
            IF ( partlist(iter)%ptype > 0 ) THEN
                partlist(iter)%xvel = partlist(iter)%xvel +            &
                                      partlist(iter)%xaccel*deltat
                partlist(iter)%yvel = partlist(iter)%yvel +            &
                                      partlist(iter)%yaccel*deltat
                partlist(iter)%xpos = partlist(iter)%xpos +            &
                                      partlist(iter)%xvel*deltat
                partlist(iter)%ypos = partlist(iter)%ypos +            &
                                      partlist(iter)%yvel*deltat
                xsmallst = MIN(xsmallst,partlist(iter)%xpos)
                xlargest = MAX(xlargest,partlist(iter)%xpos)
                ysmallst = MIN(ysmallst,partlist(iter)%ypos)
                ylargest = MAX(ylargest,partlist(iter)%ypos) 
            END IF
        END DO
        
        RETURN
    END SUBROUTINE snglleap
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                         S N G L L E A P                          !
    !==================================================================!

END MODULE leap2d
!======================================================================!
!                        E N D    M O D U L E :                        !
!                             L E A P 2 D                              !
!======================================================================!

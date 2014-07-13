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
!                           P A R T I C L E                            !
!======================================================================!
!                                                                      !
!   PURPOSE:    Provides derived type "particle" and associated        !
!               methods for use in charge-free 2D n-body simulations   !
!   CONTAINS:   getdist (function), mergeprt (subroutine), prtcl       !
!               (derived type)                                         !
!   EXTERNALS:  None                                                   !
!                                                                      !
!----------------------------------------------------------------------!
MODULE particle
    USE physvals
    IMPLICIT NONE

    !==================================================================!
    !                      B E G I N    T Y P E :                      !
    !                            P R T C L                             !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Describes a particle in an n-body simulation.  !
    !                   Particle may have type, -1 = fictitious (used  !
    !                   for tree code), 0 = null, 1 = dust, 2 = star.  !
    !                   All values in SI/MKS, unless otherwise         !
    !                   indicated.                                     !
    !                                                                  !
    !------------------------------------------------------------------!
    TYPE prtcl
        INTEGER ::  ptype = 1
        REAL(8) ::  density = 1.d3,                                    &
                    lumin = 0.d0,                                      &
                    mass = 0.d0,                                       &
                    xaccel = 0.d0,                                     &
                    xpos = 0.d0,                                       &
                    xvel = 0.d0,                                       &
                    yaccel = 0.d0,                                     &
                    ypos = 0.d0,                                       &
                    yvel = 0.d0
    END TYPE prtcl
    !==================================================================!
    !                        E N D    T Y P E :                        !
    !                            P R T C L                             !
    !==================================================================!
    
    ! Permissions of subroutines and functions
    PUBLIC  ::  getdist, mergeprt, prtcl
    
!----------------------------------------------------------------------!
CONTAINS

    !==================================================================!
    !                  B E G I N    F U N C T I O N :                  !
    !                          G E T D I S T                           !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Chris MacMackin                                !
    !   WRITTEN:        February 2014                                  !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        To find the distance between two particles.    !
    !                                                                  !
    !   ARGUMENTS:     *prtcl1  , first of the particles to find the   !
    !                             distance between.                    !
    !                  *prtcl2  , second of the particles to find the  !
    !                             distance between.                    !
    !   RETURNS:        Distance between the two particles passed as   !
    !                   arguments.                                     !
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
    REAL(8) FUNCTION getdist  ( prtcl1  , prtcl2   )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(prtcl), INTENT(IN) ::  prtcl1,                            &
                                    prtcl2
    !------------------------------------------------------------------!
    
        ! FIXME: Potential for round-off errors here
        getdist = SQRT((prtcl1%xpos - prtcl2%xpos)**2.d0 +             &
                       (prtcl1%ypos - prtcl2%ypos)**2.d0)
        
        IF ( ( getdist /= 0 ) .AND. ( getdist/                         &
           MAX(SQRT(prtcl1%xpos**2.d0 + prtcl1%ypos**2.d0),            &
           SQRT(prtcl2%xpos**2.d0 + prtcl2%ypos**2.d0)) < 1.d-12 ) )   &
           THEN
            WRITE( 6,1000)
        END IF
    
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('GETDIST : Warning: Possible loss of precision ',  &
                    'while calculation particle separation.')
    !------------------------------------------------------------------!
        
        RETURN
    END FUNCTION getdist
    !==================================================================!
    !                   E N D     F U N C T I O N :                    !
    !                          G E T D I S T                           !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         M E R G E P R T                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Merges the two particle passed as arguments.   !
    !                   Takes the sum of the two masses, set position  !
    !                   to centre of mass, and find velocity using     !
    !                   conservation of momentum. Sets one particle to !
    !                   null. If one particle is a star then merged    !
    !                   particle will also be a star.                  !
    !                                                                  !
    !   ARGUMENTS:     *prtcl1  , first of the particles to be merged. !
    !                  *prtcl2  , second of the particles to be        !
    !                             merged.                              !
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
    SUBROUTINE mergeprt ( prtcl1  , prtcl2   )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(prtcl), INTENT(INOUT)  ::  prtcl1,                        &
                                        prtcl2
        
        ! Other variables:
        INTEGER ::  newtype
        REAL(8) ::  newlum,                                            &
                    newmass,                                           &
                    newxpos,                                           &
                    newxvel,                                           &
                    newypos,                                           &
                    newyvel
    !------------------------------------------------------------------!
        
        ! If one or more particles is null type then do not merge
        IF ( ( prtcl1%ptype < 1 ) .OR. ( prtcl2%ptype < 1 ) ) THEN
            WRITE( 6,1000)
            RETURN
        END IF
        
        ! Get the new values for the merged particles:
        ! Get the type of the merged particle
        newtype = MAX(prtcl1%ptype, prtcl2%ptype)
        
        ! Set the new luminosity to the sum of the two old luminosities.
        ! This is highly contrived and nonphysical, and should be 
        ! corrected in future. However, for the simulations that this
        ! code is currently used for there will only be one luminous 
        ! particle and any dust particles falling into it will have a 
        ! tiny mass by comparison. Thus we can expect the luminosity to
        ! stay constant (the result which this subroutine will produce).
        newlum = prtcl1%lumin + prtcl2%lumin
        
        ! Combine the masses of the particles
        newmass = prtcl1%mass + prtcl2%mass
        
        ! Set position to centre of mass
        newxpos = (prtcl1%mass*prtcl1%xpos + prtcl2%mass*prtcl2%xpos)/ &
                  newmass
        newypos = (prtcl1%mass*prtcl1%ypos + prtcl2%mass*prtcl2%ypos)/ &
                  newmass
        
        ! Set the velocity of the new particle using conservation of 
        ! momentum
        newxvel = (prtcl1%mass*prtcl1%xvel + prtcl2%mass*prtcl2%xvel)/ &
                  newmass
        newyvel = (prtcl1%mass*prtcl1%yvel + prtcl2%mass*prtcl2%yvel)/ &
                  newmass
        
        ! Set the particles to the new values
        prtcl1%ptype = newtype
        prtcl1%lumin = newlum
        prtcl1%mass = newmass
        prtcl1%xpos = newxpos
        prtcl1%xvel = newxvel
        prtcl1%ypos = newypos
        prtcl1%yvel = newyvel
        
        prtcl2%ptype = 0
        
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('MERGEPRT: Warning: One or more particle is null ',&
                    'or fictitious.',/,                                &
                    'MERGEPRT: No merger performed.')
    !------------------------------------------------------------------!
    
        RETURN
    END SUBROUTINE mergeprt
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                         M E R G E P R T                          !
    !==================================================================!



    !==================================================================!
    !                  B E G I N    F U N C T I O N :                  !
    !                          T O M E R G E                           !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Chris MacMackin                                !
    !   WRITTEN:        February 2014                                  !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        To check whether two paticles should be merged !
    !                   before the next time-step.                     !
    !                                                                  !
    !   ARGUMENTS:     *prtcl1  , first of the particles being         !
    !                             checked.                             !
    !                  *prtcl2  , second of the particles being        !
    !                             checked.                             !
    !   RETURNS:        A logical that is true if the particles will   !
    !                   merge, false otherwise.                        !
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
    LOGICAL FUNCTION tomerge  ( prtcl1  , prtcl2  , sunpart  )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(prtcl), INTENT(IN) ::  prtcl1,                            &
                                    prtcl2,                            &
                                    sunpart
        
        ! Other variables:
        REAL(8) ::  hillrad,                                           &
                    rad1,                                              &
                    rad2,                                              &
                    sep,                                               &
                    vel
    !------------------------------------------------------------------!

        tomerge = .FALSE.
        sep = getdist(prtcl1,prtcl2)
        rad1 = ((3.d0*prtcl1%mass)/(4*pi*prtcl1%density))**(1.d0/3.d0)
        rad2 = ((3.d0*prtcl2%mass)/(4*pi*prtcl2%density))**(1.d0/3.d0)
        hillrad = (prtcl2%mass/sunpart%mass)**(1.d0/3.d0)*             &
                  getdist(prtcl2,sunpart)
        vel = SQRT((prtcl1%xvel - prtcl2%xvel)**2.d0 + (prtcl1%xvel -  &
              prtcl2%xvel)**2.d0)
        
        ! Check if particles in physical contact
        IF ( sep < (rad1 + rad2) ) THEN
            tomerge = .TRUE.
            RETURN
        END IF

        ! Check if particles within each other's Hill radius and
        ! travelling slowly enough to merge
        IF ( ( sep < hillrad ) .AND. ( vel < 2.d0*gconst*prtcl2%mass/  &
           ABS(sep**2.d0 - rad2**2.d0) ) ) THEN
            tomerge = .TRUE.
            RETURN
        END IF
        
        RETURN
    END FUNCTION tomerge
    !==================================================================!
    !                   E N D     F U N C T I O N :                    !
    !                          T O M E R G E                           !
    !==================================================================!


END MODULE particle
!======================================================================!
!                        E N D    M O D U L E :                        !
!                           U N T I T L E D                            !
!======================================================================!

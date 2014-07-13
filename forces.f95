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
!                             F O R C E S                              !
!======================================================================!
!                                                                      !
!   PURPOSE:    Contains subroutines for calculating the acceleration  !
!               of particles due to the various forces acting on them. !
!   CONTAINS:   gravity (subroutine), poyntrob (subroutine), radpress  !
!               (subroutine)                                           !
!   EXTERNALS:  particle (module), physvals (module)                   !
!                                                                      !
!----------------------------------------------------------------------!
MODULE forces
    USE physvals
    USE particle
    IMPLICIT NONE

CONTAINS

    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                          G R A V I T Y                           !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Chris MacMackin                                !
    !   WRITTEN:        February 2014                                  !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        To find the acceleration of prtcl1 due to the  !
    !                   gravitational force of prtcl2. Newtonian       !
    !                   gravity only.                                  !
    !                                                                  !
    !   ARGUMENTS:     *prtcl1  , the particle whose acceleration is   !
    !                             to be found.                         !
    !                  *prtcl2  , the particle whose gravity is        !
    !                             causing the acceleration.            !
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
    SUBROUTINE gravity  ( prtcl1  , prtcl2   )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(prtcl), INTENT(IN)     ::  prtcl2
        TYPE(prtcl), INTENT(INOUT)  ::  prtcl1
        
        ! Other variables:
        REAL(8) ::  dist,                                              &
                    ttlgrav
    !------------------------------------------------------------------!
        IF ( prtcl1%ptype < 1 ) THEN
            WRITE( 6,1000)
            RETURN
        ELSE IF ( prtcl2%ptype == 0 ) THEN
            WRITE( 6,1010)
            RETURN
        END IF
        
        ! TODO: Could subtraction done here be causing errors?
        ! (catastrophic cancellation?)
        dist = getdist(prtcl1,prtcl2)
        ttlgrav = gconst*prtcl2%mass/dist**2.d0
        prtcl1%xaccel = prtcl1%xaccel + ttlgrav * (prtcl2%xpos -       &
                        prtcl1%xpos)/dist
        prtcl1%yaccel = prtcl1%yaccel + ttlgrav * (prtcl2%ypos -       &
                        prtcl1%ypos)/dist
    
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('GRAVITY : Warning: Attempting to calculate ',     &
                    'force on null or fictitious',/,                   &
                    'GRAVITY : particle. Aborted.')
        1010 FORMAT('GRAVITY : Warning: Attempting to calculate ',     &
                    'gravity from null particle.',/,                   &
                    'GRAVITY : Aborted.')
    !------------------------------------------------------------------!
    
        RETURN
    END SUBROUTINE gravity
    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                          G R A V I T Y                           !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         P O Y N T R O B                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Chris MacMackin                                !
    !   WRITTEN:        February 2014                                  !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        To find the acceleration of prtcl1 due to the  !
    !                   Poynting-Robertson effect from prtcl2, a star  !
    !                   particle. Assumes the particle absorbs all     !
    !                   radiation (is a perfect black body) and both   !
    !                   particles are spheres.                         !
    !                                                                  !
    !   ARGUMENTS:     *prtcl1  , the particle whose acceleration is   !
    !                             to be found.                         !
    !                  *prtcl2  , the star particle whose radiation    !
    !                             is causing the acceleration.         !
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
    SUBROUTINE poyntrob ( prtcl1  , prtcl2   )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(prtcl), INTENT(IN)     ::  prtcl2
        TYPE(prtcl), INTENT(INOUT)  ::  prtcl1
        
        ! Other variables:
        REAL(8) ::  dist,                                              &
                    partrad,                                           &
                    starrad,                                           &
                    radforce,                                          &
                    sizedist,                                          &
                    tanforce,                                          &
                    radvel,                                            &
                    tanvel
    !------------------------------------------------------------------!
    
        IF ( ( prtcl1%ptype < 1 ) .OR. ( prtcl2%ptype /= 0 ) ) RETURN
        
        ! Get radial and tengential components of Poynting-Robertson
        ! effect, then break into x and y components
        dist = getdist(prtcl1,prtcl2)
        radvel = ((prtcl2%xvel - prtcl1%xvel)*(prtcl2%xpos -           &
                 prtcl1%xpos) + (prtcl2%yvel - prtcl1%yvel)*           &
                 (prtcl2%ypos - prtcl1%ypos))/dist
        tanvel = ((prtcl2%xvel - prtcl1%xvel)*(prtcl2%ypos -           &
                 prtcl1%ypos) + (prtcl2%yvel - prtcl1%yvel)*           &
                 (prtcl2%xpos - prtcl1%xpos))/dist
        partrad = (3.d0*prtcl1%mass/(4.d0*prtcl1%density*pi))**(1.d0/  &
                  3.d0)
        starrad = (3.d0*prtcl2%mass/(4.d0*prtcl2%density*pi))**(1.d0/  &
                  3.d0)
        sizedist = starrad/dist
        tanforce = -partrad**2.d0*prtcl2%lumin/(4.d0*cvel*starrad**2)* &
                   tanvel/cvel*(8.d0/3.d0 - 3.d0*(1.d0 -               &
                   sizedist**2.d0)**(1.d0/2.d0) + (1 -                 &
                   sizedist**2.d0)**(3.d0/2.d0)/3.d0)
        radforce = partrad**2.d0*prtcl2%lumin/(4.d0*cvel*starrad**2)*  &
                   (sizedist**2.d0 - radvel/cvel*(8.d0/3.d0 - 2*(1.d0  &
                   - sizedist**2.d0)**(1.d0/2.d0) - 2.d0/3.d0*(1 -     &
                   sizedist**2.d0)**(3.d0/2.d0)))
        
        prtcl1%xaccel = prtcl1%xaccel + (radforce*(prtcl2%xpos -       &
                        prtcl1%xpos) - tanforce*(prtcl2%ypos -         &
                        prtcl1%ypos))/(dist*prtcl1%mass)
        prtcl1%yaccel = prtcl1%yaccel + (radforce*(prtcl2%ypos -       &
                        prtcl1%ypos) + tanforce*(prtcl2%xpos -         &
                        prtcl1%xpos))/(dist*prtcl1%mass)
        
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('POYNTROB: Warning: Attempting to calculate ',     &
                    'force on null or fictitious',/,                   &
                    'POYNTROB: particle. Aborted.')
        1010 FORMAT('POYNTROB: Warning: Attempting to calculate ',     &
                    'Poynting-Robertson force',/,                      &
                    'POYNTROB: from non-radiating particle. Aborted.')
    !------------------------------------------------------------------!
    
        RETURN
    END SUBROUTINE poyntrob
    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         P O Y N T R O B                          !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         R A D P R E S S                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Chris MacMackin                                !
    !   WRITTEN:        February 2014                                  !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        To find the acceleration of prtcl1 due to the  !
    !                   radiation pressure from prtcl2, a star         !
    !                   particle. Assumes the particle absorbs all     !
    !                   radiation (is a perfect black body) and both   !
    !                   particles are spheres.                         !
    !                                                                  !
    !   ARGUMENTS:     *prtcl1  , the particle whose acceleration is   !
    !                             to be found.                         !
    !                  *prtcl2  , the star particle whose radiation    !
    !                             pressure is causing the              !
    !                             acceleration.                        !
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
    SUBROUTINE radpress ( prtcl1  , prtcl2   )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(prtcl), INTENT(IN)     ::  prtcl2
        TYPE(prtcl), INTENT(INOUT)  ::  prtcl1
        
        ! Other variables:
        REAL(8) ::  dist       ,                                       &
                    ttlrad
    !------------------------------------------------------------------!
    
        IF ( ( prtcl1%ptype < 1 ) .OR. ( prtcl2%ptype /= 0 ) ) RETURN
        
        ! Get the total force from radiation pressure then divide it
        ! into x and y components
        dist = getdist(prtcl1,prtcl2)
        ttlrad = prtcl2%lumin*(3.d0/(4.d0*pi*prtcl1%mass*              &
                 prtcl1%density))**(2.d0/3.d0)/(4.d0*dist**2.d0)
        prtcl1%xaccel = prtcl1%xaccel + ttlrad * (prtcl2%xpos -        &
                        prtcl1%xpos)/dist
        prtcl1%yaccel = prtcl1%yaccel + ttlrad * (prtcl2%ypos -        &
                        prtcl1%ypos)/dist
    
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('RADPRESS: Warning: Attempting to calculate ',     &
                    'force on null or fictitious',/,                   &
                    'RADPRESS: particle. Aborted.')
        1010 FORMAT('RADPRESS: Warning: Attempting to calculate ',     &
                    'radiation pressure from',/,                       &
                    'RADPRESS: non-radiating particle. Aborted.')
    !------------------------------------------------------------------!
    
        RETURN
    END SUBROUTINE radpress
    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         R A D P R E S S                          !
    !==================================================================!
    
END MODULE forces
!======================================================================!
!                        E N D    M O D U L E :                        !
!                             F O R C E S                              !
!======================================================================!

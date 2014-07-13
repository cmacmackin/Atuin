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
!                            R A N D G E N                             !
!======================================================================!
!                                                                      !
!   PURPOSE:    Contains functions for generating (pseudo)random       !
!               numbers.                                               !
!   CONTAINS:   normrand (function), unirand (function), unirand       !
!               (function)                                             !
!   EXTERNALS:  None                                                   !
!                                                                      !
!----------------------------------------------------------------------!
MODULE randgen
    IMPLICIT NONE
    
    ! Variable declarations:
    LOGICAL, PRIVATE    ::  initd = .FALSE.

CONTAINS

    !==================================================================!
    !                  B E G I N    F U N C T I O N :                  !
    !                         N O R M R A N D                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         John Burkardt                                  !
    !   WRITTEN:        Unknown                                        !
    !   MODIFICATIONS:  August, 2013                                   !
    !                   March, 2014                                    !
    !                       Adapted so that only need to pass an       !
    !                       argument when initializing function, by    !
    !                       Chris MacMackin.                           !
    !                                                                  !
    !   PURPOSE:        This routine returns a normal random deviate   !
    !                   with mean 0.0 and standard deviation of 1.0.   !
    !                   Modified version of code available from        !
    !                   http://people.sc.fsu.edu/~jburkardt/f_src/     !
    !                   uniform/uniform.html, distributed under the    !
    !                   GNU LGPL.                                      !
    !                                                                  !
    !   ARGUMENTS:      None                                           !
    !   RETURNS:        A normally distributed random double precision !
    !                   real value.                                    !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    REAL(8) FUNCTION normrand ( init     )
        IMPLICIT NONE

        ! Input and output variables:
        INTEGER(4), INTENT(IN), OPTIONAL    ::  init  

        ! Other variables:
        INTEGER             ::  seed
        REAL(8)             ::  r1,                                    &
                                r2,                                    &
                                x
        REAL(8), PARAMETER  ::  pi = 3.141592653589793d0
    !------------------------------------------------------------------!

        IF ( PRESENT(init) .EQV. .TRUE. ) THEN
            r1 = unirand (init)
        ELSE IF ( initd .EQV. .FALSE. ) THEN
            IF ( seed == 0 ) seed = 1
            r1 = unirand(seed)
        ELSE
            r1 = unirand()
        END IF
        r2 = unirand()
        
        x = SQRT(-2.d0*LOG(r1))*COS(2.d0*pi*r2)
        normrand = x

        RETURN
    END FUNCTION normrand
    !==================================================================!
    !                   E N D     F U N C T I O N :                    !
    !                         N O R M R A N D                          !
    !==================================================================!



    !==================================================================!
    !                  B E G I N    F U N C T I O N :                  !
    !                          U N I R A N D                           !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         John Burkardt                                  !
    !   WRITTEN:        Unknown                                        !
    !   MODIFICATIONS:  May, 2007                                      !
    !                   March, 2014                                    !
    !                       Adapted so that only need to pass an       !
    !                       argument when initializing function, by    !
    !                       Chris MacMackin.                           !
    !                                                                  !
    !   PURPOSE:        This routine returns a uniform random deviate  !
    !                   between 0.0 and 1.0. Modified version of code  !
    !                   available from http://people.sc.fsu.edu/       !
    !                   ~jburkardt/f_src/uniform/uniform.html,         !
    !                   distributed under the GNU LGPL.                !
    !                                                                  !
    !   ARGUMENTS:     *init    , an integer which provides initial    !
    !                             seed for first call of this          ! 
    !                             function. Seed must not be equal to  !
    !                             zero. Passing this argument will     !
    !                             reset the random number sequence.    !
    !                             (OPTIONAL)                           !
    !   RETURNS:        A new pseudorandom variate, strictly between 0 !
    !                   and 1.                                         !
    !                                                                  !
    !   EXTERNALS:      None                                           !
    !                                                                  !
    !------------------------------------------------------------------!
    REAL(8) FUNCTION unirand  ( init )
        IMPLICIT NONE

        ! Input and output variables:
        INTEGER(4), INTENT(IN), OPTIONAL    ::  init  

        ! Other variables:
        INTEGER(4)              ::  k
        INTEGER(4), PARAMETER   ::  i4_huge = 2147483647
        INTEGER(4), SAVE        ::  seed = 0
    !------------------------------------------------------------------!
    
        IF ( PRESENT(init) .EQV. .TRUE. ) THEN
            seed = init
            initd = .TRUE.
        END IF
        
        IF ( seed == 0 ) THEN
            WRITE( 0,1000)
            STOP
        END IF

        k = seed / 127773
        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        IF ( seed < 0 ) THEN
            seed = seed + i4_huge
        END IF

        unirand = REAL(seed) * 4.656612875d-10

    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('UNIRAND : ERROR: Value of seed equal to zero.')
    !------------------------------------------------------------------!

        RETURN
    END FUNCTION unirand 
    !==================================================================!
    !                   E N D     F U N C T I O N :                    !
    !                          U N I R A N D                           !
    !==================================================================!

END MODULE randgen    
!======================================================================!
!                        E N D    M O D U L E :                        !
!                            R A N D G E N                             !
!======================================================================!

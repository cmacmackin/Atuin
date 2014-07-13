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
!                           P H Y S V A L S                            !
!======================================================================!
!                                                                      !
!   PURPOSE:  Contains physical constants which can be easily imported !
!             into other programs. This is an evolving module, with    !
!             new constants being added as they are needed. All values !
!             are SI, unless otherwise indicated.                      !
!   CONTAINS: None                                                     !
!                                                                      !
!----------------------------------------------------------------------!
MODULE physvals
    IMPLICIT NONE

    REAL(8), PARAMETER  ::  auval = 1.4959787066d11,                   &
                            d_earth = 5.514d3,                         &
                            d_sol = 1.409d3,                           &
                            cvel = 2.99792458d8,                       &
                            gconst = 6.6742867d-11,                    &
                            hbar = 1.05457162853d-34,                  &
                            l_sol = 3.8395d26,                         &
                            m_earth = 5.9736d24,                       &
                            m_elec = 9.10938215d-31,                   &
                            m_sol = 1.9891d30,                         &
                            pi = 3.1415926535897932d0,                 &
                            q_elec = 1.60217648740d-19,                &
                            secsyear = 3.15576d7

END MODULE physvals    
!======================================================================!
!                        E N D    M O D U L E :                        !
!                           P H Y S V A L S                            !
!======================================================================!

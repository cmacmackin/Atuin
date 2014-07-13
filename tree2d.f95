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
!                             T R E E 2 D                              !
!======================================================================!
!                                                                      !
!   PURPOSE:    Provides data structure and subroutines for tree code  !
!               algorithms used in 2D n-body simulations. Currently    !
!               designed for systems with few stars (eg: not globular  !
!               clusters or galaxies).                                 !
!   CONTAINS:   allaccel (subroutine), nodeacc (subroutine), treenode  !
!               (derived type), updtall (subroutine), updtnode         !
!               (subroutine)                                           !
!   EXTERNALS:  pdstrbtn (module), particle (module)                   !
!                                                                      !
!----------------------------------------------------------------------!
MODULE tree2d
    USE particle
    USE pdstrbtn
    USE physvals
    IMPLICIT NONE

    !==================================================================!
    !                      B E G I N    T Y P E :                      !
    !                         T R E E N O D E                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        A node in the tree used to track particles.    !
    !                   May be a "leaf" containing 0 or 1 particles,   !
    !                   with no child nodes, or a "branch" containing  !
    !                   4 child nodes which may be branches or leaves  !
    !                   in their own right.                            !
    !                                                                  !
    !------------------------------------------------------------------!
    TYPE treenode
        INTEGER                                 ::  numprtcl = 0
        LOGICAL                                 ::  hasstar = .FALSE.
        REAL(8)                                 ::  xmin = 0.d0,       &
                                                    xmax = 0.d0,       &
                                                    xsize = 0.d0,      &
                                                    ymin = 0.d0,       &
                                                    ymax = 0.d0,       &
                                                    ysize = 0.d0
        TYPE(treenode), DIMENSION(:,:), POINTER ::  children => NULL()
        TYPE(prtcl), POINTER                    ::  contents => NULL()
    END TYPE treenode
    !==================================================================!
    !                        E N D    T Y P E :                        !
    !                         T R E E N O D E                          !
    !==================================================================!
    
    ! Variable declarations
    INTEGER, PUBLIC                 ::  subnumx = 2,                   &
                                        subnumy = 2
    REAL(8), PUBLIC                 ::  absdvel = 1.5d4,               &
                                        approxif = 1.d0,               &
                                        chngdtif = 5.d0,               &
                                        margin = 1.d-1,                &
                                        reducedt = 5.d-1,              &
                                        velchng = 1.d-1
    TYPE(treenode), SAVE, PRIVATE   ::  topnode
    
    ! Permissions
    PUBLIC  ::  updtall, allaccel
    PRIVATE ::  treenode, nodeacc, rmbelow, updtnode
    
!----------------------------------------------------------------------!
CONTAINS

    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                          A L L A C C E L                         !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Iterates through all non-null particles and    !
    !                   computes acceleration experienced by each, as  !
    !                   well as checking for mergers. This subroutine  !
    !                   is the public interface for the nodeacc        !
    !                   subroutine.                                    !
    !                                                                  !
    !   ARGUMENTS:     *getaccel, the function which computes the      !
    !                             value of the acceleration between    !
    !                             two particles.                       !
    !                                                                  !
    !                   Arguments marked with an asterisk (*) are ones !
    !                   for which the user must specify a value when   !
    !                   the subroutine is called. Those without an     !
    !                   asterisk (*) are used only to return a value   !
    !                   to the caller.                                 !
    !                                                                  !
    !   EXTERNALS:      getaccel (subroutine)                          !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE allaccel ( getaccel )     
        IMPLICIT NONE
    
        ! Give compiler information about function used to compute value
        ! of acceleration between particles
        INTERFACE
            SUBROUTINE getaccel ( prtcl1  , prtcl2   )
                USE particle
                TYPE(prtcl), INTENT(IN)     ::  prtcl2
                TYPE(prtcl), INTENT(INOUT)  ::  prtcl1
            END SUBROUTINE getaccel
        END INTERFACE
        
        INTEGER                 ::  iter
        REAL(8)                 ::  accel,                             &
                                    dt,                                &
                                    mindt,                             &
                                    vel
        TYPE(prtcl), POINTER    ::  mrgwith
    !------------------------------------------------------------------!
        
        mindt = HUGE(mindt)
        numnnull = 0
        energy = 0.d0
        angmom = 0.d0
        xlinmom = 0.d0
        ylinmom = 0.d0
        xcom = 0.d0
        ycom = 0.d0
        totmass = 0.d0
        
        ! Iterate through non-null particles and compute acceleration.
        ! Also calculate what time-step to use for next iteration and
        ! check if particles should be merged.
        DO iter = 1, ttlprtcl
            IF ( partlist(iter)%ptype > 0 ) THEN
                IF ( ASSOCIATED(mrgwith) .EQV. .TRUE. ) NULLIFY(mrgwith)
                numnnull = numnnull + 1
                partlist(iter)%xaccel = 0.d0
                partlist(iter)%yaccel = 0.d0
                
                ! Merge particle as many times as necessary
                DO 
                    CALL nodeacc(topnode,partlist(iter),getaccel,      &
                                 mrgwith)
                    IF ( ASSOCIATED(mrgwith) .EQV. .TRUE. ) THEN
                        CALL mergeprt(partlist(iter),mrgwith)
                        NULLIFY(mrgwith)
                        CALL updtall
                    ELSE
                        EXIT
                    END IF
                END DO
        
                ! TODO: Might not have quite fixed v=0 case properly
                ! Time-step update
                accel = SQRT(partlist(iter)%xaccel**2.d0 +             &
                        partlist(iter)%yaccel**2.d0)
                vel = SQRT(partlist(iter)%xvel**2.d0 +                 &
                      partlist(iter)%yvel**2.d0)
                IF ( accel == 0.d0 ) THEN
                    dt = HUGE(dt)
                ELSE IF ( vel == 0.d0 ) THEN
                    dt = absdvel/accel
                ELSE
                    dt = vel*velchng/accel
                END IF
                mindt = MIN(mindt,dt)
                
                ! Update energy, angular/liner momentum and centre of 
                ! mass trackers
                angmom = angmom + partlist(iter)%mass*                 &
                         (partlist(iter)%xpos*partlist(iter)%yvel -    &
                         partlist(iter)%ypos*partlist(iter)%xvel)
                energy = energy + 5.d-1*partlist(iter)%mass*           &
                         (partlist(iter)%xvel**2.d0 +                  &
                         partlist(iter)%yvel**2.d0)
                xlinmom = xlinmom + partlist(iter)%mass*               &
                          partlist(iter)%xvel
                ylinmom = ylinmom + partlist(iter)%mass*               &
                          partlist(iter)%yvel
                xcom = xcom + partlist(iter)%mass*partlist(iter)%xpos
                ycom = ycom + partlist(iter)%mass*partlist(iter)%ypos
                totmass = totmass + partlist(iter)%mass
            END IF
        END DO
        
        ! Update time-step, if necessary
        IF ( mindt < deltat ) THEN
            deltat = reducedt*mindt
            nchngdt = nchngdt + 1
            WRITE( 6,1000) simtime/secsyear
        ELSE IF ( mindt > chngdtif*deltat ) THEN
            deltat = reducedt*mindt
            nchngdt = nchngdt + 1
            WRITE( 6,1010) simtime/secsyear
        END IF

        xcom = xcom/totmass
        ycom = ycom/totmass
        
    !------------------------------------------------------------------!
    !                      Write format statements                     !
    !------------------------------------------------------------------!
        1000 FORMAT('ALLACCEL: Decreased time-step size at simulation',&
                    ' time ',1PG12.5,' years.')
        1010 FORMAT('ALLACCEL: Increased time-step size at simulation',&
                    ' time ',1PG12.5,' years.')
    !------------------------------------------------------------------!
    
        RETURN
    END SUBROUTINE allaccel
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                          A L L A C C E L                         !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                          N O D E A C C                           !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Recursively finds the contribution of the      !
    !                   particles contained within the current node to !
    !                   the acceleration experienced by the particle   !
    !                   in the argument list. Also checks if the       !
    !                   particle needs to merge with another.          !
    !                                                                  !
    !   ARGUMENTS:     *curnode , the node whose contribution to the   !
    !                             acceleration is being calculated.    !
    !                  *curpart , the particle whose acceleration is   !
    !                             being calculated.                    !
    !                  *getaccel, the function which computes the      !
    !                             value of the acceleration between    !
    !                             two particles.                       !
    !                   mrgwith , a prtcl type pointer which will be   !
    !                             allocated if curpart needs to merge. !
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
    RECURSIVE SUBROUTINE nodeacc  ( curnode , curpart , getaccel,      &
                                    mrgwith  )
        IMPLICIT NONE
        
        ! Give compiler information about function used to compute value
        ! of acceleration between particles
        INTERFACE
            SUBROUTINE getaccel ( prtcl1  , prtcl2   )
                USE particle
                TYPE(prtcl), INTENT(IN)     ::  prtcl2
                TYPE(prtcl), INTENT(INOUT)  ::  prtcl1
            END SUBROUTINE getaccel
        END INTERFACE
        
        ! Input and output variables:
        TYPE(treenode), INTENT(IN)          ::  curnode
        TYPE(prtcl), TARGET, INTENT(INOUT)  ::  curpart
        TYPE(prtcl), POINTER                ::  mrgwith
        
        ! Other variables
        INTEGER ::  iter
    !------------------------------------------------------------------!
    
        ! If current node is empty, do nothing
        IF ( ASSOCIATED(curnode%contents) .EQV. .FALSE. ) THEN
            RETURN
        ! If current node is a leaf and doesn't contain the same 
        ! particle as the one for which the acceleration is being  
        ! found, check whether particles will merge and, if not, 
        ! calculate the acceleration
        ELSE IF ( curnode%numprtcl == 1 ) THEN
            IF ( ASSOCIATED(curnode%contents,curpart) .EQV. .FALSE. )  &
               THEN
                IF ( tomerge(curpart,curnode%contents,partlist(1))     &
                  .EQV. .TRUE. ) THEN
                    !TODO: rather than having to do a complete update of
                    ! the tree, would I be able to just do one here?
                    mrgwith => curnode%contents
                ELSE
                    CALL getaccel(curpart,curnode%contents)
                    energy = energy - 5.d-1*(gconst*curpart%mass*      &
                             curnode%contents%mass/getdist(curpart,    &
                             curnode%contents))
                END IF
            END IF
            RETURN
        ! If fictitious particle in the current node is sufficiently far
        ! away then compute accelerations using fictitious particle 
        ! (if doesn't contain a star)
        ! FIXME: Might need to change how this works for rectangular 
        ! nodes
        ELSE IF (  MAX(curnode%xsize,curnode%ysize)/                   &
                 getdist(curnode%contents,curpart) < approxif ) THEN
            CALL getaccel(curpart,curnode%contents)
            energy = energy - 5.d-1*(gconst*curpart%mass*              &
                     curnode%contents%mass/getdist(curpart,            &
                     curnode%contents))
            RETURN
        ! Otherwise repeat this process with child nodes
        ELSE
            DO iter = 0, (subnumx*subnumy - 1)
                CALL nodeacc(curnode%children(iter/subnumx + 1,        &
                             MOD(iter,subnumy) + 1),curpart,getaccel,  &
                             mrgwith)
                IF ( ASSOCIATED(mrgwith) .EQV. .TRUE. ) EXIT
            END DO
            RETURN
        END IF
        
    END SUBROUTINE nodeacc
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                          N O D E A C C                           !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                           R M B E L O W                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Recursively removes all children of the node   !
    !                   in the argument list.                          !
    !                                                                  !
    !   ARGUMENTS:     *curnode , the node whose children (if they     !
    !                   exist) are to be removed.                      !
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
    RECURSIVE SUBROUTINE rmbelow ( curnode  )
        IMPLICIT NONE
        
        ! Input and output variables:
        TYPE(treenode), INTENT(INOUT)   ::  curnode
        
        ! Other variables:
        INTEGER ::  iter
    !------------------------------------------------------------------!
        
        IF ( ASSOCIATED(curnode%children) ) THEN
            ! Remove any children of the child nodes
            DO iter = 0, (subnumx*subnumy - 1)
                CALL rmbelow(curnode%children(iter/subnumx + 1,        &
                             MOD(iter,subnumy) + 1))
            END DO
            ! Remove the child nodes
            DEALLOCATE(curnode%children)
            DEALLOCATE(curnode%contents)
        END IF
        
        IF ( ASSOCIATED(curnode%contents) .EQV. .TRUE. )               &
            NULLIFY(curnode%contents)
        
    END SUBROUTINE rmbelow
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                           R M B E L O W                          !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                           U P D T A L L                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Dynamically updates the entire tree structure  !
    !                   so that there is only one particle per "leaf". !
    !                   This subroutine is the public interface for    !
    !                   the updtnode subroutine (which would be far    !
    !                   easier to call incorrectly).                   !
    !                                                                  !
    !   ARGUMENTS:      None                                           !
    !                                                                  !
    !   EXTERNALS:      updtnode (subroutine)                          !
    !                                                                  !
    !------------------------------------------------------------------!
    SUBROUTINE updtall        
        IMPLICIT NONE
        
        ! Declare and initialize variables
        INTEGER                         ::  iter
        INTEGER, DIMENSION(ttlprtcl,2)  ::  cnddts
        LOGICAL                         ::  resizex,                   &
                                            resizey
    !------------------------------------------------------------------!

        ! Initialize and update variables
        FORALL(iter = 1:ttlprtcl) cnddts(iter,1) = iter
        cnddts(:,2) = 0
        topnode%numprtcl = numnnull
        
        ! Prevent divide by zero errors
        IF ( xlargest == xsmallst ) xlargest = xlargest + TINY(xlargest)
        IF ( ylargest == ysmallst ) ylargest = ylargest + 1.d0
        
        ! Deactivate old tree
        CALL rmbelow(topnode)
        IF ( ASSOCIATED(topnode%contents) ) DEALLOCATE(topnode%contents)
        
        ! Get size of tree
        topnode%xmax = xlargest + margin*(xlargest - xsmallst)
        topnode%xmin = xsmallst - margin*(xlargest - xsmallst)
        topnode%xsize = topnode%xmax - topnode%xmin
        topnode%ymax = ylargest + margin*(ylargest - ysmallst)
        topnode%ymin = ysmallst - margin*(ylargest - ysmallst)
        topnode%ysize = topnode%ymax - topnode%ymin
        resizex = .TRUE.
        resizey = .TRUE.

        ! Update tree
        CALL updtnode(cnddts,topnode,numnnull,resizex,resizey)

        RETURN
    END SUBROUTINE updtall
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                           U P D T A L L                          !
    !==================================================================!



    !==================================================================!
    !                B E G I N    S U B R O U T I N E :                !
    !                         U P D T N O D E                          !
    !==================================================================!
    !                                                                  !
    !   AUTHOR:         Christopher MacMackin                          !
    !   WRITTEN:        February, 2014                                 !
    !   MODIFICATIONS:  None                                           !
    !                                                                  !
    !   PURPOSE:        Dynamically adapts the tree structure so that  !
    !                   there is only one particle per "leaf". Each    !
    !                   branch node contains a fictitious particle     !
    !                   which is located at the centre of mass of all  !
    !                   of the real particles contained in the child   !
    !                   nodes.                                         !
    !                                                                  !
    !   ARGUMENTS:     *cndddts , the particles contained within the   !
    !                             parent node. This is a 2D array,     !
    !                             with the array index of the particle !
    !                             in the first column and the second   !
    !                             column contains a number, set to 0   !
    !                             if the particle has not been         !
    !                             assigned to another (sister) node    !
    !                             already and 1 if it has been         !
    !                             assigned to a node.                  !
    !                  *curnode , the node for which to find the mass  !
    !                             and location of the particle or      !
    !                             fictitious particle.                 !
    !                  *numabove, the number of particles positioned   !
    !                             within the parent node.              !
    !                  *resizex , indicates whether the size of this   !
    !                             node has been adjusted in the x-     !
    !                             dimension.                           !
    !                  *resizey , indicates whether the size of this   !
    !                             node has been adjusted in the y-     !
    !                             dimension.                           !
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
    RECURSIVE SUBROUTINE updtnode ( cnddts  , curnode , numabove,      &
                                    resizex , resizey  )
        IMPLICIT NONE
        
        ! Input and output variables:
        INTEGER, INTENT(IN)                     ::  numabove
        LOGICAL, INTENT(IN)                     ::  resizex,           &
                                                    resizey
        INTEGER, INTENT(INOUT), DIMENSION(:,:)  ::  cnddts
        TYPE(treenode), INTENT(INOUT)           ::  curnode
        
        ! Other variables:
        INTEGER, DIMENSION(numabove,2)  ::  innode
        INTEGER                         ::  iter,                      &
                                            num
        LOGICAL                         ::  newkids
        TYPE(treenode), POINTER         ::  child
    !------------------------------------------------------------------!
        ! Initialize vtestariables and pointers
        curnode%hasstar = .FALSE.
        curnode%numprtcl = 0
        innode(:,2) = 0
        iter = 0
        num = 0
        newkids = .FALSE.
        
        ! Find which particles are contained in current node
        DO iter = 1, numabove
            num = cnddts(iter,1)
            IF ( ( cnddts(iter,2) == 0 ) .AND. ( partlist(num)%xpos >= &
                 curnode%xmin ) .AND. ( partlist(num)%xpos <           &
                 curnode%xmax ) .AND. ( partlist(num)%ypos >=          &
                 curnode%ymin ) .AND. ( partlist(num)%ypos <           &
                 curnode%ymax ) .AND. ( partlist(num)%ptype > 0 ) ) THEN
                curnode%numprtcl = curnode%numprtcl + 1
                innode(curnode%numprtcl,1) = num
                cnddts(iter,2) = 1
                IF ( partlist(num)%ptype == 2 ) THEN
                    curnode%hasstar = .TRUE.
                END IF
            END IF
        END DO
        
        ! If leaf, assign particle to it (if one is present)
        IF ( curnode%numprtcl == 0 ) THEN
            ! If node used to hold a particle, remove it
            IF ( ASSOCIATED(curnode%contents) .EQV. .TRUE. ) THEN
                NULLIFY(curnode%contents)
            END IF
            ! If node has child nodes, remove them and the node's 
            ! fictitious particle
            IF ( ASSOCIATED(curnode%children) .EQV. .TRUE. ) THEN
                CALL rmbelow(curnode)
                NULLIFY(curnode%contents)
            END IF
            RETURN
        ELSE IF ( curnode%numprtcl == 1 ) THEN
            ! If node has child nodes, remove them and node's fictitious
            ! particle
            IF ( ASSOCIATED(curnode%children) .EQV. .TRUE. ) THEN
                CALL rmbelow(curnode)
                DEALLOCATE(curnode%contents)
            END IF
            curnode%contents => partlist(innode(1,1))
            RETURN
        END IF
        
        ! If not already child nodes, create them and a fictitious 
        ! particle. Otherwise, reset the fictitious particle
        IF ( ASSOCIATED(curnode%children) .EQV. .FALSE. ) THEN
            ALLOCATE(curnode%children(subnumx,subnumy),curnode%contents)
            newkids = .TRUE.
            curnode%contents%ptype = -1
        ELSE
            curnode%contents%mass = 0.d0
            curnode%contents%xpos = 0.d0
            curnode%contents%ypos = 0.d0
            curnode%contents%ptype = -1
        END IF
        
        ! Initialize child nodes, recursively sort particles into them,
        ! and use result to find fictitious particle for current node
        DO iter = 0, (subnumx*subnumy - 1)
            ! Associate pointer with child node
            child => curnode%children(iter/subnumy + 1,MOD(iter,       &
                     subnumx) + 1)
            
            ! Set boundaries and sizes of child nodes if necessary
            IF ( ( newkids .EQV. .TRUE. ) .OR. ( resizex .EQV.         &
                 .TRUE. ) ) THEN
                child%xsize = curnode%xsize/DBLE(subnumx)
                IF ( iter/subnumx > 0 ) THEN
                    child%xmin = curnode%children(iter/subnumx,        &
                                 MOD(iter,subnumy)+1)%xmax
                ELSE
                    child%xmin = curnode%xmin
                END IF
                child%xmax = child%xmin + child%xsize
            END IF
            IF ( ( newkids .EQV. .TRUE. ) .OR. ( resizey .EQV.         &
                 .TRUE. ) ) THEN
                child%ysize = curnode%ysize/DBLE(subnumy)
                IF ( MOD(iter,subnumy) > 0 ) THEN
                    child%ymin = curnode%children(iter/subnumx+1,      &
                                 MOD(iter,subnumy))%ymax
                ELSE
                    child%ymin = curnode%ymin
                END IF
                child%ymax = child%ymin + child%ysize
            END IF
            
            ! Recursively fill child node
            CALL updtnode(innode(1:curnode%numprtcl,:),child,          &
                          curnode%numprtcl,resizex,resizey)

            ! Get values for fictitious particle
            IF ( ASSOCIATED(child%contents) .EQV. .TRUE. ) THEN
                curnode%contents%mass = curnode%contents%mass +        &
                                        child%contents%mass
                curnode%contents%xpos = curnode%contents%xpos +        &
                                        child%contents%xpos*           &
                                        child%contents%mass
                curnode%contents%ypos = curnode%contents%ypos +        &
                                        child%contents%ypos*           &
                                        child%contents%mass
            END IF
        END DO
        curnode%contents%xpos = curnode%contents%xpos/                 &
                                curnode%contents%mass
        curnode%contents%ypos = curnode%contents%ypos/                 &
                                curnode%contents%mass
        NULLIFY(child)
            
        RETURN
    END SUBROUTINE updtnode
    !==================================================================!
    !                  E N D    S U B R O U T I N E :                  !
    !                         U P D T N O D E                          !
    !==================================================================!

END MODULE tree2d
!======================================================================!
!                        E N D    M O D U L E :                        !
!                             T R E E 2 D                              !
!======================================================================!

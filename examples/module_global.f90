!======================================================================================================================
!> @file        module_global.f90
!> @brief       This is for an example case of p3ta.
!> @author      
!>              - Kiha Kim (k-kiha@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>
!> @date        June 2019
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2019 Kiha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

!>
!> @brief       Module for global parameters
!>
module module_global
    implicit none
    double precision, parameter :: PI = acos(-1.d0)
    
    ! PHYSICAL PARAMETERS
    double precision :: Pr  !< Prandtl number
    double precision :: Ra  !< Reyleigh number

    ! ITERATION STEPS 
    integer :: Tmax             !< Maximum number of iteration steps
    
    ! COMPUTATIONAL SIZE FOR SPACE AND TIME DISCRETIZATIONS
    !> @{ Grid numbers in each direction
    integer :: nx,ny,nz 
    !> @}
    !> @{ Grid numbers minus 1 in each direction
    integer :: nxm,nym,nzm
    !> @}
    !> @{ Grid numbers plus 1 in each direction
    integer :: nxp,nyp,nzp
    !> @}
    double precision :: dt                  !< Length of time step
    double precision :: dtStart             !< Initial dt
    double precision :: tStart              !< Initial simulation time
    
    ! DOMAIN SIZE FOR THE PHYSICAL PROBLEMS
    !> @{ Lengths of the physical domain
    double precision :: lx, ly, lz
    !> @}
    !> @{ Discretized grid lengths of the physical domain
    double precision :: dx, dy, dz
    !> @}
    
    ! BOUNDARY CONDITIONS OF THE HOT AND COLD WALLS AND OTHER PROPERTIES
    double precision :: theta_cold              !< Boundary temperature of cold wall
    double precision :: theta_hot               !< Boundary temperature of hot wall
    double precision :: alphaG                  !< Thermal expansion coefficient x gravitational acceleration
    double precision :: nu                      !< Kinematic viscosity
    double precision :: Ct                      !< Thermal diffusivity

    contains 
    !>
    !> @brief       Module for global parameters
    !> @param       np_dim      Number of MPI processes in 3D topology
    !>
    subroutine module_global_inputpara(np_dim)
        implicit none
        integer, intent(out) :: np_dim(0:2)

        integer :: npx, npy, npz   ! Variables to read number of processes in 3D topology

        ! Namelist variables for file input
        namelist /meshes/ nx, ny, nz
        namelist /procs/ npx, npy, npz
        namelist /time/ tmax

        open (unit = 1, file = "PARA_INPUT.inp")
            read (1, meshes)
            read (1, procs)
            read (1, time)
        close (1)

        np_dim(0) = npx
        np_dim(1) = npy
        np_dim(2) = npz

        ! PHYSICAL PARAMETERS
        Pr = 5.0d0; Ra = 2.d+2

        ! COMPUTATIONAL SIZE FOR SPACE AND TIME DISCRETIZATIONS
        nx = nx+1; ny = ny+1; nz = nz+1
        nxm = nx-1; nym = ny-1; nzm = nz-1
        nxp = nx+1; nyp = ny+1; nzp = nz+1

        dtStart = 5.0D-3; tStart = 0.d0

        ! DOMAIN SIZE FOR THE PHYSICAL PROBLEMS
        lx = 1.0d0; ly = 1.0d0; lz = 1.0d0

        ! BOUNDARY CONDITIONS OF THE HOT AND COLD WALLS AND OTHER PROPERTIES
        theta_cold = -1.d0; theta_hot = 2.d0 + theta_cold
        alphaG = 1.d0; nu = 1.d0/sqrt(Ra/(alphaG*Pr*ly**3.*(theta_hot-theta_cold)))
        Ct = nu/Pr
        
    end subroutine module_global_inputpara

end module module_global
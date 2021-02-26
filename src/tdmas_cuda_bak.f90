!>
!> @brief       Solve many tridiagonal systems of equations using the Thomas algorithm.
!>              First index indicates the number of independent many tridiagonal systems to use vectorization.
!>              Second index indicates the row number in the tridiagonal system .
!> @param       a       Coefficient array in lower diagonal elements
!> @param       b       Coefficient array in diagonal elements
!> @param       c       Coefficient array in upper diagonal elements
!> @param       d       Coefficient array in the right-hand side terms
!> @param       n1      Number of rows in each process, size of the tridiagonal matrix N divided by nprocs
!> @param       n2      Number of tridiagonal systems per process
!>
subroutine tdma_many_cuda(a, b, c, d, n1, n2, nthds)

    use omp_lib

    implicit none

    integer, intent(in) :: n1,n2,nthds
    double precision, intent(inout) :: a(n1,n2,nthds), b(n1,n2,nthds), c(n1,n2,nthds), d(n1,n2,nthds)
    
    integer :: i,j,ti
    double precision, allocatable, dimension(:) :: r

    allocate(r(1:n2))

    do ti=1,nthds
        do j=1,n2
            d(1,j,ti)=d(1,j,ti)/b(1,j,ti)
            c(1,j,ti)=c(1,j,ti)/b(1,j,ti)
        enddo

        do j=1,n2
            do i=2,n1
                r(j)=1.d0/(b(i,j,ti)-a(i,j,ti)*c(i-1,j,ti))
                d(i,j,ti)=r(j)*(d(i,j,ti)-a(i,j,ti)*d(i-1,j,ti))
                c(i,j,ti)=r(j)*c(i,j,ti)
            enddo
        enddo

        do j=1,n2
            do i=n1-1,1,-1
                d(i,j,ti)=d(i,j,ti)-c(i,j,ti)*d(i+1,j,ti)
            enddo
        enddo
    enddo

    deallocate(r)

end subroutine tdma_many_cuda

!>
!> @brief       Solve many cyclic tridiagonal systems of equations using the Thomas algorithm.
!>              First index indicates the number of independent many tridiagonal systems to use vectorization.
!>              Second index indicates the row number in the tridiagonal system.
!> @param       a       Coefficient array in lower diagonal elements
!> @param       b       Coefficient array in diagonal elements
!> @param       c       Coefficient array in upper diagonal elements
!> @param       d       Coefficient array in the right-hand side terms
!> @param       n1      Number of rows in each process, size of the tridiagonal matrix N divided by nprocs
!> @param       n2      Number of tridiagonal systems per process
!>
subroutine tdma_cycl_many_cuda(a, b, c, d, n1, n2, nthds)

    use omp_lib

    implicit none

    integer, intent(in) :: n1,n2,nthds
    double precision, intent(inout) :: a(n1,n2,nthds), b(n1,n2,nthds), c(n1,n2,nthds), d(n1,n2,nthds)
    
    integer :: i,j,ti
    double precision, allocatable, dimension(:,:) :: e
    double precision, allocatable, dimension(:) :: rr

    allocate(e(1:n1,1:n2),rr(1:n2))

    do ti=1,nthds
        do j=1,n2
            e(:,j)  =0.0d0
            e(2,j)  = -a(2,j,ti)
            e(n1,j) = -c(n1,j,ti)
        enddo

        do j=1,n2
            d(2,j,ti)=d(2,j,ti)/b(2,j,ti)
            e(2,j)   =e(2,j)   /b(2,j,ti)
            c(2,j,ti)=c(2,j,ti)/b(2,j,ti)
        enddo

        do j=1,n2
            do i=3,n1
                rr(j)=1.d0/(b(i,j,ti)-a(i,j,ti)*c(i-1,j,ti))
                d(i,j,ti)=rr(j)*(d(i,j,ti)-a(i,j,ti)*d(i-1,j,ti))
                e(i,j)=rr(j)*(e(i,j)-a(i,j,ti)*e(i-1,j))
                c(i,j,ti)=rr(j)*c(i,j,ti)
            enddo
        enddo

        do j=1,n2
            do i=n1-1,2,-1
                d(i,j,ti)=d(i,j,ti)-c(i,j,ti)*d(i+1,j,ti)
                e(i,j)=e(i,j)-c(i,j,ti)*e(i+1,j)
            enddo
        enddo

        do j=1,n2
            d(1,j,ti)=(d(1,j,ti)-a(1,j,ti)*d(n1,j,ti)-c(1,j,ti)*d(2,j,ti))/(b(1,j,ti)+a(1,j,ti)*e(n1,j)+c(1,j,ti)*e(2,j))
        enddo

        do j=1,n2
            do i=2,n1
                d(i,j,ti) = d(i,j,ti) + d(1,j,ti)*e(i,j)
            enddo
        enddo
    enddo

    deallocate(e,rr)

end subroutine tdma_cycl_many_cuda
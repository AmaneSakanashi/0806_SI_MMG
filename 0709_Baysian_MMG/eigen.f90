module eigen
  implicit none
  public :: mydger, mydsyr, mydgemm, mydgemv, myddot, mydsyev
  private :: pythag, tred2, tqli, assert_int_range, assert_int_lower, assert_int_upper
contains

  subroutine mydsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
    implicit none
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: n, lwork, lda
    double precision, intent(inout) :: a(lda, n) ! Fixed on Oct. 4th, 2017
    double precision, intent(out) :: w(n), work(lwork)
    integer, intent(out) :: info
    ! default lwork = 3 * n - 1
    call assert_int_lower( lda, n )
    call assert_int_lower( lwork, 3 * n - 1 )    
    call tred2(a(1:n, 1:n), w, work, n)
    call tqli(w, work, a(1:n, 1:n), n, info)
    ! BUG Fixed on Jan. 4th, 2017. Thanks to Mr. Miyagi. (info parameter added)
  end subroutine mydsyev

  subroutine mydgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character, intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    double precision, intent(in) :: alpha,beta    
    double precision, intent(in) :: a(:,:), b(:,:) 
    double precision, intent(inout) :: c(:,:)
    integer :: i, j
    if ( transa .eq. 'T') then
       if ( transb .eq. 'T') then
          do j = 1, n
             do i = 1, m
                c(i, j) = c(i, j) * beta + alpha * myddot( a(1:k, i), b(j, 1:k) )
             end do
          end do
       else
          do j = 1, n
             do i = 1, m
                c(i, j) = c(i, j) * beta + alpha * myddot( a(1:k, i), b(1:k, j) )
             end do
          end do
       end if
    else
       if ( transb .eq. 'T') then
          do j = 1, n
             do i = 1, m
                c(i, j) = c(i, j) * beta + alpha * myddot( a(i, 1:k), b(j, 1:k) )
             end do
          end do
       else
          do j = 1, n
             do i = 1, m
                c(i, j) = c(i, j) * beta + alpha * myddot( a(i, 1:k), b(1:k, j) )
             end do
          end do
       end if
    end if
  end subroutine mydgemm

  subroutine mydgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    implicit none    
    character, intent(in) :: trans
    integer, intent(in) :: lda,incx,incy,m,n
    double precision, intent(in) :: alpha, beta    
    double precision, intent(in) :: a(lda,n), x(:)
    double precision, intent(inout) :: y(:)
    integer :: i
    if ( trans .eq. 'N') then
       do i = 1, m
          y(i) = y(i) * beta + alpha * myddot( a(i, 1:n), x(1:n) )
       end do
    else
       do i = 1, n
          y(i) = y(i) * beta + alpha * myddot( a(1:m, i), x(1:m) )
       end do
    end if
  end subroutine mydgemv

  double precision function myddot(dx,dy)
    implicit none
    double precision, intent(in) :: dx(:),dy(:)
    myddot = dot_product(dx,dy)
    return
  end function myddot

  subroutine mydger(m,n,alpha,x,incx,y,incy,a,lda)
    implicit none
    double precision, intent(in) :: alpha
    integer, intent(in) :: incx,incy,lda,m,n
    double precision, intent(in) :: x(:),y(:)
    double precision, intent(inout) :: a(:,:)
    integer :: j    
    do j = 1, n
       a(1:m, j) = a(1:m, j) + (alpha * y(j)) * x(1:m)
    end do
  end subroutine mydger

  subroutine mydsyr(uplo,n,alpha,x,incx,a,lda)
    double precision, intent(in) :: alpha
    integer, intent(in) :: incx,lda,n
    character, intent(in) :: uplo
    double precision, intent(in) :: x(:)
    double precision, intent(inout) :: a(:, :)
    call mydger(n,n,alpha,x,1,x,1,a,lda)
  end subroutine mydsyr

  double precision function pythag( a, b )
    implicit none    
    double precision, intent(in) :: a, b

    double precision :: absa, absb

    absa = abs(a)
    absb = abs(b)
    if (absa > absb) then
       pythag = absa * sqrt(1.0d0 + (absb/absa) ** 2)
    else
       pythag = absb * sqrt(1.0d0 + (absa/absb) ** 2)
    end if
  end function pythag
  
  subroutine tred2( a, d, e, n )
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: a(n, n)
    double precision, intent(out) :: d(n)
    double precision, intent(out) :: e(n)

    integer :: l, j, i
    double precision :: scale, hh, h, g, f
    double precision :: gg(n)
    double precision :: agg(n, n)    

    do i = n, 2, -1
       l = i - 1
       h = 0.0d0
       if (l > 1) then
          scale = sum(abs(a(i, 1:l)))
          if (scale == 0.0d0) then
             e(i) = a(i, l)
          else
             a(i, 1:l) = a(i, 1:l) / scale
             h = sum(a(i, 1:l)**2)
             f = a(i, l)
             g = -sign(sqrt(h), f)
             e(i) = scale * g
             h = h - f * g
             a(i, l) = f - g
             a(1:l, i) = a(i, 1:l) / h
             do j = 1, l
                e(j) = (myddot(a(j, 1:j), a(i, 1:j)) &
                     + myddot(a(j+1:l, j), a(i, j+1:l))) / h
             end do
             f = myddot(e(1:l), a(i, 1:l))
             hh = f / (h + h)
             e(1:l) = e(1:l) - hh * a(i, 1:l)
             do j = 1, l
                a(j, 1:j) = a(j, 1:j) - a(i, j) * e(1:j) - e(j) * a(i, 1:j)
             end do
          end if
       else
          e(i) = a(i, l)
       end if
       d(i) = h
    end do
    d(1) = 0.0d0
    e(1) = 0.0d0
    do i=1, n
       l = i - 1
       if (d(i) /= 0.0d0) then
          gg(1:l) = 0.0d0 ! bug fixed on Jan 18, 2018          
          call mydgemv('T', l, l, 1.0d0, a(1:l, 1:l), l, a(i, 1:l), 1, 0.0d0, gg(1:l), 1)
          !gg(1:l) = matmul(a(i, 1:l), a(1:l, 1:l))
          do j = 1, l
             agg(1:l, j) = a(1:l, i) * gg(j)  
          end do
          a(1:l, 1:l) = a(1:l,1:l) - agg(1:l, 1:l)
       end if
       d(i) = a(i, i)
       a(i, i) = 1.0d0
       a(i, 1:l) = 0.0d0
       a(1:l, i) = 0.0d0
    end do
  end subroutine tred2
    
  subroutine tqli(d, e, z, n, info)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: d(:), e(:)
    double precision, intent(inout) :: z(:, :)
    integer, intent(out) :: info

    integer i, iter, l, m, ndum
    double precision, parameter :: max_iter = 30
    double precision :: b, c, dd, f, g, p, r, s
    double precision :: ff(n)

    info = 0
    ndum = n
    e(:) = eoshift(e(:), 1)
    do l = 1, n
       iter = 0
       iterate: do
          do m = l, n-1
             dd = abs(d(m)) + abs(d(m+1))
             if (abs(e(m)) + dd == dd) exit
          end do
          if (m == l) exit iterate
          if (iter == max_iter) then
             info = info + 1
             exit iterate
          end if
          iter = iter + 1
          g = (d(l+1) - d(l)) / (2.0d0 * e(l))
          r = pythag(g, 1.0d0)
          g = d(m) - d(l) + e(l) / (g + sign(r, g))
          s = 1.0d0
          c = 1.0d0
          p = 0.0d0
          do i = m-1, l, -1
             f = s * e(i)
             b = c * e(i)
             r = pythag(f, g)
             e(i+1) = r
             if (r == 0.0) then
                d(i+1) = d(i+1) - p
                e(m) = 0.0d0
                cycle iterate
             end if
             s = f / r
             c = g / r
             g = d(i+1) - p
             r = (d(i) - g) * s + 2.0d0 * c * b
             p = s * r
             d(i+1) = g + p
             g = c * r - b
             ff(1:n) = z(1:n, i+1)
             z(1:n, i+1) = s * z(1:n, i) + c * ff(1:n)
             z(1:n, i) = c * z(1:n, i) - s * ff(1:n)
          end do
          
          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0d0
       end do iterate
    end do
  end subroutine tqli

  subroutine assert_int_lower ( var, lower )
    implicit none
    integer, intent(in) :: var, lower
    if (lower > var) then
       print *, 'lower bound violation.'
    end if
  end subroutine assert_int_lower

  subroutine assert_int_upper ( var, upper )
    implicit none
    integer, intent(in) :: var, upper
    if (var > upper) then
       print *, 'upper bound violation.'
    end if
  end subroutine assert_int_upper

  subroutine assert_int_range ( var, lower, upper )
    implicit none
    integer, intent(in) :: var, lower, upper
    if (lower > var) then
       print *, 'lower bound violation.'
    else if (var > upper) then
       print *, 'upper bound violation.'
    end if
  end subroutine assert_int_range
  
end module eigen

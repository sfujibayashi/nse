subroutine locate(xx,n,x,j)
  implicit NONE
  integer :: j,n
  real(8) :: x,xx(n)
  integer :: jl,jm,ju
  
  jl=0                 !Initialize lower
  ju=n+1               !and upper limits.
10 if(ju-jl.gt.1) then  !If we are not yet done,
     jm=(ju+jl)/2      !compute a midpoint,
     if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
        jl=jm          !and replace either the lower limit
     else
        ju=jm          !or the upper limit, as appropriate.
     endif
     goto 10           !Repeat until
  endif                !the test condition 10 is satisfied.
  if(x.eq.xx(1)) then  !Then set the output
     j=1
  elseif(x.eq.xx(n)) then
     j=n-1
  else
     j=jl
  endif
  return               !and return.
end subroutine locate

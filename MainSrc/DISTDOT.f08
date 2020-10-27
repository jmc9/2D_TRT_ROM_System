real*8 function distdot(n,x,ix,y,iy)
integer n, ix, iy
real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)

CALL DDOT(n,x,ix,y,iy)

END FUNCTION

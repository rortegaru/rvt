      PROGRAM xqromb
C     driver for routine qromb
      REAL PIO2
      PARAMETER(PIO2=1.5707963)
      REAL a,b,fint,func,s,cl68_integrand,zlow,zup,result
      EXTERNAL fint,func,cl68_integrand
      a=0.0
      b=PIO2
      zlow=6
      zup=10
      write(*,'(1x,a)') 'Integral of FUNC computed with QROMB'
      write(*,'(1x,a,f10.6)') 'Actual value of integral is',
     *     fint(b)-fint(a)
c     call qromb(func,a,b,s)
      call qromb(cl68_integrand,zlow,zup,result)
      write(*,'(1x,a,f10.6)') 'Result from routine QROMB is',result
      END

      REAL FUNCTION func(x)
      REAL x
      func=(x**2)*(x**2-2.0)*sin(x)
      END

      REAL FUNCTION fint(x)
C     integral of FUNC
      REAL x
      fint=4.0*x*((x**2)-7.0)*sin(x)-((x**4)-14.0*(x**2)+28.0)*cos(x)
      END


*----------------- BEGIN CL68_INTEGRAND -----------------------------
      REAL FUNCTION cl68_integrand(z)
* Dates: 06/09/95 - Written by D.M. Boore.  See 7/11/82 notes for
*                   stochastic model, with 6/9/95 addition that uses
*                   a variable transformation to remove the sqrt
*                   singularity at the origin.
*        01/03/95 - Name changed from cl_int to cl68_integrand
*        02/05/00 - Made changes suggested by C. Mueller to avoid
*                   possible numerical problem when "an" is large

      real xi, an, y, z
      PARAMETER (xi=1.0)
      PARAMETER (an=1.0)
*      cl68_integrand = 2.0*(1.0-(1.0-xi*exp(-z**2))**an)  ! original
      y = 1.0 * alog(1.0-1.0*exp(-z**2))   ! Mueller modification
      cl68_integrand = 2.0*(1.0-exp(y))  ! Mueller modification

      return
      end
*----------------- END CL68_INTEGRAND -----------------------------
      


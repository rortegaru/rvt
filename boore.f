
c-----
c	David Boore's smsim code extracts -- these are used since they
c	create the DEFINED random time series with reasonable 
c	displacements
c-----
C*----------------- BEGIN Acc2VD -----------------------------
      subroutine acc2vd(acc, npts, dt, rmv_trnd, v0, d0, vel, dis)


C* Compute velocity and displacement time series from acceleration,
C* assuming that the acceleration
C* is represented by straight lines connecting the digitized values.

C* Dates: 02/09/99 - Written by D.M. Boore
C*        01/07/00 - Added initial velocity and displacements (v0, d0).
C*                   This also requires the addition of a linear trend in
C*                   displacement.
C*                   Also, Bill Joyner, Chris Stephens, and I considered how to
C*                   handle the first point.  Assuming v0 = 0, BAP uses a 
C*                   trapezoidal rule such that the first vel(1) = 0.5*dt*a(1),
C*                   but in acc2vd, vel(1) = 0.  In effect, BAP starts the
C*                   integration at -dt rather than 0.0, with the result that
C*                   because vel(1) depends on dt, vel(1) will change if the
C*                   digitization interval is changed.  We agreed that this is 
C*                   not correct.

      real acc(npts), vel(npts), dis(npts)
      logical rmv_trnd
      double precision cumv, cumd, a1, a2, v1, 
     1 ddt, ddt_2, ddtdt_6

      if (rmv_trnd) then      
C* remove trend first (straight line between first and last points)
C* Note: acc is replaced with detrended time series
         call rmvtrend(acc, npts)
      endif

C* compute velocity and displacement, using analytical formulas based
C* on representing the acceleration as a series of straightline segments.

      ddt     = dble(dt)
      ddt_2   = dble(dt/2)
      ddtdt_6 = dble(dt**2/6)

      cumv = 0.0
      cumd = 0.0

      vel(1) = sngl(cumv) + v0
      dis(1) = sngl(cumd) + d0
      do j=2,npts
        a1 = acc(j-1)
        a2 = acc(j)
        v1 = vel(j-1)
        cumv = cumv + (a1 + a2)*ddt_2
        vel(j) = sngl(cumv) + v0
        cumd = cumd + v1*ddt + (2.0*a1 + a2)*ddtdt_6
        dis(j) = sngl(cumd) + d0  ! no linear trend neede; it's include in v1
      end do

      return
      end
C*----------------- END Acc2VD -----------------------------

C* ----------------------------- BEGIN RMVTREND ----------------
      subroutine rmvtrend(y, n)

C* Removes a straightline fit to first and last points, replacing
C* the input array with the detrended array

C* Dates: 02/09/99 - written by D. Boore


      real y(*)

      y1 = y(1)
      y2 = y(n)
      slope = (y2 - y1)/float(n-1)

      do i = 1, n
        y(i) = y(i) - (y1 + slope*float(i-1))
      end do

      return
      end
C* ----------------------------- END RMVTREND ----------------
C*----------------- BEGIN GET_ACC -----------------------------
      subroutine getacc(acc)
C* Dates:  06/07/95 - Written by D.M. Boore
C*         06/08/95 - Normalize by sqrt(avg square amp)
C*         06/09/95 - Get r, amag from common, not parameter list
C*         08/08/95 - Removed dc from noise segment before going to
C*                    the frequency domain (in get_acc).
C*         08/08/95 - Changed things so that remove dc from noise sample before
C*                    applying the window (in get_acc).
C*         10/17/95 - Remove dc or not from random series, depending on a flag
C*                    (irmvdc) that is passed through the include statement
C*         11/14/95 - Added call to const_am0_gsprd
C*         11/16/95 - Replaced 'shape = ...' statement with call to spect_amp
C*         12/06/95 - Major simplification in exponential window, removing
C*                    tapers (the window itself provides tapers) and correcting
C*                    a bug so that the window can be used for a duration
C*                    other than the duration used to determine the window 
C*                    parameters.
C*         12/22/95 - Removed rmean real part of work only
C*         01/22/96 - Use REALFT rather than FORK
C*         02/10/99 - Transfer choice of npts from rvtdsubs.for to this
C*                    subroutine.  This also includes pad extension if a low-cut
C*                    filter is used.  
C*                    In all cases 
C*                    npts*dt >= max(tshift, tfpad) + dur_fctr * nmotion * dt 
C*                               + tfpad 
C*         03/12/99 - Suggest increasing dt if the time series is too long.
C*         08/03/99 - Added uniform deviate if desired

c-----
c	CHANGES BY RBH AUG 01 2001
c	indent, comments, Boore's is in lower case, RBH uppercase
c-----

	INTEGER N2
	PARAMETER (N2=131072)
	real work(N2)
	real acc(N2)
	REAL windbox, windexp

	include 'smsim.fi'

c-----
c	Set spectral parameters:
c 	call this now because need corner frequencies
c                            for window durations
c-----	
	call spect_scale()  
c-----
c	Fill an array with the proper duration of windowed Gaussian noise:
c-----


	do 1000 i = 1, N2
		acc(i)=0.
		work(i) = 0.
 1000	continue

c-----
c	Modifications in order to use another window 5/30/1995
c	Calculate the number of points in the noise sample:
c-----

c----RBH
c	GET DURATION FOR THE PATH
c-----
	CALL GETDUR(PATHDUR,R)

	if(indxwind .eq. 0) then        
c-----
c	 	BOX WINDOW
c-----
		durex = dursource(w_fa, fa, w_fb, fb)+
     1			PATHDUR
		nmotion = durex/dt
		ntaper = ifix(taper*nmotion)  
c-----
c		Increase duration of motion by tapers 
c		(figure 1/2 front and back):
c		taper is fractional, not percent
c-----
		nmotion = nmotion + ntaper  

	else                           
c-----
c		EXPONENTIAL WINDOW
c-----
		durex = dursource(w_fa, fa, w_fb, fb)+
     1			PATHDUR
		nmotion = 2.0 * durex/dt   
c-----
c	doubled following Boore (1983, p. 1869)
c-----
		tmotion = nmotion * dt
		tw = tmotion * twdtmotion
	endif

c-----
c	Calculate nstart and nstop, depending on tshift and if a
c	low-cut filter is used.
c-----
	if (fcut .eq. 0.0) then
		tfpad = 0.0
	else
c-----
c		(Converse, USGS OFR 92-296A, p. 2-3)
c----- 
		tfpad = 1.5 * (norder/2) / fcut 
	endif

	if (tfpad .gt. tshift) tshift = tfpad
	nstart = tshift/dt + 1
	nstop = nstart + nmotion

c-----
c	Calculate npts, depending on tshift, dur_fctr, and if a
c	low-cut filter is used.
c	compute smallest power of two for which the resulting duration 
c		is greater
c	than or equal to tsimdur:
c-----

	tsimdur = nstart*dt + dur_fctr * nmotion * dt + tfpad

	npts = 2.0**ifix(alog10(tsimdur/dt + 1.0)/alog10(2.0))
	if ( (npts-1)*dt .lt. tsimdur) npts = 2 * npts
	if (npts .gt. N2) then
        	write(*,'(2x,a,i5,a,i5,a)') 
     1       	' *** FATAL ERROR: NPTS (',
     1                                    npts, 
     1       	') > NDIMEN (', 
     1                                     ndimen, 
     1       	'); QUITTING!'
   	     write(*, '(2x,a)') ' TRY INCREASING DT'
       		stop
	endif

c-----
c	Generate the Gaussian white noise with zero mean, unit variance:
c-----

	if (iran_type .eq. 0) then
		do 1010 i = nstart, nstop
			work(i) = gasdev(iseed)
 1010		continue
	else
		do 1020 i = nstart, nstop
			work(i) = ran1(iseed) - 0.5
 1020		continue
	endif
	if (irmvdc .eq. 0) then
		rmean = 0.0
	else
		call mean(work, nstart, nstop, rmean)
	endif

c-----
c	Window the white noise:
c-----

	if(indxwind .eq. 0) then        
c-----
c		BOX WINDOW
c-----

		do 1040 i = nstart, nstop
			work(i) = windbox(i,nstart,nstop,ntaper)
     1			* (work(i) - rmean)
 1040		continue

	else                           
c-----
c		EXPONENTIAL WINDOW
c-----

		do 1050 i = nstart, nstop
			t = float(i-nstart) * dt
			work(i)=windexp(t,tw,eps_wind,eta_wind,new_mr)
     1				* (work(i) - rmean)
 1050		continue
	endif

c-----
c	Transform to frequency domain:
c-----
	call realft(work,npts,+1)

c-----
c	Find sqrt of squared spectral amplitude for normalization:
c-----
	call avrlft(work,npts,avgsqamp) 
	sqrt_avgsqamp = sqrt(avgsqamp)

c-----
c	Get frequency independent parameter:
c-----
	call const_am0_gsprd()     
c-----
c	returns freq_indep_factor through common
c-----
	iaorins = 1             
c-----
c	ground motion
c-----
	idva = 2                
c-----
c	acceleration
c-----

	df = 1.0/(float(npts)*dt)
	do 1060 i = 2,npts/2
		f = (i-1) * df
		sfact = spect_amp(f)/ sqrt_avgsqamp
		work(2*i-1) = sfact * work(2*i-1)
		work(2*i)   = sfact * work(2*i)
 1060	continue
	work(1) = 0.0              
c-----
c	OK for acceleration spectrum, but not for
c	displacement spectrum
c-----
	fnyq = 1.0/(2.0*dt)
	work(2) = spect_amp(fnyq) * work(2)/ sqrt_avgsqamp

c-----
c	Transform back to time domain:
c-----
	call realft(work,npts,-1)

c-----
c	Apply FFT normalization:
c-----
	afact = 2 * df
	do 1070 i = 1, npts
		acc(i) = afact * work(i)
 1070	continue
	return
	end
C*----------------- END GET_ACC -----------------------------
*----------------- BEGIN MEAN -----------------------------
      subroutine mean(work, nstart, nstop, rmean)
*  Dates: 01/22/96 - Written by D. Boore
      real work(*)
* find mean of the array:
      sum = 0.0
      do i = nstart, nstop
        sum = sum + work(i)
      end do
      rmean = sum/float(nstop-nstart+1)
      return
      end
*----------------- END MEAN -----------------------------
*------------------- BEGIN AVGSQ_REALFT -------------------------------
* Dates: 01/22/96 - Written by D.M. Boore
      subroutine avrlft( s, npts, avgsqamp)
      real s(*)
      sum=0.
      do j = 2, npts/2                   !don't include the dc or Nyquist values
        sum=sum + s(2*j-1)**2 + s(2*j)**2 ! odd, even = real, imag spect
      end do
      avgsqamp = sum/float(npts/2 - 1)
      return
      end
*  ------------------- END AVGSQ_REALFT -------------------------------
*  ------------------- BEGIN WIND_BOX -------------------------------
      function windbox( i, nstart, nstop, ntaper)
c
c applies cosine tapered window.
c unit amplitude assumed
c
c written by D. M. Boore
c
c latest revision: 9/26/95
      real windbox

      windbox = 0.0
      if ( i .lt. nstart .or. i. gt. nstop) return
      windbox = 1.0
      if ( i .ge. nstart+ntaper .and. i .le. nstop-ntaper ) return
c
      pi = 4.0 * atan(1.0)
c
      dum1 = (nstop+nstart)/2.0
      dum2 = (nstop-nstart-ntaper)/2.0
c
      windbox = 0.5 * (1.0 - sin( pi*
     *  ( abs(float(i)-dum1) - dum2 ) /float(ntaper) ) )
      return
      end
*  ------------------- END WIND_BOX -------------------------------

*  ------------------- BEGIN WIND_EXP ---------------------------
      function  windexp( t, tw, epswind, etawind, new_mr)
c
c     apply Sargoni and Hart (1974) window, with parameters
c     tw, eps (fraction of tw to reach peak), and
c     eta ( fraction of peak ampl. at tw).  See Boore (BSSA, 1983,
*     p. 1869).  Note that t can be larger than tw.

* Dates:
*         05/30/95 - Initially written by Dave Boore, based on a routine
*                    by G. Atkinson but with significant structural 
*                    changes and the taper to be a raised cosine taper
*         12/06/95 - Removed the taper and simplified the parameters in the
*                    calling list.  Also added a switch so that the window
*                    shape coefficients are only computed once.  This assumes
*                    that the shape coefficients will be the same for
*                    any run of the program... a good assumption.
*         12/28/95 - Used new_mr, set from driver, to control computation
*                    of b,c,a; before I used "firstcall", set in this 
*                    subprogram, but this gave an error if the driver
*                    looped over multiple values of m and r.

      real windexp
      logical new_mr
      save a, b, c

      if (new_mr) then
        b = -epswind * alog(etawind)/
     :      (1. + epswind*(alog(epswind)-1.))
        c = b/(epswind*tw)
        a = (exp(1.0)/(epswind*tw))**b
        new_mr = .false.
      endif

      windexp = 0.
      if( t .lt. 0.0) return
c
c     Apply Sargoni and Hart window.
c
      windexp = a*t**b * exp(-c*t)

      return
      end
*------------------- END WIND_EXP ---------------------------
*----------------- BEGIN SPECT_AMP -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
      function spect_amp(f)
C-----
C	RBH COMMON
C-----
	COMMON/BUTFLT/FL,FH,NBUT,DDTT,DOBUT,IPOW
	REAL FL, FH, DDTT
	INTEGER NBUT,IPOW
	LOGICAL DOBUT
	COMPLEX S, ZAMP

      include 'smsim.fi'

      tempf = 1.0
      if ( idva .ne. 0 ) tempf = (twopi*f)**float(idva)
	CALL GETSIT(SITEAMP, F)
      spect_amp =   freq_indep_factor * ( tempf ) * 
     :      buttrlcf(f, fcut, norder) * 
     :      spect_shape(f, fa, fb, pf, pd, am0b_m0, numsource) * 
     :      SITEAMP *
     :      dimin(f)
                                 ! Could save some multiplications
                                 ! by removing freq_indep_factor outside the
                                 ! function, but at the cost of possible
                                 ! confusion.

	IF(IAORINS.EQ.1)THEN
		H = 1
	ELSE IF(IAORINS.EQ.2)THEN
c-----
c	converts from displacement 
c	to pseudo velocity response for IPOW = 0
c
c	IPOW = 0 is pseudo velocity of mass -> PSV
c	IPOW = 1 is velocity of mass ->SV
c	IPOW = 2 is acceleration of mass ->Sa
c-----
		wn = twopi * fosc      
		h = harmoscf( f, fosc, damp, idva )
		IF(IPOW.EQ.0)THEN
			H = H * WN
		ELSE IF(IPOW.EQ.1)THEN
			H = H * TWOPI * F
		ELSE IF(IPOW.EQ.2)THEN
			H = H * TWOPI * F * TWOPI * F
		ENDIF
	ELSE IF(IAORINS.EQ.3)THEN
		ZAMP = CMPLX(1.0,0.0)
		FAC = 6.2831853*F*DDTT
C-----
C		     i omega dt
C		z = e
C-----
		S = CMPLX(COS(FAC),SIN(FAC))
		ZAMP = CMPLX(1.0,0.0)
		CALL ZPASSF(NBUT,2,.FALSE.,FL,FH,S,ZAMP,DDTT,1)
		CFAC = CABS(ZAMP)
		H = CFAC
	ENDIF
      spect_amp = h * spect_amp   ! spect_amp contains the spectral amplitude.

      return
      end
*----------------- END SPECT_AMP -----------------------------

*----------------- BEGIN CONST_AM0_GSPRD -----------------------------
* Dates: 11/14/95 - Written by D.M. Boore
*        07/02/99 - Added magnitude-dependent "depth" from Atkinson
*                   and Silva, which required adding some parameters to
*                   the passed arguments in gsprd
*        06/08/00 - Moved computation of const into this subroutine

      subroutine const_am0_gsprd()

      include 'smsim.fi'

* Define constant, for r=r_ref(km).Usually set r_ref = 1.0 km.  Be careful
* if a different value or different units are used.  In particular, using 
* different units will require a change in the scaling factor of 1.e-20 below

      const=prtitn*rtp*fs*(1.e-20)/(4.*pi*rho*beta**3*r_ref)

c-----RBH
c
C      freq_indep_factor = const*am0*
C     :  gsprd(r, r_ref, nsprd_segs, rlow, slope, numsource, amag)
C*                         (const from Get_Params, am0 from Spect_Scale) 
C-----
	CALL GETGEO(gr,r)
      freq_indep_factor = const*am0*
     :  gr
*                         (const from Get_Params, am0 from Spect_Scale) 

      return
      end
*----------------- END CONST_AM0_GSPRD -----------------------------

*----------------- BEGIN GSPRD -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
*        07/02/99 - Added magnitude-dependent "depth" from Atkinson
*                   and Silva, which required adding some parameters to
*                   the passed arguments
*        06/05/00 - Added some explanation of r
*        06/08/00 - Make gsprd nondimensional through the use of r_ref, which 
*                   now appears in the definition of variable const
*                   in const_am0_gsprd

      function gsprd(r,r_ref,nsprd_segs,rlow,slope,numsource,amag)
      real r_ref, rlow(*), slope(*), geff(10)
* Note that generally r = hypocentral distance.  For Atkinson and Silva 
* (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
* their paper; below their eq. 4), so that rmod is, in effect, accounting
* source depth twice 
      rmod = r
      if (numsource .eq. 9) then            ! Atkinson and Silva (2000)
        deff = 10.0**(-0.05 + 0.15 * amag)
        rmod = sqrt(r**2 + deff**2)
      endif
      geff(1) = r_ref/rlow(1)  ! usually set r_ref = 1.0 km.  Be careful
                               ! if a different value or different units are
                               ! used.  In particular, using different units
                               ! will require a change in the scaling factor
                               ! of 1.e-20 used in the definition of const in
                               ! const_am0_gsprd

      do i = 2, nsprd_segs
        geff(i) = geff(i-1)*(rlow(i)/rlow(i-1))**slope(i-1)
      end do
      if (rmod .le. rlow(1)) then
        j = 1
      else if (rmod .ge. rlow(nsprd_segs)) then
        j = nsprd_segs
      else
        call locate(rlow, nsprd_segs, rmod, j)
      endif
      gsprd = (geff(j)) * (rmod/rlow(j))**slope(j)

      return
      end
*----------------- END GSPRD -----------------------------
*----------------- BEGIN DIMIN -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
*        07/02/99 - Added modification to r required by Atkinson
*                   and Silva (1999)
*        06/05/00 - Substitute c_q for beta in akappaq and add comments
*                   about r

      function dimin(f)
      real dimin, mag
      include 'smsim.fi'

* Note that generally r = hypocentral distance.  For Atkinson and Silva 
* (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
* their paper; below their eq. 4), so that rmod is, in effect, accounting
* source depth twice 
      rmod = r
      if (numsource .eq. 9) then            ! Atkinson and Silva (2000)
        deff = 10.0**(-0.05 + 0.15 * amag)  
        rmod = sqrt(r**2 + deff**2)
      endif

      akappaq = rmod/(c_q*q(f))

      mag = amag    
      dimin = exp( -pi*(kappa_f(mag) + akappaq) * f)/
     :   sqrt( 1. + (f/fm)**8.0 )
      

      return
      end
*----------------- END DIMIN -----------------------------
*----------------- BEGIN KAPPA_F -----------------------------
* Dates: 02/28/97 - Written by D.M. Boore
      function kappa_f(mag)
      real mag
      include 'smsim.fi'
  
      kappa_f = akappa + dkappadmag*(mag-amagkref)

      return
      end
*----------------- END KAPPA_F -----------------------------
        
*----------------- BEGIN Q -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
*        12/14/95 - Added check for equal transition frequencies
      function q(f) 
      logical firstcall 
      save qt1, qt2, st
      include 'smsim.fi'
      firstcall=.true.
      q = 9999.0
      if (f .eq. 0.0) return
        
      if (firstcall) then
        qt1 = qr1*(ft1/fr1)**s1
        qt2 = qr2*(ft2/fr2)**s2
        st = 0.0
        if (ft1 .ne. ft2) then
          st = alog10(qt2/qt1)/alog10(ft2/ft1)
        endif
        firstcall = .false.
      endif
      
      if ( f .le. ft1) then
        q = qr1*(f/fr1)**s1
      else if ( f .ge. ft2) then
        q = qr2*(f/fr2)**s2
      else
        q = qt1*(f/ft1)**st
      endif

      return
      end
*----------------- END Q -----------------------------
*  ------------------- BEGIN BUTTRLCF -------------------------------
* Dates: 06/07/95 - Written by D.M. Boore
      function buttrlcf( f, fcut, norder)
c
c Computes the response of an norder, bidirectional
* high-pass Butterworth filter.  This is the filter
* used by the AGRAM processing (the equation was
* provided by April Converse).

* Modification: 3/27/89 - created by modifying HiPassF

      real buttrlcf
      buttrlcf = 1.0
      if ( fcut.eq.0.0 ) return

      buttrlcf = 0.0

      if ( f .eq. 0.0) return

      buttrlcf = 1.0/ (1.0+(fcut/f)**(2.0*norder))

* Note: Because this filter is intended to simulate the
* two-pass, zero-phase (acausal) Butterworth filter used in
* accelerogram processing, buttrlcf = 1/2 when f = fcut, not 1/sqrt(2) as in
* a single-pass, causal Butterworth filter.

      return
      end
*  ------------------- END BUTTRLCF -------------------------------
*----------------- BEGIN SPECT_SHAPE -----------------------------
* Source Displacement Spectrum
* Dates: 06/07/95 - Written by D.M. Boore
*        11/15/95 - changed source types
*        12/02/96 - added Atkinson and Silva (model 4) (I did this earlier)
*                   and Haddon (model 5)
*        02/28/97 - Added Atkinson's new version of the source spectra
*                   derived from Atkinson and Silva (this will appear
*                   in Atkinson and Boore... Evaluation of source models...).
*                   (model 6)
*        06/24/97 - added Boatwright and Choy scaling (model 7).
*        07/21/97 - Added Joyner ENA model (model 8; the spectral shape is
*                   the same as his model 2, but I because the corner frequency
*                   relations are new I have to repeat the shape with the new 
*                   model number).
*        09/02/98 - Renamed model 6 to AB98-Ca to be consistent with usage
*                   in Tables 3 and 4 in Atkinson and Boore, BSSA, v. 88, 
*                   p. 917--934.
*        02/16/99 - Added Atkinson and Silva, 1999, as model 9
*        06/05/00 - Changed "AS99" to "AS2000" because the paper was published
*                   in 2000 (BSSA 90, 255--274)

      function spect_shape(f, fa, fb, pf, pd, am0b_m0, numsource)
      real spect_shape
      goto (1, 2, 3, 4, 5, 6, 7, 8, 9, 10), numsource

      write(*, '(a, i5, a)') ' !!!!!! numsource = ',
     : numsource, ', .ne. legal value; quitting !!!!!!'
      stop
 
* Single corner frequency:
1     sb = 1.0
      sa = 1.0/( 1.0 + (f/fa)**pf )**pd
      go to 900

* Joyner model
2     sb = 1.0/ ( 1.0 + (f/fb)**2 )**0.25
      sa = 1.0/ ( 1.0 + (f/fa)**2 )**0.75
      go to 900

* Atkinson 1993 model
3     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

* Atkinson & Silva 1996 (same format as Atkinson 1993) model
4     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

* Haddon (see 12/02/96 notes; approximate spectra in Fig. 10 of
* Haddon's paper in BSSA, v. 86, p. 1312)
5     pda = 1.0/8.0
      pdb = 1.0/8.0
      pfa = 1/pda
      pfb = 1/pdb
      sa = 1.0/( 1.0 + (f/fa)**pfa )**pda
      sb = 1.0/( 1.0 + (f/fb)**pfb )**pdb
      go to 900

* AB98-Ca (Atkinson & Boore 1998) (same format as Atkinson 1993) model
6     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

* Boatwright and Choy (this is the functional form used by 
*  Boore and Atkinson, BSSA 79, p. 1761)
7     sa = 1.0
      if (f .ge. fa) sa = fa/f
      sb = 1.0/sqrt( 1.0 + (f/fb)**2 )
      go to 900 

* Joyner model (to be used with his ENA two-corner model)
8     sb = 1.0/ ( 1.0 + (f/fb)**2 )**0.25
      sa = 1.0/ ( 1.0 + (f/fa)**2 )**0.75
      go to 900

* 
* AS2000 (Atkinson and Silva, 2000, BSSA 90, 255--274) 
* (same format as Atkinson 1993) model
9     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

* For customized relation:
10    sb = 1.0
      sa = 1.0
      go to 900

 
900   continue
      spect_shape = sa*sb

      return
      end
*----------------- END SPECT_SHAPE -----------------------------
*----------------- BEGIN HARMOSCF -----------------------------
      function harmoscf( f, fosc, damp, idva)
c harmonic oscillator displacement response. 
* idva = 0 for response to displacement
* idva = 2 for response to acceleration
* idva = 1 returns 0 for the response
* The response is normalized to be unity in the flat portion.

* Written by D. Boore, 12/01/83
c latest modification:  3/25/89
*                       7/30/00 - Changed dum for both cases of idva
*                                 (see notes from same date).

      pi = 4.0*atan(1.0)
      twopi = 2.0 * pi

      if (idva .eq. 0) dum = (twopi*f)**2/twopi**2
      IF (IDVA .EQ. 1) DUM = 0.0
      if (idva .eq. 2) dum =          1.0/twopi**2

      harmoscf = dum/sqrt( ( fosc*fosc - f*f )**2
     * + ( 2.0*f*fosc*damp )**2 )
      return
      end
*----------------- END HARMOSCF -----------------------------
*----------------- BEGIN DURSOURCE -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
      function dursource(w_fa, fa, w_fb, fb)
      real dursource
      dursource = w_fa/fa + w_fb/fb
      return
      end
*----------------- END DURSOURCE -----------------------------
*----------------- BEGIN SPECT_SCALE -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
*        06/05/96 - Added modified Atkinson & Silva scaling
*        12/02/96 - added Haddon scaling (see spect_shape.for comments)
*        02/28/97 - added Atkinson and Boore scaling 
*                   (see spect_shape.for comments)
*        06/24/97 - added Boatwright and Choy scaling
*        07/10/97 - changed A93, AS96, AS97 scaling to be constant
*                   stress for M < Mc, where Mc is the magnitude for which
*                   am0b_m0 = 1.0.  The single corner model for smaller
*                   magnitudes is determined so that the high frequency
*                   level is matched at M = Mc.
*        07/21/97 - Added Joyner 2-corner model for ENA, as specified 
*                   in his notes prepared for the SSHAC workshop (published 
*                   in SSHAC report, NUREG/CR-6372, p. B-303--B-305, 1997).
*        08/06/97 - I discovered that Joyner erroneously fit vertical spectral
*                   amplitudes.  He provided a revised model, fitting the
*                   horizontal amplitudes.  I changed model 8 accordingly.
*        09/02/98 - Renamed model 6 to AB98-Ca to be consistent with usage
*                   in Tables 3 and 4 in Atkinson and Boore, BSSA, v. 88, 
*                   p. 917--934.
*        02/16/99 - Added Atkinson and Silva (1999)
*        06/05/00 - Changed "AS99" to "AS2000" because the paper was published
*                   in 2000 (BSSA 90, 255--274)

	subroutine spect_scale()

	include 'smsim.fi'

	am0 = 10.**(1.5*amag + 16.05)
	am0b_m0 = 0.0

	IF(NUMSOURCE.EQ.1)then
c-----
c	Single corner frequency:
c-----
		stress = stressc*10.0**(dlsdm*(amag-amagc))
		fa = (4.906e+06) * beta * (stress/am0)**(1.0/3.0)
		fb = fa
		w_fa = 1.0
		w_fb = 0.0

	ELSE IF(NUMSOURCE.EQ.2)then
c-----
c	Joyner scaling:
c-----
		am0c = 10.0 ** ( 1.5*amagc + 16.05 )
		stress = stressc*10.0**(dlsdm*(amag-amagc))
		rat = stress/am0
		dum = 4.906e+06
		if ( am0 .gt. am0c ) rat = stress/am0c
		fb = ( dum*beta ) * ( fbdfa )**(3./4.) * 
     1			( rat )**(1./3.)
		fa = ( dum*beta ) * (stress)**(1./3.) * 
     1			(am0c)**(1./6.) *
     *			( fbdfa )**(-0.25) * ( am0 )**(-0.5)
		if ( am0 .lt. am0c ) fa = fb / fbdfa

	ELSE IF(NUMSOURCE.EQ.3)then
c-----
c		Atkinson 93 scaling:
c-----
		if (amag .gt. 4.0) then
			fa = 10.0**(2.41 - 0.533 * amag)
			fb = 10.0**(1.43 - 0.188 * amag)      
			am0b_m0 = 10.0**(2.52 - 0.637 * amag)
		else
c-----
c			 fa = fb for M = 2.84
c-----
			fb = 10.0**(2.678 - 0.5 * amag)
			fa = fb
			am0b_m0 = 1.0
		endif
		w_fa = 0.5
		w_fb = 0.0
	ELSE IF(NUMSOURCE.EQ.4)then

c-----
c		Atkinson and Silva 96 scaling, with am0b_m0 
c		modified by D. Boore on 6/04/96
c-----
		if (amag .gt. 4.6) then
			fa = 10.0**(2.181 - 0.496 * amag)
			fb = 10.0**(1.778 - 0.302 * amag)   
c-----
c		DMB's fitting of spctral ratios
c-----
			am0b_m0 = 10.0**(3.440 - 0.746 * amag)  
c-----
c		 in Atkinson & Silva preprint
c	       am0b_m0 = 10.0**(2.764 - 0.623 * amag) 
c-----
		else
c-----
c			 fa = fb for M = 2.08
c-----
			fb = 10.0**(2.689 - 0.5 * amag)
			fa = fb
			am0b_m0 = 1.0
		endif
	ELSE IF(NUMSOURCE.EQ.5)then

c-----
c		Haddon scaling:
c-----
		fa = 10.0**(0.3 - (1.5/3)*(amag-4.0))
		fb = 10.0**(1.4 - (1.5/3)*(amag-4.0))  
		w_fa = 0.5
		w_fb = 0.0
c-----
c		 fa < fb for all M
c-----
	ELSE IF(NUMSOURCE.EQ.6)then

c-----
c		AB98-Ca (Atkinson and Boore 98 scaling, 
c		based on fitting a model to the 
c		Atkinson and Silva 1997 Fourier amplitude spectra 
c		for California; see Atkinson and Boore, BSSA, v. 88, 
c		p. 917--934).
c-----
		if (amag .gt. 4.8) then
			fa = 10.0**(2.181 - 0.496 * amag)
			fb = 10.0**(1.308 - 0.227 * amag)      
			am0b_m0 = 10.0**(3.223 - 0.670 * amag)
		else
c-----
c			 fa=fb for M = 3.25
c-----
			fb = 10.0**(2.617 - 0.5 * amag)
			fa = fb
			am0b_m0 = 1.0
		endif
		w_fa = 0.5
		w_fb = 0.0
	ELSE IF(NUMSOURCE.EQ.7)then

c-----
c		Boatwright and Choy (this is not from Boore and 
c		Atkinson, BSSA 79, p. 1761;  it is based on new 
c		fits by DMB on 6/24/97 to data in Boat.&Choy, 1992 BSSA.
c		See BC_BSSA.WQ1 in \haddon subdirectory and 
c		handwritten notes on  yellow sheets.
c		except set fa=fb=constant stress scaling for M<=5.3)
c-----
		fa = 10.0**(3.409 - 0.681 * amag)
		fb = 10.0**(1.495 - 0.319 * amag)
		if (amag .le. 5.3) then
				fa = 0.634*10.0**(0.5*(5.3 - amag)) 
c-----
c			 0.634 = 10^(logfa+logfb)/2 at M5.3
c-----
				fb = fa
		endif
		w_fa = 0.5
		w_fb = 0.0
	ELSE IF(NUMSOURCE.EQ.8)THEN
c-----
c		Joyner ENA scaling:
c-----
		fa = 10.0**(2.312 - 0.5 * amag)
		fb = fa
	ELSE IF(NUMSOURCE.EQ.9)THEN
c-----
c		AS2000 -- Atkinson and Silva, 2000 scaling, based on fitting 
c		a point source model to finite source calculations, 
c		with constraints on various modeling studies, 
c		with modification for very small magnitude (when eps = 1) 
c-----
		if (amag .gt. 2.4) then
			fa = 10.0**(2.181 - 0.496 * amag)
			fb = 10.0**(2.41  - 0.408 * amag)     
			am0b_m0 = 10.0**(0.605 - 0.255 * amag)
		else
c-----
c			fa=fb for M = 3.25
c-----
			fb = 10.0**(1.431 - 0.5 * (amag - 2.4))
			fa = fb
			am0b_m0 = 1.0
		endif
		w_fa = 0.5
		w_fb = 0.0

c-----
c	 For customized scaling:
c-----
	ELSE IF(NUMSOURCE.EQ.10)THEN
		fa = 0.0
		fb = 0.0
	ELSE
      write(*, '(a, i5, a)') ' !!!!!! numsource = ',
     : numsource, ', .ne. legal value; quitting !!!!!!'
      stop
	ENDIF

	end
*----------------- END SPECT_SCALE -----------------------------
*----------------- BEGIN GET_MOTION -----------------------------
      subroutine get_motion(gmsim)

* Dates: 06/06/95 - written by D.M. Boore, adapted from main.for, which see for
*                   code for printing out much more info.
*        06/07/95 - Changed to two term asymptotic expression for output.  The
*                   'exact' causes problems.
*        06/08/95 - Changed anz to have a minimum value of 1.33 (as mentioned
*                   by Toro in section 3, vol. 3, appendices B & C, Seismic
*                   Hazard Methodology for Nuclear Facilities in the Eastern
*                   United States, EPRI, April 30, 1985, eq. B-35, p. B-80.)
*        06/08/95 - Changed names of variables such as 'anclee' to more
*                   meaningful names, and also put the sqrt 2 into the proper
*                   place rather than carrying it around with rms (in 
*                   'afact' in subroutine main on the VAX).
*        06/09/95 - Changed method of computing 'exact' solution to a numerical
*                   integration.
*        10/17/95 - Eliminated exact series evaluation, as well as the switch
*                   from asymptotic to "exact"; I also eliminated the commented
*                   statements that included Toro's "clumping" correction...
*                   the source code for that has been retained in smsmclmp.for.
*        12/14/95 - Pass freq20, durex through common block in smsim.fi
*        12/19/95 - Added zup to exact_numrcl (as of 1/3/96, cl68_numrcl_int)
*        12/25/95 - Add Herrmann's integration of C&L-H eq. 6.4
*        12/26/95 - Pass eps_rv and ane through smsim.fi common rather than
*                   through the parameter list
*        12/28/95 - Changed variable names to indicate that equation 6.8 of
*                   Cartwright and Longuet-Higgins is being used.
*        12/30/95 - Remove Herrmann's integration of C&L-H eq. 6.4
*        01/14/99 - Add calculation of Arias intensity
*        01/17/99 - Added a variable "osc_crrctn" that controls the way that
*                   duration is calculated in computing response spectra:
*                   osc_crrctn = 0: original (Boore, 1983); no correction; 
*                                   no longer used (supplanted by 
*                                   osc_crrctn = 1), included here for
*                                   completeness
*                   osc_crrctn = 1: Boore and Joyner (1984); used up to 
*                                   now in smsim
*                   osc_crrctn = 2: Liu&Pezeshk (1996)'s empirical model in
*                                   which L. J. Liu by replaced alpha with 
*                                   the spectral shape factor k and n = 2 
*                                   instead n= 3 in Boore and Joyner (1984)
*        03/13/99 - Rearranged "if then" statements for osc_crrctn
*                   such that osc_crrctn appears in the order 0, 1, 2
*        02/05/00 - Used Chuck Mueller's suggestion to change code in 
*                   cl68_integrand to avoid possible numerical problems.

      character e_or_a*1
      real amom0,amom1,amom2,amom4,anz,ane
      include 'smsim.fi'



c
c compute moments and frequencies, bandwidth factors, rms.
c
C
      amom0=amom_rv(0)
      amom1=amom_rv(1)
      amom2=amom_rv(2)
      amom4=amom_rv(4)

c-----RBH
      ARG=ABS(1.0-AMOM1*AMOM1/(AMOM0*AMOM2))
      if(abs(arg).lt.1.0e-20) arg=1.0e-20
      deltay=sqrt(arg)                    ! see MAIN, used in V-M or T&U calcs
      xi = amom2/sqrt(amom0*amom4)
      arg=(1.-xi*xi)
      if(abs(arg).lt.1.0e-20) arg=1.0e-20
      eps_rv=sqrt(arg)
      freq20=sqrt(amom2/amom0)/(2.0*pi)
      freq42=sqrt(amom4/amom2)/(2.0*pi)

c because the acceleration is a transient, the durations for
c computing rms and for determining n may be different.
c durex Duration of excitation) is used for computing N. 
c The rms is computed using
c trms, which is determined as durex for
c regular time series and durex + tosc * ( rf/(rf+avib) ) for the
c oscillator output, where tosc is the time for an oscillator
c to decay to 1/e and rf = (fosc * durex)**3. avib is an adjustable
c parameter, currently set to 1/3. The reason for using the oscillator
c duration is that the random vibration theory gives the ratio
c of peak to rms. If the spectral energy is spread over a number of
c cycles as in an oscillator, the "local" rms is smaller than
c if contained just within durex. 
c In effect, we are forced to apply a fixup in an
c attempt to get around this basic limitation.

c determine duration of excitation.
c----RBH
c	GET DURATION FOR THE PATH
c-----
	CALL GETDUR(PATHDUR,R)
	durex = dursource(w_fa, fa, w_fb, fb)+PATHDUR

c-----
C	determine duration of rms
c-----
	trms = durex

	if(iaorins.eq.2)then
		if(osc_crrctn.eq.0)then
c-----
C	no correction
c-----
			avib = 1.
			rf = 0.0
		else if (osc_crrctn .eq. 1) then
c-----
c			correction for Harmonic oscillator
c-----
c			use BJ correction
c-----
			avib = 1./3.
			rf = (fosc * durex)**3
		else
c-----
c			avib, rf are modified by Liu & Pezeshk 
c			(in paper submitted to BSSA in 1999) as follows:
c-----
			avib = deltay * sqrt(twopi)
			rf = (fosc * durex)**2
		endif
		tosc = 1.0/(twopi*damp*fosc)
		trms = trms + tosc * ( rf/(rf+avib) )
	else if(iaorins.eq.3)then
c-----
c		correction for narrow Butterworth filtered velocity
c-----
C		trms = trms + 2.0 * perosc
		trms = trms + 1.0 * perosc
	endif

	rms=sqrt(amom0/trms)

	ane = 2.0*freq42 * durex
	anz = 2.0*freq20 * durex
	if (ane .le. 1.0) ane = 1.002
	if (anz .le. 1.33) anz = 1.33

c-----
c	Compute Arias intensity (only meaningful if pga: iaorins = 1 and idva = 2)
c-----
	if (iaorins .eq. 1 .and. idva .eq. 2) then
c-----
c		 acceleration of gravity, assuming acc units of cm/s^2
c-----
		g = 980.0    
		arias_fctr = pi/(2.0*g)
		arias_rv = arias_fctr * amom0
	else
		arias_fctr = 0.0
		arias_rv = arias_fctr * amom0
	endif

c-----
c	factor of 2.0 because consider positive maxima
c	and negative minima.
c	ane is an estimate of the total number of extrema.
c	anz ( also=ane*sqrt(1.0-eps_rv*eps_rv) ) is an estimate of the
c	number of positive and negative zero crossings.
c
c
c	compute Cartwright and Longuet-Higgins estimates of peak/rms:
c
c-----
	pk_rms_cl_1 = sqrt(2.0*alog(anz))              ! 'cl' = Cart. & L-H
	pk_rms_cl_2 = pk_rms_cl_1+0.5772/pk_rms_cl_1
	pk_rms_cl_eq68 = cl68_numrcl_int( ane, xi, zup)
	e_or_a = 'e'                  ! might be used in a print statement
	pk_cl_eq68 = rms * pk_rms_cl_eq68    ! eq. 6.8 of C&L-H
	pk_cl_1 = rms * pk_rms_cl_1          ! 1 term asymptotic
	pk_cl_2 = rms * pk_rms_cl_2          ! 2 term asymptotic
	
	ETA05 = QINV(0.05,EPS_RV,ANE)
	ETA95 = QINV(0.95,EPS_RV,ANE)
	PK05_CL_EQ68 = RMS * ETA05    
	PK95_CL_EQ68 = RMS * ETA95    
	gmsim = pk_cl_eq68

c-----
c	that is all
c-----
      return
      end
*----------------- END GET_MOTION -----------------------------

*----------------- BEGIN CL68_NUMRCL_INT -----------------------------
      function cl68_numrcl_int( an_in, xi_in, zup)
* Numerical integration of eq. 6.8 in Cartwright and Longuet-Higgins
* Dates: 06/09/95 - Written by D.M. Boore, and tested using CHK_INT.
*                   I also plotted the integrand for typical values
*                   of an, xi, and found that it decays strongly to zero
*                   by a value of 5 for the variable.  I use 10 as an upper
*                   limit, which should be much more than enough.  The 
*                   integration routines are such, however, that I could 
*                   probably use a much larger m=number with little extra
*                   time.
*        12/19/95 - Added zup to exact_numrcl
*        12/28/95 - Name changed from exact_numrcl to cl68_numrcl
*        01/03/96 - Name changed from cl68_numrcl to cl68_numrcl_int
*        03/13/99 - On the advice of R. Herrmann, substituted qromb for
*                   qmidpnt.  Numerical tests indicate that both give the
*                   same answers, but qromb evaluates the function at the 
*                   endpoints (zlow and zup), whereas qmidpnt does not (and
*                   is appropriate for an improper integral for which
*                   the integrand cannot be evaluated right at the endpoints;
*                   this is not the case here).

      external cl68_integrand

      common /clint/ xi, an

      an = an_in
      xi = xi_in      
      zlow = 0.0 ! where is zup? It is read in as a parameter in the input file.
c     zup_step = 0.1
c     do while (cl68_integrand(zup) .eq. 0.0)
c      zup = zup - zup_step
c     end do
c     zup = zup + zup_step

      call qromb(cl68_integrand,zlow,zup,result)

      cl68_numrcl_int = result/sqrt(2.0)

      return
      end
*----------------- END CL68_NUMRCL_INT -----------------------------

*----------------- BEGIN CL68_INTEGRAND -----------------------------
      real function  cl68_integrand(z)
* Dates: 06/09/95 - Written by D.M. Boore.  See 7/11/82 notes for
*                   stochastic model, with 6/9/95 addition that uses
*                   a variable transformation to remove the sqrt
*                   singularity at the origin.
*        01/03/95 - Name changed from cl_int to cl68_integrand
*        02/05/00 - Made changes suggested by C. Mueller to avoid
*                   possible numerical problem when "an" is large

      common /clint/ xi, an
*      cl68_integrand = 2.0*(1.0-(1.0-xi*exp(-z**2))**an)  ! original
      y = an * alog(1.0-xi*exp(-z**2))   ! Mueller modification
      cl68_integrand = 2.0*(1.0-exp(y))  ! Mueller modification

      return
      end
*----------------- END CL68_INTEGRAND -----------------------------
      
*----------------- BEGIN AMOM_RV -----------------------------
      function amom_rv(i)
c Dates: 06/06/95 - Written by D.M. Boore, patterned after AMOMI, which 
*                   see for more detailed history.
*        11/14/95 - Obtain eps_int from get_params and pass through common

	COMMON/BUTFLT/FL,FH,NBUT,DDTT,DOBUT,IPOW
	REAL FL, FH, DDTT
	INTEGER NBUT,IPOW
	LOGICAL DOBUT

      external derivs

      include 'smsim.fi'

      h1 = 0.1
      hmin = 0.0
      imom = i       ! keep param i in parameter list rather than imom;
                     ! imom is passed through the rv common block

      result = 0.0
c-----
c	Use Numerical Recipes routine to integrate
c
c	odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs)
c		integrate from x1 to x2
c		Note for Butterworth filter the limits are changed
c-----
	if(DOBUT)THEN
		call odeint(result, 1, FH/10.0, FL*10.0, eps_int, h1, hmin,
     1            nok, nbad, derivs)
	ELSE
		call odeint(result, 1, 0.0, fup, eps_int, h1, hmin,
     1            nok, nbad, derivs)
	ENDIF


      amom_rv = 2.0 * result

      return
      end
*----------------- END AMOM_RV -----------------------------

*----------------- BEGIN DERIVS -----------------------------
* Dates: 06/07/95 - Written by D.M. Boore
      subroutine derivs(freq, y, dydf)

      include 'smsim.fi'

      f = freq
      if (freq .eq. 0.0) f = 0.001
      w = twopi * f
      
      a = spect_amp(f)

      if(imom .eq. 0) then
        dydf=a*a
      else
        dydf=a*a*w**imom
      endif

      return
      end
*----------------- END DERIVS -----------------------------


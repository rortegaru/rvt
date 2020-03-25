c-----
c	this file contains the various filter responses
c	that describe the measured motion. There are
c	two versions: time domain filters and complex frequency domain filters
c	
c	zpass	- time domain butterworth filter. 8'th order high pass a
c		fc/sqrt(2), 8'th order lowpass at fc*sqrt(2)
c	        subroutine zpass(x,nbut,ift,zphase,fl,fh,dt,n)
c			x     input time series
c			nbut  filter order
c			ift   0 low pass, 1 high pass, 2 band pass
c			zphase - .true. if zero phase
c			fl    lower corner of bandpass
c			fh    upper corner of bandpass
c	zpassf  - frequency domain filter
c	        subroutine zpassf(nbut,ift,zphase,fl,fh,z,zamp,dt,NPTS)
c			nbut  filter order
c			ift   0 low pass, 1 high pass, 2 band pass
c			zphase - .true. if zero phase
c			fl    lower corner of bandpass
c			fh    upper corner of bandpass
c			z	Cplx	- 0 + i omega
c			zamp	Cplx	- complex filter response
c			dt	R	- sampling interval
c			npts	I	- number of points -> get nyquist
c	psrv    - time domain response spectra filter
c		subroutine psrv(fc,damp,psv,x,a,sd,sv,sa,dt,n)
c			fc	R*4	- oscillator filter prequency
c			damp	R*4	- damping
c			psv	R*4	- Pseudo velocity spectrum
c			x	R*4	- ground acceleration
c			a	R*4	- output of oscillator
c			sd	R*4	- maximum displacement of oscillator
c			sv	R*4	- maximum velocity of oscillator
c			sa	R*4	- maximum accleration of oscillator
c			dt	R*4	- sample interval
c			n	I*4	- number of samples
c-----
c	psrvf   - frequency domain response spectra filter
c		subroutine psrv(fc,damp,z,zamp
c			fc	R*4	- oscillator filter prequency
c			damp	R*4	- damping
c			z	Cplx	- 0 + i omega
c			zamp	Cplx	- complex filter response
c-----
        subroutine zpass(x,nbut,ift,zphase,fl,fh,dt,n)
c-----
c	x	R*4	- input time series
c	nbut	I*4	-
c	ift	I*4	- 0 low pass		^
c			  1 high pass		^
c			  2 band pass
c	zphase	L	- zero phase if .true.
c	fl	R*4	- low pass corner frequency
c	fh	R*4	- high pass corner frequency measured down from
c				nyquist
c-----	
c	This procedure uses a Butterworth filter. The high pass is
c	based on a low pass design.
c
c	Low Pass  H(s)  = 1/ [ (s/A) + 1 ]   
c
c	Bilinear Transform
c
c	A ->   (2/T) tan (AT/2), where T is the sampling interval
c
c	s->  (2/T) [ 1 - z^-1]/[1 + z^-1]
c
c	Thus (s/A) -> (1/a)[ 1 - z^-1]/[1 + z^-1] where a = tan(AT/2)
c
c	H(z) = a( 1 + z^-1) / [ (1+a) + (a-1)z^-1 ]
c	The low pass recursive filter is thus
c
c	y(n) = { a [ x(n) + x(n-1) ] - (a-1) y(n-1) } / (1 + a)
c
c	High pass H(s) = 1/ [ (A/s) + 1 ]
c	H(z) = ( 1 - z^-1) / [ (1+a) + (a-1)z^-1 ]
c	     = (1/a)( 1 - z^-1) / [ (1+1/a) - (1-1/a)z^-1 ]
c	Note that
c	1/a = cot(AT/2) = tan(pi - AT/2) = tan(A'T/2) where A' = wN - a,
c		and wN is the Nyquist
c
c	Thus a low pass equation can be made into a High pass by
c	1) changing the sign of the z^-1 term, and
c	2) using A' = Wn -a as the 'low pass corner'
c-----
        dimension x(n)
	logical zphase
        fnyq=1./(2.*dt)
        tupi=2.*3.141592654
	ict = 9999
c-----
c	iflag	0 - single pass filter, either low or high pass
c	iflag	1 - double pass to make bandpass
c-----
	if(ift.eq.0)then
		a=1.
		w0=tupi*fl
        	iflag=0
	else if(ift.eq.1)then
		a=-1.
		w0=tupi*(fnyq-fh)
        	iflag=0
	else if(ift.eq.2)then
		a=1.
		w0=tupi*fl
		w1=tupi*(fnyq-fh)
		iflag=1
	endif

  25    continue
        iq=(1+(-1)**nbut)/2
c-----
c	iq = 1 if nbut is even
c	iq = 0 if nbut is odd
c-----
        n2b=nbut/2
c-----
c	run through the second order butterworths
c-----
	do 35 i=1,n2b
		call cs2(ct,w0,nbut,iq,i)
		ict=i
		call bttr2(x,a,ct,w0,ict,n,dt)
  35	continue
	if(iq.eq.0)then
		call bttr1(x,a,w0,n,dt)
	endif
	if(zphase)then
		call revers(x,n)
		do 45 i=1,n2b
			call cs2(ct,w0,nbut,iq,i)
			ict=i
			call bttr2(x,a,ct,w0,ict,n,dt)
  45		continue
		if(iq.eq.0) then
			call bttr1(x,a,w0,n,dt)
		endif
		call revers(x,n)
	endif
        if(iflag.eq.0) go to 60
c-----
c	now perform a high pass filter
c-----
        a=-1.
        w0=w1
        iflag=0
        go to 25
  60    continue
        return
        end

        subroutine bttr2(x,a,b,w,ict,n,dt)
c-----
c	2 pole filter s^2 + s b + w^2
c
c	x	R*4	- time series to be filtered, returned in place
c	a	R*4	- 1 for LP, -1 for HP peak response = 1
c	b	R*4	- filter coefficient
c	w	R*4	- filter coefficient
c	ict	I*4	- if HP, then later stages always start at 0,0
c	n	I*4	- number of points in time series
c-----
        dimension x(n)
	za = 2.0*a
c-----
c	frequency warping for bilinear - note that wp = warped(w0)*(2/T)
c-----
        call set(w,wp,ak,dt)
        wp2=wp*wp
        b0=wp2
        b1=1.+b*ak+wp2
        b2=za*(wp2-1.)
        b3=1.-b*ak+wp2
	if(a.le.0.0)then
        	yt1=0.0
        	yt2=0.0
        	if(ict.le.1) then
        		xt1=x(1)
        		xt2=x(1)
		else
			xt1=0.0
        		xt2=0.0
		endif
	else
		yt1=x(1)
		yt2=x(1)
		xt1=x(1)
		xt2=x(1)
	endif
c-----
c	yt3 = y(n)
c	yt2 = y(n-1)
c	yt1 = y(n-2)
c	xt2 = x(n-1)
c	xt1 = x(n-2)
c-----
	do 10 i=1,n
		yt3=(b0*(x(i)+za*xt2+xt1)-b2*yt2-b3*yt1)/b1
		xt1=xt2
		xt2=x(i)
		yt1=yt2
		yt2=yt3
c-----
c	save storage by using original array for output
c-----
		x(i)=yt3
  10	continue
	return
	end

        subroutine bttr1(x,a,w0,n,dt)
c-----
c	single pole filter
c-----
c	x	R*4	- time series to be filtered, output is
c				stored in the same array
c	a	R*4	-  1 low pass
c			  -1 high pass
c	w0	R*4	- angular corner frequency of filter
c-----
c	Lowpass				Highpass
c
c	  w0				  s
c	_____				_____
c	(s+w0)				(s+w0)
c
c	w0(1 + z^-1)			(2/T)(1 - z^-1)
c	-----------------------		---------------------------
c	(2/T +w0) - z^-1(2/T -w0) 	(2/T +w0) - z^-1(2/T -w0)
c-----
        dimension x(n)
c-----
c	frequency warping for bilinear - note that wp = warped(w0)*(2/T)
c-----
        call set(w0,wp,ak,dt)
c-----
c	(2/T)b0 == (2/T + w0)	(2/T)b1 == (sign) ( a - 2/T) 
c-----
        b0=wp+1.0
        b1=a*(wp-1.0)
c-----
c	If high pass then do not worry about initial DC
c	If low pass, then preserve the DC offset
c-----
	xt = x(1)
	if(a.le.0.0)then
		yt = 0.0
	else
		yt = x(1)
	endif
	do 10 i=1,n
		yt=(wp*(x(i)+a*xt)-b1*yt)/b0
		xt=x(i)
		x(i)=yt
   10	continue
	return
	end

	subroutine cs2(ct,w0,nbut,iq,i)
c-----
c	define the real value of the complex pole
c	corresponding to the butterworth filter of order nbut
c
c	For an n'th order Butterworth filter there are
c	n poles in the s < 0 complex plane. The poles
c	are symmetric about the negative real s-axis
c
c	The angular separation of poles is pi/n 
c	If n is odd, one pole is on the negative axis
c	if n is even, the first is 0.5 the separation away
c
c	This is part of the second order filter response, e.g., a
c	filter of the form
c
c	s^2 + 2 s cos(theta) w0 + w0^2
c
c	or since the poles are in complex conjugate pairs, the
c	filter denominator is  (s + a)(s + a*) = s^2 + 2s Re(a) + aa*
c
c-----
c	ct	R*4	- $
c	w0	R*4	- 2 w0 cos(theta
c	i	I*4	- the number of the pole
c	iq	I*4	- 0 number of poles is even
c			  1 number of poles is odd
c-----
		tupi = 6.2831853
        	ct=2.*w0*cos((2*i-iq)*tupi/(4.0*nbut))
        return
        end

        subroutine set(w,wp,ak,dt)
c-----
c	perform frequency warping required for bilinear filter
c-----
c	w	R*4	- omega
c	wp	R*4	- factor such that OMEGA = (2/dt)*wp
c	ak	R*4	-
c-----
        wp=tan(w*dt/2.)
        ak=wp/w
        return
        end

	subroutine revers(x,n)
c-----
c	reverse the time series x
c
c	x	R*4	- time series of n elements
c	n	I*4	- length of the time series
c-----
	dimension x(n)
c-----
c	define a temporary buffer since the initial array may be very large
c	this will make the exchange more efficient
c-----
	parameter (NBUF=500)
	dimension buf(NBUF)
	m=NBUF
c-----
c	simple case
c-----
	if(n.eq.1)return
c-----
c	the general case
c-----
	nn=n
	k=n/m
	mr=n-k*m
c	mr = mod(k,n)
	jj=1
	do 40 i=1,k+1
		if(i.eq.k+1) m=mr
		if(m.eq.0) go to 40
		do 10 j=1,m
			buf(j)=x(n-m+j)
   10		continue
		if(n.ne.m)then
			nn=nn-m
			if(nn.ne.0) then
				do 20 j=1,nn
					x(n-j+1)=x(n-m-j+1)
   20				continue
			endif
		endif
		do 30 j=1,m
			x(jj)=buf(m-j+1)
			jj=jj+1
   30		continue
   40	continue
	return
	end


        subroutine zpassf(nbut,ift,zphase,fl,fh,z,zamp,dt,n)
c-----
c	nbut	I*4	-
c	ift	I*4	- 0 low pass		^
c			  1 high pass		^
c			  2 band pass
c	zphase	L	- zero f=phase is .true.
c	fl	R*4	- low pass corner frequency
c	fh	R*4	- high pass corner frequency measured down from
c				nyquist
c	z	Cplx	- 0 + i omega
c	zamp	Cplx	- complex filter response
c-----	
c	This procedure uses a Butterworth filter. The high pass is
c	based on a low pass design.
c
c	Low Pass  H(s)  = 1/ [ (s/A) + 1 ]   
c
c	Bilinear Transform
c
c	A ->   (2/T) tan (AT/2), where T is the sampling interval
c
c	s->  (2/T) [ 1 - z^-1]/[1 + z^-1]
c
c	Thus (s/A) -> (1/a)[ 1 - z^-1]/[1 + z^-1] where a = tan(AT/2)
c
c	H(z) = a( 1 + z^-1) / [ (1+a) + (a-1)z^-1 ]
c	The low pass recursive filter is thus
c
c	y(n) = { a [ x(n) + x(n-1) ] - (a-1) y(n-1) } / (1 + a)
c
c	High pass H(s) = 1/ [ (A/s) + 1 ]
c	H(z) = ( 1 - z^-1) / [ (1+a) + (a-1)z^-1 ]
c	     = (1/a)( 1 - z^-1) / [ (1+1/a) - (1-1/a)z^-1 ]
c	Note that
c	1/a = cot(AT/2) = tan(pi - AT/2) = tan(A'T/2) where A' = wN - a,
c		and wN is the Nyquist
c
c	Thus a low pass equation can be made into a High pass by
c	1) changing the sign of the z^-1 term, and
c	2) using A' = Wn -a as the 'low pass corner'
c-----
	logical zphase
	complex z, zamp
        fnyq=1./(2.*dt)
        tupi=2.*3.141592654
	ict = 9999
c-----
c	iflag	0 - single pass filter, either low or high pass
c	iflag	1 - double pass to make bandpass
c-----
	if(ift.eq.0)then
		a=1.
		w0=tupi*fl
        	iflag=0
	else if(ift.eq.1)then
		a=-1.
		w0=tupi*(fnyq-fh)
        	iflag=0
	else if(ift.eq.2)then
		a=1.
		w0=tupi*fl
		w1=tupi*(fnyq-fh)
		iflag=1
	endif

  25    continue
        iq=(1+(-1)**nbut)/2
c-----
c	iq = 1 if nbut is even
c	iq = 0 if nbut is odd
c-----
        n2b=nbut/2
c-----
c	run through the second order butterworths
c-----
	do 35 i=1,n2b
		call cs2(ct,w0,nbut,iq,i)
		ict=i
		call bttr2f(a,ct,w0,ict,n,dt,z,zamp)
  35	continue
	if(iq.eq.0)then
		call bttr1f(a,w0,n,dt,z,zamp)
	endif
	if(zphase)then
		do 45 i=1,n2b
			call cs2(ct,w0,nbut,iq,i)
			ict=i
			call bttr2f(a,ct,w0,ict,n,dt,conjg(z),zamp)
  45		continue
		if(iq.eq.0) then
			call bttr1f(a,w0,n,dt,conjg(z),zamp)
		endif
	endif
        if(iflag.eq.0) go to 60
c-----
c	now perform a high pass filter
c-----
        a=-1.
        w0=w1
        iflag=0
        go to 25
  60    continue
        return
        end

        subroutine bttr2f(a,b,w,ict,n,dt,z,zamp)
c-----
c	2 pole filter s^2 + s b + w^2
c
c	x	R*4	- time series to be filtered, returned in place
c	a	R*4	- 1 for LP, -1 for HP peak response = 1
c	b	R*4	- filter coefficient
c	w	R*4	- filter coefficient
c	ict	I*4	- if HP, then later stages always start at 0,0
c	n	I*4	- number of points in time series
c-----
	complex z, zamp
	za = 2.0*a
c-----
c	frequency warping for bilinear - note that wp = warped(w0)*(2/T)
c-----
        call set(w,wp,ak,dt)
        wp2=wp*wp
        b0=wp2
        b1=1.+b*ak+wp2
        b2=za*(wp2-1.)
        b3=1.-b*ak+wp2
	zamp = zamp * b0 * ( 1.0 + za*conjg(z)+conjg(z)*conjg(z))/
     1		( b1 + b2*conjg(z) + b3*conjg(z)*conjg(z))
	return
	end

        subroutine bttr1f(a,w0,n,dt,z,zamp)
c-----
c	single pole filter
c-----
c	x	R*4	- time series to be filtered, output is
c				stored in the same array
c	a	R*4	-  1 low pass
c			  -1 high pass
c	w0	R*4	- angular corner frequency of filter
c-----
c	Lowpass				Highpass
c
c	  w0				  s
c	_____				_____
c	(s+w0)				(s+w0)
c
c	w0(1 + z^-1)			(2/T)(1 - z^-1)
c	-----------------------		---------------------------
c	(2/T +w0) - z^-1(2/T -w0) 	(2/T +w0) - z^-1(2/T -w0)
c-----
	complex z, zamp
c-----
c	frequency warping for bilinear - note that wp = warped(w0)*(2/T)
c-----
        call set(w0,wp,ak,dt)
c-----
c	(2/T)b0 == (2/T + w0)	(2/T)b1 == (sign) ( a - 2/T) 
c-----
        b0=wp+1.0
        b1=a*(wp-1.0)
c-----
c	If high pass then do not worry about initial DC
c	If low pass, then preserve the DC offset
c-----
	zamp = zamp * wp*(1.0 + a*conjg(z))/(b0 + b1*conjg(z))
	return
	end


	subroutine psrv(fc,damp,psv,x,a,sd,sv,sa,dt,n)
c-----
c	fc	R*4	- oscillator filter prequency
c	damp	R*4	- damping
c	psv	R*4	- Pseudo velocity spectrum
c	x	R*4	- ground acceleration
c	a	R*4	- output of oscillator
c	sd	R*4	- maximum displacement of oscillator
c	sv	R*4	- maximum velocity of oscillator
c	sa	R*4	- maximum accleration of oscillator
c	dt	R*4	- sample interval
c	n	I*4	- number of samples
c-----
	dimension x(n),a(n)
	real*4 psv,sd,sv,sa
	real*8 t0, tp, ts, tfirst
c-----
c		Damping = 0.05 e.g., 5%
c-----
		tupi=2.0*3.1415927
		dmp=damp
c-----
c		baseline correction?
c		call pcn01(x,v,t)
c-----
		call pcn03(n,dt,dmp,fc,a,x,1.,sd,sv,sa,fs,td,tv)
c-----
c		get time of maxima in terms of time after origin time
c-----
		w=tupi*fc
		psv=w*sd
	return
	end

c ***********************************************
c
c       subroutine for calculation of spectra from earthquake record
c               digitized at equal time intervals
c                       pcn03
c       Nigam and Jennings, Caltech publication, Pasadena, CA, June 1968
c
c ************************************************
	subroutine pcn03(n,del,dmp,fc,oa,ga,sf,sd,sv,sa,fs,td,tv)
c-----
c	n	I*4	- number of points in time series
c	del	R*4	- data sampling interval in seconds
c	dmp	R*4	- array of damping, e.g., 0.05 for 5%
c	fc	R*4	- array of filter frequencies
c	oa	R*4	- array of oscillator output
c	ga	R*4	- array of accelerations
c	sf	R*4	- scale factor
c	sd	R*4	- oscillator displacment
c	sv	R*4	- oscillator velocity
c	sa	R*4	- oscillator acceleration
c	fs	R*4	- Fourier amplitude spectra which is
c				last time value
c	td	R*4	- time of maximum displacment
c	tv	R*4	- time of maximum velocity
c-----
	real*4 ga(n), oa(n)
	dimension x(3),g(2),a(2,2),b(2,2),ty(3)
c-----
c	run input through Single Degree of Freedom oscillator,
c	but keep processing to see if maximum occurs after the
c	end of the time series because of the impulse response
c	of the filter
c-----
c	loop over damping coefficient
c-----
		d=dmp
c-----
c		loop over oscillator frquency
c-----
			p=1./fc
			w=2.*3.141592654/p
c-----
c			choice of interval of integration , which is
c			smaller than original delta t
c-----
			delp=p/20.
			l=del/delp+1.-1.e-05
			delt=del/float(l)
c-----
c			computaion of matrices a and b 
c				for this frequency and damping
c-----
			call pcn04(d,w,delt,a,b)
c-----
c			initialization
c-----
			x(1)=0.0
			x(2)=0.0
			x(3)=0.0
			dmax=0.0
			vmax=0.0
			amax=0.0
			i=1
			dw=2.*w*d
			w2=w**2
			ia=2.*p/delt+1.e-05
c-----
c			computation of oscillator response
c-----
			gap1 = ga(2)
			gam1 = ga(1)
			oa(1) = 0.0
   7			continue
			sl=(gap1-gam1)/float(l)
			do 6 m=1,l
				g(1)=(gam1+sl*float(m-1))*sf
				g(2)=(gam1+sl*float(m))*sf
				ty(1)=a(1,1)*x(1)+a(1,2)*x(2)
     1					+b(1,1)*g(1)+b(1,2)*g(2)
				ty(2)=a(2,1)*x(1)+a(2,2)*x(2)
     1					+b(2,1)*g(1)+b(2,2)*g(2)
				ty(3)=-(dw*ty(2)+w2*ty(1))
c-----
c				monitoring the max. values 
c-----
				time=i*del
				if(abs(ty(1)).gt.abs(dmax)) then
					tdmax=time
					dmax=ty(1)
				endif
				x(1)=ty(1)
				if(abs(ty(2)).gt.abs(vmax)) then
					tvmax=time
					vmax=ty(2)
				endif
				x(2)=ty(2)
				if(abs(ty(3)).gt.abs(amax)) then
					tamax=time
					amax=ty(3)
				endif
				x(3)=ty(3)
    6			continue
c-----
c			test for end of integration
c-----
			i=i+1
			if(i.le.n)oa(i) = x(1)
c			if(abs(fc - 8.0).lt.0.001)write(9,*)i,ga(i),x(1)

			if(i.eq.n)then
				vend = x(2)
			endif
			gam1 = gap1
			if(i.lt.n)then
				gap1 = ga(i+1)
			else if(i.ge.n)then
				gap1 = gam1
			endif
			if(i.lt. (n+ia))go to 7
			if(d.le.1.0e-03) then
				fs=abs(vend)
			else
				fs = 0.0
			endif
			sd=abs(dmax)
			sv=abs(vmax)
			sa=abs(amax)
			td=tdmax
			tv=tvmax
			ta=tamax
	return
	end


c ******************************************************
c
c       subroutine for comutation of matrices a and b
c                       pcn04
c     Nigam and Jennings, Caltech publication, Pasadena, CA, June 1968
c
c ******************************************************
	subroutine pcn04(d,w,delt,a,b)
c-----
c	Define matrix elements for recursive solution of
c	
c	x'' + 2dwx' +wwx = a
c-----
c	( x  ) =  ( a11    a12 ) ( x  ) + ( b11  b12 )(a(i)   )
c	( x' ) =  ( a21    a22 ) ( x' ) + ( b21  b22 )(a(i+1) )
c	     i+1                       i
c-----	
c	d	R*4	- damping value
c	w	R*4	- oscillator natural frequency
c	delt	R*4	- sampling interval
c	a	R*4	- a matrix
c	b	R*4	- b matrix
c-----
	implicit double precision (a-h,o-z)
	real*4 a(2,2),b(2,2), d, w, delt
		dw = d*w
		w2 = w*w
		w3 = w*w*w
		wd = w*dsqrt(1.0d+00-d*d)
		a0=dexp(-dw*delt)
		if(wd*delt .lt. 0.00001)then
			swd = delt
		else
			swd = dsin(wd*delt)/wd
		endif
		cwd = dcos(wd*delt)
		wdswd = wd* dsin(wd*delt)

		a(1,1) = a0 *(dw*swd + cwd)
		a(1,2) = a0 * swd
		a(2,1) = - a0 * w2 * swd
		a(2,2) = a0 *(cwd - dw * swd)

		fac = (2*d*d - 1.0d+00)/(w2*delt)
		w2i = 1.0/w2
		dow = d / w
		dow3 = 2.0*d/(w3*delt)

		dum = a0 * ( (fac + dow)*swd + (dow3 + w2i)*cwd)
     1			- dow3
		b(1,1) = sngl(dum)
		diff = a0 * ( dow * swd + w2i * cwd) + w2i

		dum = - a0 * ( fac * swd + dow3*cwd)
     1			- w2i + dow3
		b(1,2) = sngl(dum)

		dum = a0 *( ( fac + dow)*( cwd - dw*swd)
     1				- (dow3 + w2i)*(wdswd + dw*cwd) )
     2			+ 1.0/(w2*delt)
		b(2,1) = sngl(dum)

		dum = - a0* ( fac * ( cwd - dw*swd)
     1				- dow3*(wdswd + dw*cwd) )
     2				- 1.0/(w2*delt)
		b(2,2) = sngl(dum)

	return
	end
c  **********************************
c      subroutine for parabolic base line correction
c                      pcn01
c    Nigam and Jennings, Caltech publication, Pasadena, CA, June 1968
c  **********************************
        subroutine pcn01(x,v,t)
        common /dtn/dt,n
        dimension x(n),v(n),t(n)
        v0=0.0
        be1=0.0
        be2=0.0
        be3=0.0
        t(1)=0.0
        v(1)=v0
        m=n-1
        dt2=dt*dt
        do 1 i=1,m
        t(i+1)=t(i)+dt
        t2=t(i)*t(i)
        t3=t2*t(i)
        t4=t(i)*t(i+1)
        t5=t(i+1)*t(i+1)
        t6=t(i+1)*t5
        v(i+1)=v(i)+0.5*dt*(x(i)+x(i+1))
        be1=be1+0.5*v(i)*dt*(t(i)+t(i+1))+0.04166667*(x(i)*(3.*t(i)
     1  +0.5*t(i+1))+x(i+1)*(t(i)+3.*t(i+1)))*dt2
        be2=be2+0.3333333*v(i)*dt*(t2+t4+t5)+0.01666667*dt2*(x(i)*
     1  (4.*t2+7.*t4+9.*t5)+x(i+1)*(t2+3.*t4+6.*t5))
        be3=be3+0.25*v(i)*dt*(t3+t2*t(i+1)+t(i)*t5+t6)+
     1  0.008333333*dt2*(x(i)*(0.5*t3+9.*t2*t(i+1)+12.*t(i)*t5+14.*t3)
     2  +x(i+1)*(t3+3.*t2*t(i+1)+6.*t(i)*t5+10.*t6))
   1    continue
        yo=t(n)
        beta1=be1/yo**3
        beta2=be2/yo**4
        beta3=be3/yo**5
        ca1=300.*beta1-900.*beta2+630.*beta3
        ca2=-1800.*beta1+5760.*beta2-4200.*beta3
        ca3=1890.*beta1-6300.*beta2+4725.*beta3
        c0=ca1
        c1=ca2/yo
        c2=ca3/yo**2
        do 5 i=1,n
   5    x(i)=x(i)-c0-c1*t(i)-c2*t(i)*t(i)
        return
        end

	subroutine psrvf(fn,damp,s,zamp)
c-----
c	single degree of freedom response to acceleration
c-----
c			fn	R*4	- oscillator filter prequency
c			damp	R*4	- damping
c			s	Cplx	- 0 + i omega
c			zamp	Cplx	- complex filter response
c-----
	real fn, damp
	complex s, zamp

	wn = 6.2831853*fn
	zamp = cmplx(-1.0,0.0)/(s*s + 2.0*damp*wn*s + wn*wn)
	return
	end

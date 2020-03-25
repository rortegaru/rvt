	subroutine rvsim(niter,kseed,  ddt, 
     1		dofilt, dofourier,fn,zeta,fmax,instyp,
     1		doftbar,peak,peak05,peak95,xmommg,
     1		sigma,vssrc,rhsrc,jsrc,dist,q0,eta,velq,kappa)
c-----
c	niter	I	- number of iterations
c	iseed	I	- random  number seed 
c	dt	R	- sampling interval
c	acc	R	- acceleration array
c	vel	R	- velocity array
c	dis	R	- displacment array
c	dofilt	
c	dofourier
c	fn
c	zeta
c	fmax
c	instyp
c	doftbar
c	peak
c	dist	R	- distance
c-----
	integer niter, kseed
	real  ddt
	logical dofilt, dofourier, doftbar
	real fn, zeta, fmax
	integer instyp
	real kappa


	

	complex  s, zamp
	real cfac
C-----
C	RBH COMMON
C-----
	COMMON/BUTFLT/FL,FH,NBUT,DDTT,DOBUT,IPOW
	REAL FL, FH, DDTT
	INTEGER NBUT,IPOW
	LOGICAL DOBUT
c-----
c	common to boore's routines
c		indxwind = 1  exponential
c		iran_type = 0 Gaussian
c		fcut = fmax
c-----
	include 'smsim.fi'

	iran_type = 0
c-----
c	common/magdist
c-----
	amag = xmommg
	r = dist
c-----
c	source_scale_params
c	common /source_scale_params/ stressc, dlsdm, amagc,  
c    :                       stress, fa, fb, am0, am0b_m0, fbdfa
c-----
	stressc = sigma
	dlsdm = 0.0
	fbdfa = 4.0
	amagc = 7.0
c-----
c      common /const_params/ rho, beta, prtitn, rtp, fs, pi, twopi, 
c     :                      const, iaorins, idva, freq_indep_factor
c-----
	beta = vssrc
	rho = rhsrc
	prtitn = 0.707
	rtp = 0.55
	fs = 2.0
	r_ref = 1.0
	iaorins = 1
	idva = 2
	pi = 3.1415927
	twopi = 6.2831853
c-----
c      common /site_dimin_params/ fm, akappa, dkappadmag, amagkref
c-----
	akappa = kappa
	dkappadmag = 0.0
	amagkref = 0.0
	fm = fmax
c-----
c      common /source_shape_params/ numsource, pf, pd
c-----
	numsource = jsrc
	pf = 2.0
	pd = 1.0
c-----
c	tdsim_params
c-----
	dt = ddt
	tshift = 5.0
	nsims = niter
	seed = kseed
	iseed = -abs(kseed)
	dur_fctr = 1.0
c-----
c	lowcut_filter__params
c-----
	fcut =  0.05
	norder = 2
c-----
c      common /q_params/ fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q
c-----
	fr1 = 1.0
	qr1 = q0
	s1 = eta
	ft1 = 1.0
	ft2 = 1.0
	fr2 = 1.0
	qr2 = q0
	s2 = eta
	c_q = velq
c-----
c      common /wind_params/indxwind,taper,twdtmotion,eps_wind,eta_wind
c-----
	indxwind = 1
	taper = 0.05
	twdtmotion = 1.0
	eps_wind = 0.2
	eta_wind = 0.05

	new_mr = .true.
c      common /rv/ fup, imom, zup, eps_int, amp_cutoff, osc_crrctn,
c     :            eps_rv, ane, anz,
c     :            pk_rms_cl_1, pk_rms_cl_2, 
c     :            pk_cl_eq68, pk_rms_cl_eq68,
c     :            pk_cl_1, pk_cl_2,
c     :            arias_rv
	zup = 10.0
	eps_int = 0.00001
	eps_int = 0.001
	amp_cutoff = 0.001
c-----
c	osc_crrctn = 1  Boore, D.M. and W. B. Joyner (1984). A note on the use
c		of random vibration theory to predict peak amplitudes of
c		transient signals, Bull. Seism. Soc. Am. 74, 2035-2039
c	osc_crrctn = 2  Liu, L. and S. Pezeshk (1999). An improvement on
c		the estimation of pseidoresponse spectra velocity using RVT 
c		method, Bull. Seism. Soc. Am. 89, 1384-1389.
c-----
	osc_crrctn = 2
c-----
c	common /rs/ perosc, fosc, damp
c-----
	fosc = fn
	damp = zeta
	if(fn.gt.0.0)then
		perosc = 1.0/fn
	endif
	
c-----
c	compute fup  is also done in rvtdsubs, but if the driver is changed c-----
	if( akappa .eq. 0.0) then
		fup = fm/amp_cutoff**0.25
	else  
		fup = 
     :   amin1(fm/amp_cutoff**0.25, -alog(amp_cutoff)/(pi*akappa))
	endif
	



c-----
c	If filter is used, set up the filter parameters
c-----
	dobut = dofilt
	if(dofilt)then
		nbut = 8
		fl = fn*1.414
		fh = fn/1.414
		ddtt = dt
c-----
c		get normalization at the center frequency fn
c-----		
		fac = 6.2831853*fn*dt
c-----
c		     i omega dt
c		z = e
c-----
		s = cmplx(cos(fac),sin(fac))
		zamp = cmplx(1.0,0.0)
		call zpassf(nbut,2,.false.,fl,fh,s,zamp,dt,1)
		cfac = cabs(zamp)
	endif

c-----
c	compute specific motions
c-----
	if(instyp.eq.1)then
c-----
c		acceleration
		idva = 2
c-----
		if(dofilt)then
			iaorins = 3
			call smrv(pga,pga05,pga95,dofilt)
			pga = pga / cfac
			pga05 = pga05 / cfac
			pga95 = pga95 / cfac
		else
			iaorins = 1
			call smrv(pga,pga05,pga95,dofilt)
		endif
		peak = pga
		peak05 = pga05
		peak95 = pga95
	else if(instyp.eq.2)then
		idva = 1
c-----
c		velocity
c-----
		if(dofilt)then
			iaorins = 3
			call smrv(pgv,pgv05,pgv95,dofilt)
			pgv = pgv / cfac
			pgv05 = pgv05 / cfac
			pgv95 = pgv95 / cfac
		else
			iaorins = 1
			call smrv(pgv,pgv05,pgv95,dofilt)
		endif
		peak = pgv
		peak05 = pgv05
		peak95 = pgv95
	else if(instyp.eq.3)then
c-----
c		displacement
c-----
		idva = 0
		if(dofilt)then
			iaorins = 3
			call smrv(pgd,pgd05,pgd95,dofilt)
			pgd = pgd / cfac
			pgd05 = pgd05 / cfac
			pgd95 = pgd95 / cfac
		else
			iaorins = 1
			call smrv(pgd,pgd05,pgd95,dofilt)
		endif
		peak = pgd
		peak05 = pgd05
		peak95 = pgd95
	else if(instyp.ge.9 .and.instyp.le.13)then
		iaorins = 2
		perosc = 1.0/fn
		fosc = fn
c-----
c		filter the acceleration and select 
c		sd, sv, sa,  psv or psa
c-----
		if(instyp.eq.9)then
c-----
c			SD
c-----
			idva = 0
			ipow = 0
			call smrv(psv,psv05,psv95,dofilt)
			peak = psv / (6.2831853*fosc)
			peak05 = psv05 / (6.2831853*fosc)
			peak95 = psv95 / (6.2831853*fosc)
		else if(instyp.eq.10)then
c-----
c			SV
c-----
			idva = 0
			ipow = 1
			call smrv(sv,sv05,sv95,dofilt)
			peak = sv
			peak05 = sv05
			peak95 = sv95
		else if(instyp.eq.11)then
c-----
c			SA
c-----
			idva = 0
			ipow = 2
			call smrv(sa,sa05,sa95,dofilt)
			peak = sa
			peak05 = sa05
			peak95 = sa95
		else if(instyp.eq.12)then
c-----
c			PSV
c-----
			idva = 0
			ipow = 0
			call smrv(psv,psv05,psv95,dofilt)
			peak = psv
			peak05 = psv05
			peak95 = psv95
		else if(instyp.eq.13)then
c-----
c			PSA
c-----
			idva = 0
			ipow = 0
			call smrv(psv,psv05,psv95,dofilt)
			peak = psv * (6.2831853*fosc)
			peak05 = psv05 * (6.2831853*fosc)
			peak95 = psv95 * (6.2831853*fosc)
		endif
	endif
C	WRITE(0,*)'PEAK,FREQ20,INSTYP,DOFILT,CFAC:',PEAK,FREQ20,INSTYP,DOFILT,CFAC
	
	return
	end

	subroutine smrv(vmax,vmax05,vmax95,dofilt)
c-----
c	vmax	R	peak value
c	dofilt	L	.true. use Butterworth
c-----
c	In common blocks
c	
c	iaorins	I	1 no instrument response
c			2 harmonic oscillator
c			3 not used
c	idva	I	2 acceleration
c			1 velocity
c			0 displacement
	real vmax
	logical dofilt
c-----
	include 'smsim.fi'
c-----
c	Begin RV -- this part follows Boore's SMSIM_RV.FOR
c-----


c-----
c	Set spectral parameters:
c-----
	call spect_scale()   

c-----
c	Get frequency-independent factor:
c-----
	call const_am0_gsprd()

c-----
c	For iaorins = 2,
C	Assume that instrument response is relative to ground
C	displacement.  In this case idva = 0, and to make sure
C	that this is so, I include the following statement:
c-----
	if (iaorins.eq.3)then
		osc_crrctn = 0
	endif

	call get_motion(vmax)  
	vmax05 = PK05_CL_EQ68
	vmax95 = PK95_CL_EQ68
	return
	end


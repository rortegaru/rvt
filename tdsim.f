	subroutine tdsim(niter,kseed,  ddt, 
     1		dofilt, dofourier,fn,zeta,fmax,instyp,
     1		doftbar,avgmx,xmommg,
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
c	damp
c	fmax
c	instyp
c	doftbar
c	avgmx
c	dist	R	- distance
c-----
	integer niter, kseed
	real  ddt
	logical dofilt, dofourier, doftbar
	real fn, damp, fmax
	integer instyp
	real kappa


	
	integer N2
	parameter (N2=131072)
	real acc(N2), vel(N2), dis(N2)
	real out(N2)

	logical rmv_trend
	complex  s, zamp
	real cfac
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
c	common rs
c-----
	damp = zeta
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



c-----
c	If filter is used, set up the filter parameters
c-----
	if(dofilt)then
		nbut = 8
		fl = fn*1.414
		fh = fn/1.414
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
c	initialize sums
c-----
	avgprv = 0.0
	prvsumsq = 0.0
	prvsimsd = 0.0
	avgpga = 0.
	avgpgv = 0.
	avgpgd = 0.
	pgacumsq = 0.
	pgvcumsq = 0.
	pgdcumsq = 0.
	pgasimsd = 0.
	pgvsimsd = 0.
	pgdsimsd = 0.

	avgmx = 0.0
	avgmxsd = 0.0
	avgmxsq = 0.0

	rmv_trend = .false.

	do 1000 isum = 1, nsims
c-----
c		create the acceleration
c-----
		call getacc(acc)  
		v0 = 0.0
		d0 = 0.0
c-----
c		get ground velocity, displacement
c-----
		 call acc2vd(acc, npts, dt, rmv_trend, v0, d0, vel, dis)  
c-----
c		get pga, pgv, pgd
c-----
C	do 1010 j=1,npts
C		t = (j-1)*dt
C	WRITE(0,'(1X,F11.6,1P3(1X,E11.4))')T,ACC(J),VEL(J),DIS(J)
 1010	continue
		call domxmn(acc, npts, depmax, depmin, depmen)
		pga = amax1( abs(depmax), abs(depmin))

		call domxmn(vel, npts, depmax, depmin, depmen)
		pgv = amax1( abs(depmax), abs(depmin))


		call domxmn(dis, npts, depmax, depmin, depmen)
		pgd = amax1( abs(depmax), abs(depmin))

		avgpga = avgpga + pga
		avgpgv = avgpgv + pgv
		avgpgd = avgpgd + pgd
		pgacumsq = pgacumsq + pga*pga
		pgvcumsq = pgvcumsq + pgv*pgv
		pgdcumsq = pgdcumsq + pgd*pgd                     
c-----
c		compute specific motions
c-----
		if(instyp.eq.1)then
			if ( dofilt)then
c-----
c			get filtered acceleration
c-----
			call zpass(acc,nbut,2,.false.,fl,fh,dt,npts)
			call domxmn(acc, npts, depmax, depmin, depmen)
			vmax = amax1( abs(depmax), abs(depmin))/cfac
			else
				vmax = pga
			endif
			avgmx = avgmx + vmax
			avgmxsq = avgmxsq + vmax*vmax
		else if(instyp.eq.2)then
			if( dofilt)then
c-----
c			get filtered velocity
c-----
			call zpass(vel,nbut,2,.false.,fl,fh,dt,npts)
			call domxmn(vel, npts, depmax, depmin, depmen)
			vmax = amax1( abs(depmax), abs(depmin))/cfac
			else
				vmax = pgv
			endif
			avgmx = avgmx + vmax
			avgmxsq = avgmxsq + vmax*vmax
		else if(instyp.eq.3)then
			if(dofilt)then
c-----
c			get filtered displacement
c-----
			call zpass(dis,nbut,2,.false.,fl,fh,dt,npts)
			call domxmn(dis, npts, depmax, depmin, depmen)
			vmax = amax1( abs(depmax), abs(depmin))/cfac
			else
				vmax = pgd
			endif
			avgmx = avgmx + vmax
			avgmxsq = avgmxsq + vmax*vmax
		else if(instyp.ge.9 .and.instyp.le.13)then
c-----
c			filter the acceleration and select 
c			sd, sv, sa,  psv or psa
c-----
			call psrv(fn,damp,psv,acc,out,sd,sv,sa,dt,npts) 
			if(instyp.eq.9)then
				vmax = sd
			else if(instyp.eq.10)then
				vmax = sv
			else if(instyp.eq.11)then
				vmax = sa
			else if(instyp.eq.12)then
				vmax = psv
			else if(instyp.eq.13)then
				vmax = 6.2831853*fn*psv
			endif
			avgmx = avgmx + vmax
			avgmxsq = avgmxsq + vmax*vmax
		endif
 1000	continue
c-----
c	compute the averages
c-----
	pgasim = avgpga/float(nsims)
	pgvsim = avgpgv/float(nsims)
	pgdsim = avgpgd/float(nsims)        
	avgmx  = avgmx /float(nsims)
c-----
c	approximate standard deviations -- to preclude problens with 
c	division by (N-1) use just N as a guide
c-----
	pgasimsd = sqrt (( pgacumsq - nsims*pgasim*pgasim)/nsims)
	pgvsimsd = sqrt (( pgvcumsq - nsims*pgvsim*pgvsim)/nsims)
	pgdsimsd = sqrt (( pgdcumsq - nsims*pgdsim*pgdsim)/nsims)
	avgvmxsd = sqrt (( avgmxsq  - nsims*avgmx*avgmx  )/nsims)
c-----
c	compute specific motions
c-----
C      write(0,'(t5,a, t16,a, t25,a, t38,a, t46,a, t60,a)')
C     :    'pgd(cm)',    'std/pgd',
C     :    'pgv(cm/s)',  'std/pgv',
C     :    'pga(cm/s2)', 'std/pga'
C      write(0,'(1p,6(1x,e10.2))')
C     :     pgdsim, pgdsimsd/pgdsim,
C     :     pgvsim, pgvsimsd/pgvsim,
C     :     pgasim, pgasimsd/pgasim             
C	per=1.0/fn
C	write(0,
C     :      '( t4,f7.3, t13,f8.3, 1p, t22,e9.2, t32,e10.2,
C     :         t43,e10.2, 1x,e10.3)')   
C     :	per,fn,avgmx,avgvmxsd/avgmx
	
	
	return
	end


	subroutine domxmn(x,npts,depmax,depmin,depmen)
c-----
c	get extremal values of the time series
c-----
	real*4 x(*)
	real*4 depmax,depmin,depmen
	integer*4 npts
	depmax = -1.0e+38
	depmin =  1.0e+38
	sum = 0.0
	do 1000 i=1, npts
		if( x(i) .gt. depmax) depmax = x(i)
		if( x(i) .lt. depmin) depmin = x(i)
		sum = sum + x(i)
 1000	continue
	if(npts.gt.0)then
		depmen = sum / npts
	else
		depmen = -12345.
	endif
	return
	end

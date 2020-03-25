	subroutine getinfo(fname)
c-----
c	read the control file to define everything about the source
c	propagation model
c	This is done by cycling through key words
c-----
c	subroutine arguments
c-----
c	fname	- Ch*80		name of the control file
c-----
	character fname*(*)
c-----
c	common blocks
c-----
	parameter(NDMAX=100)
	common/propq/q0,eta,vel
		real q0,eta,vel
	common/propg/ndist, rr(NDMAX),rpow(NDMAX)
		integer ndist
		real rr, rpow
	common/propd/ntime, rt(NDMAX),  tt(NDMAX)
		integer ntime	
		real rt, tt
	common/props/nsite, fs(NDMAX),  ss(NDMAX), kappa
		integer nsite
		real fs, ss, kappa
	common/srcfm/fmax, jsrc, sigma, betas, dens
		real fmax, sigma, betas, dens
		integer jsrc
	common/proc/comment
		character comment*80
c-----
c	local program variables
c-----
	logical iext
	character cstr*40
	integer i
c-----
c	initialize everything
c-----
	call initp()
c-----
c	determine if the control file exists
c-----
	if(fname .eq. ' ')then
		iext = .false.
	else
		inquire(file=fname,exist=iext)
	endif
	if(.not. iext )return
c-----
c	now read the control file and act on key words
c-----
	open(1,file=fname,status='old',
     1		form='formatted',access='sequential')
	rewind 1
c-----
c	now read the file and act
c-----
 1000	continue
	read(1,'(a)',end=9999)cstr
		if     (cstr(1:4).eq.'KAPP' .or. cstr(1:4).eq.'kapp')then
			read(1,*,end=9999,err=9999)kappa
		else if(cstr(1:4).eq.'QETA' .or. cstr(1:4).eq.'qeta')then
			read(1,*,end=9999,err=9999)q0, eta
		else if(cstr(1:4).eq.'FMAX' .or. cstr(1:4).eq.'fmax')then
			read(1,*,end=9999,err=9999)fmax
		else if(cstr(1:4).eq.'QVEL' .or. cstr(1:4).eq.'qvel')then
			read(1,*,end=9999,err=9999)vel
		else if(cstr(1:4).eq.'COMM' .or. cstr(1:4).eq.'comm')then
			comment = cstr
		else if(cstr(1:4).eq.'DIST' .or. cstr(1:4).eq.'dist')then
			read(1,*,end=9999,err=9999)ndist
			do 2001 i=1,ndist
				read(1,*,end=9999,err=9999)rr(i), rpow(i)
 2001			continue
			rr(1) = 1.0
		else if(cstr(1:4).eq.'DURA' .or. cstr(1:4).eq.'dura')then
			read(1,*,end=9999,err=9999)ntime
			do 2002 i=1,ntime
				read(1,*,end=9999,err=9999)rt(i), tt(i)
 2002			continue
			if(ntime .eq. 1)then
				ntime = 2
				rt(2) = 1000.0
				tt(2) = tt(1)
			endif
		else if(cstr(1:4).eq.'SITE' .or. cstr(1:4).eq.'site')then
			read(1,*,end=9999,err=9999)nsite
			do 2003 i=1,nsite
				read(1,*,end=9999,err=9999)fs(i), ss(i)
 2003			continue
			if(nsite.eq.1)then
				nsite = 2
				fs(2) = 1000.0
				ss(2) = ss(1)
			endif
		else if(cstr(1:4).eq.'SHEA' .or. cstr(1:4).eq.'shea')then
				read(1,*,end=9999)betas
		else if(cstr(1:4).eq.'DENS' .or. cstr(1:4).eq.'dens')then
				read(1,*,end=9999)dens
		else if(cstr(1:4).eq.'SIGM' .or. cstr(1:4).eq.'sigm')then
				read(1,*,end=9999)sigma
				jsrc = 1
		else if(cstr(1:4).eq.'AT93' .or. cstr(1:4).eq.'at93')then
				jsrc = 3
		else if(cstr(1:6).eq.'AB98CA' .or. cstr(1:6).eq.'ab98ca')then
				jsrc = 6
		else if(cstr(1:6).eq.'AS2000' .or. cstr(1:6).eq.'as2000')then
				jsrc = 9
		endif
	go to 1000
 9999	continue
	close (1)
	return
	end


	subroutine getdur(dur,dist)
c-----
c	Obtain the duration as a function of distance
c
c	dur	Real	- propagation duration in seconds (rreturned)
c	dist	Real	- distance in km
c-----
	real dur, dist
	parameter(NDMAX=100)
	common/propd/ntime, rt(NDMAX),  tt(NDMAX)
		integer ntime	
		real rt, tt

	integer i
	real p
c-----
c	search through the list
c-----
	do 1000 i=1,ntime-1
		if(dist.ge.rt(i) .and.dist.le.rt(i+1))then
			p = (rt(i+1) - dist)/(rt(i+1) - rt(i))
			dur = p*tt(i) + (1.0-p)*tt(i+1)
		endif
 1000	continue
	return
	end

	subroutine getsit(v,freq)
c-----
c	Obtain the site amplification as a function of frequency
c
c	v	Real	- Site amplification
c	freq	Real	- frequency in Hz
c-----
	real v, freq
	parameter(NDMAX=100)
	common/props/nsite, fs(NDMAX),  ss(NDMAX), kappa
		integer nsite
		real fs, ss, kappa

	integer i
	real p
c-----
c	search through the list
c-----
	if(freq .lt. fs(1))then
		v = ss(1)
	else if(freq.gt.fs(nsite))then
		v = ss(nsite)
	else
		do 1000 i=1,nsite-1
			if(freq.ge.fs(i) .and.freq.le.fs(i+1))then
				p = (fs(i+1) - freq)/(fs(i+1) - fs(i))
				v = p*ss(i) + (1.0-p)*ss(i+1)
			endif
1000		continue
	endif
	return
	end

	subroutine getgeo(gr, dist)
c-----
c	obtain the G(r) geometrical spreading function
c----
	real gr, dist
	parameter(NDMAX=100)
	common/propg/ndist, rr(NDMAX),rpow(NDMAX)
		integer ndist
		real rr, rpow
	integer i
c-----
c	this is little complicated since we must carry through everything
c	we assume gr = 1.0 at 1.0 km
c-----
	r = dist
	if(r.eq.0.0)r = 1
	gr = 1.0
	if(ndist.eq.1)then
		gr = gr * (r/1.0)**rpow(1)
	else
		do 1000 i=1,ndist-1
			if(r.ge.rr(i))then
				if( r.ge.rr(i+1))then
					gr = gr * (rr(i+1)/rr(i))**rpow(i)
				else if(r.lt.rr(i+1))then
					gr = gr * (r/rr(i))**rpow(i)
				endif
			endif
 1000		continue
		if(r.gt.rr(ndist))then
			gr = gr * (r/rr(ndist))**rpow(ndist)
		endif
	endif
	return
	end

	subroutine initp()
c-----
c	initialize everything for the case of no control file
c	or a limited input
c-----
c-----
c	common blocks
c-----
	parameter(NDMAX=100)
	common/propq/q0,eta,vel
		real q0,eta,vel
	common/propg/ndist, rr(NDMAX),rpow(NDMAX)
		integer ndist
		real rr, rpow
	common/propd/ntime, rt(NDMAX),  tt(NDMAX)
		integer ntime	
		real rt, tt
	common/props/nsite, fs(NDMAX),  ss(NDMAX), kappa
		integer nsite
		real fs, ss, kappa
	common/srcfm/fmax, jsrc, sigma, betas, dens
		real fmax, sigma, betas, dens
		integer jsrc
	common/proc/comment
		character comment*80

	q0 = 10000.0
	eta = 0.0
	vel = 3.5
	ndist = 1
	rr(1) = 1.0
	rpow(1) = -1.0
	ntime = 2
	rt(1) = 0.0
	rt(2) = 10000.0
	tt(1) = 0.0
	tt(2) = 0.0
	nsite = 2
	fs(1) = 0.0
	fs(2) = 10000.0
	ss(1) = 1.0
	ss(2) = 1.0
	kappa = 0.0
	fmax = 100.0
	jsrc = 1
	sigma = 1.0
	betas = 3.5
	rhos = 2.7
	comment ='Default parameters'
	return
	end

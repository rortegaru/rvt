c------
c	Cartwright and Longuet-Higgens for bounds
c-----

	function qinv(prob,epsilon,xn)
c-----
c	obtain inverse of qetamx(epsilon,eta,xn) = prob
c	using interval halving
c-----
	eta0 = -20.0
	f0  = qetamx(epsilon,eta0,xn)-prob
	eta1 = 40.0
	f1  = qetamx(epsilon,eta1,xn)-prob
	niter = 0
	itmx = 25
 1000	continue
		eta2 = 0.5*(eta0 + eta1)
		f2  = qetamx(epsilon,eta2,xn)-prob
		if(sign(1.0,f2) .eq. sign(1.0,f0))then
			f0 = f2
			eta0 = eta2
		else
			f1 = f2
			eta1 = eta2
		endif
		niter = niter + 1
		if(niter.gt.itmx)go to 9999
	go to 1000
 9999	continue
	eta = 0.5*(eta0+eta1)
	qinv = eta
	return
	end

	function qetamx(epsilon,eta,xn)
c----
c	integrate formula 6.2
c-----
	if(xn.le.1.0)then
		xnn = 1.002
	else
		xnn = xn
	endif
C		y = xnn*alog((1.0 - qclh(epsilon,eta)))
C		qetamx = exp(y)
	qetamx = ( 1.0 - qclh(epsilon,eta))**xnn
	return
	end

	function qclh(epsilon,eta)
c-----
c	Cartwright and Longuet-Higgins 5.2
c-----
	if(epsilon.eq.1.0)then
		eps = 0.999
	else
		eps = epsilon
	endif
	eroot = sqrt(1.0 - eps**2)
	efac = 0.5*eta*eta
	fac = eroot * expchk(-efac)
	if(eta.gt.0.0)then
		if(epsilon.eq.0.0)then
			qclh = fac
		else
			qclh = qq(eta/eps) + fac*pp(eta*eroot/eps)
		endif
	else
		if(eps.eq.0.0)then
			qclh = 1.0
		else
			qclh = qq(eta/eps) + fac*pp(eta*eroot/eps)
		endif
	endif
	return
	end

	function expchk(x)
c-----
c	for negative exponents, verify to avoid underflow
c-----
	if(x.lt.-88.0)then
		expchk = 0.0
	else
		expchk = exp(x)
	endif
	return
	end

	function pp(xval)
c-----
c	26.2.19 Abromowitz and Stegun
c-----
		if(xval.lt.0)then
			x = - xval
		else
			x = xval
		endif
		pp = 1.0 -0.5*(
     1			1.0 +x*(0.0498673470 
     2			+x*(0.0211410061
     3			+x*(0.0032776263
     4			+x*(0.0000380036
     5			+x*(0.0000488906
     6			+x*(0.0000053830)))))))**(-16)
		if(xval.lt.0)then
			pp = 1.0 - pp
		endif
	return
	end

	function qq(x)
c-----
c
c-----
		qq = 1 - pp(x)
	return
	end


f1 = open('data11.dat','w')
f1.write('#Name\tTeff\tlog(g)\t[M/H]\tVturb\tPeriod\tP+err\tP-err\tTc\ta/R*\ta+err\ta-err\tEccentricity\tOmega\tRp/R*\tr+err\tr-err\tTc-err\tRA\tDEC\n')

a = True

while a:
	name = input('Enter the name of object: ')
	ra = input('Enter RA: ')
	dec = input('Enter DEC: ')
	teff = input('Enter the effective temperature of the star: ')
	lg = input('Enter the log(g) of the star: ')
	mh = input('Enter the metallicity: ')
	vturb = input('Enter the micro-turbulant velocity: ')
	p = input('Enter the period: ')
	pperr = input('Enter the +err in period: ')
	pnerr = input('Enter the -err in period: ')
	tc = input('Enter the Tc: ')
	tcerr = input('Enter the error in Tc: ')
	a = input('Enter a/R*: ')
	aperr = input('Enter the +err in a/R*: ')
	anerr = input('Enter the -err in a/R*: ')
	ecc = input('Enter the eccentricity of the planet: ')
	omega = input('Enter omega: ')
	r = input('Enter Rp/R*: ')
	rperr = input('Enter the +err in Rp/R*: ')
	rnerr = input('Enter the -err in Rp/R*: ')
	f1.write(name + '\t' + teff + '\t' + lg + '\t' + mh + '\t' + vturb + '\t' + p + '\t' + pperr + '\t' + pnerr + '\t' + tc + '\t' + a + '\t' + aperr + '\t' + anerr + '\t' + ecc + '\t' + omega + '\t' + r + '\t' + rperr + '\t' + rnerr + '\t' + tcerr + '\t' + ra + '\t' + dec + '\n')
	b = input('If you still wish to continue then press any key, otherwise type NO')
	if b == 'NO':
		a = False
	else:
		continue

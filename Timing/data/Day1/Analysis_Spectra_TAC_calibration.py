# Tac calibration
# File1: Spectra_TAC_calibration.root
# File2: calibr_TAC_2.root

#Peaks: 
#13ns
#(2235 +- 20) ch

#17ns
#(4950 +- 30) ch

#21ns
#(7555 +- 30) ch 

#25ns
#(10390 +- 30) ch

#29ns
#(13080 +- 40) ch

import numpy
import scipy.optimize as op
t = numpy.array([13, 17, 21, 25, 29])
ch = numpy.array([2235, 4950, 7555, 10390, 13080])
s = numpy.array([20, 30, 30, 30, 40])

def f(x,a,b):
	return a*x+b

popt, pcov = op.curve_fit(f, t, ch, sigma = s)

print("channel = a * t + b")
print("[a,b] = ", popt)
print("cov: ", pcov)

# channel = a * t + b
# ('[a,b] = ', array([  677.20357145, -6578.52261946]))
# ('cov: ', array([[   18.21537828,  -340.21625909],
#        [ -340.21625909,  6920.47269011]]))

# SPECTRA ANALYSIS
# Fit of Compton Edges
#Fit function: gauss + costant

#Peak 511
#chi/cdf = 116.6/131
 # FCN=116.586 FROM MIGRAD    STATUS=CONVERGED     162 CALLS         163 TOTAL
 #                     EDM=2.68033e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
 #  EXT PARAMETER                                   STEP         FIRST   
 #  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
 #   1  p0           5.88485e+03   7.63332e+01   5.14765e-02  -1.70591e-05
 #   2  p1           3.49333e+03   5.59523e+00   6.30375e-03   8.70258e-05
 #   3  p2           5.09088e+02   9.20098e+00   5.19157e-03   3.30287e-04
 #   4  p3           1.92130e+03   7.51779e+01   3.32337e-02   2.02921e-05

#Peak 1275
#chi/ndf = 378.3/282
 # FCN=378.279 FROM MIGRAD    STATUS=CONVERGED      57 CALLS          58 TOTAL
 #                     EDM=7.07964e-14    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   0.9 per cent
 #  EXT PARAMETER                                   STEP         FIRST   
 #  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
 #   1  p0           1.57780e+03   6.35556e+00  -1.64792e-05  -3.51732e-09
 #   2  p1           1.13404e+04   8.01616e+00   4.55914e-05   7.72766e-08
 #   3  p2           7.60665e+02   6.51891e+00  -3.76095e-05   1.90247e-08
 #   4  p3           6.36067e+01   3.36907e+00   1.40714e-05  -4.86313e-09













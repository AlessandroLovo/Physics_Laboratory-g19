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

# SPECTRA ANALYSYS, without background (work well)
# DETECTOR 1
# peak 1
#   1  Constant     6.35370e+03   8.87729e+00   6.59111e-02   3.45940e-05
#   2  Mean         3.43230e+03   3.19438e+00   1.16178e-02   5.97373e-05
#   3  Sigma        7.04745e+02   3.15905e+00   8.01438e-06   1.99382e-01
# peak 2
#   1  Constant     1.45289e+03   4.19882e+00   2.31654e-02  -2.34716e-04
#   2  Mean         1.13456e+04   7.88697e+00   2.08253e-02  -1.32533e-04
#   3  Sigma        7.71133e+02   8.27999e+00   1.47403e-05   6.60048e-02

# DETECTOR 2
# peak 1
#    1  Constant     4.20344e+03   6.33246e+00   5.58663e-02   1.13715e-05
#   2  Mean         4.28586e+03   3.97333e+00   1.77435e-02   1.68264e-05
#   3  Sigma        8.68150e+02   3.81496e+00   9.25942e-06   6.38889e-02
# peak 2
#   1  Constant     1.18921e+03   3.33149e+00   1.89897e-02  -5.45447e-06
#   2  Mean         1.33434e+04   6.14427e+00   2.13684e-02   3.26011e-06
#   3  Sigma        8.44520e+02   6.43728e+00   1.15997e-05   4.05863e-03








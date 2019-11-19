# TR1
   1  Constant     4.45658e+03   1.19642e+01   4.85013e-02  -1.77470e-06
   2  Mean         2.11258e+03   6.65762e-01   2.86280e-03  -3.18462e-04
   3  Sigma        2.02150e+02   6.85377e-01   3.36796e-06   3.77568e-01

# TR2
# chi/ndf = 166/61
#  FCN=166.795 FROM MIGRAD    STATUS=CONVERGED      44 CALLS          45 TOTAL
#                      EDM=1.83921e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
#   EXT PARAMETER                                   STEP         FIRST   
#   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
#    1  Constant     4.74726e+03   1.28458e+01   6.51901e-02  -1.25900e-05
#    2  Mean         1.34000e+03     fixed    
#    3  Sigma        1.94361e+02   6.01840e-01   7.04270e-06  -1.71325e-01

# TR3
   1  Constant     4.39535e+03   1.14710e+01   5.65922e-02  -2.50118e-08
   2  Mean         2.04678e+03   4.59919e-01   3.00540e-03  -1.55691e-06
   3  Sigma        2.00755e+02   4.69555e-01   4.16688e-06  -1.12895e-03


# TR4
   1  Constant     4.87560e+03   1.27054e+01   6.13723e-02  -3.40513e-06
   2  Mean         1.68339e+03   5.05689e-01   3.03713e-03   3.31880e-05
   3  Sigma        1.93573e+02   5.45573e-01   3.83491e-06  -1.63798e-01


# TR5
   1  Constant     4.96034e+03   1.19373e+01   2.51369e-01   1.38413e-07
   2  Mean         2.16586e+03   4.15646e-01   1.12617e-02  -7.67345e-06
   3  Sigma        1.94198e+02   3.64463e-01   6.98733e-06  -4.92415e-01

# TR6
   1  Constant     4.83558e+03   1.38085e+01   4.46037e-02   1.71003e-06
   2  Mean         1.66075e+03   9.24678e-01   2.64873e-03  -4.09938e-05
   3  Sigma        1.90012e+02   1.06343e+00   3.72324e-06   5.25918e-02

# TR7
   1  Constant     4.80276e+03   1.25077e+01   5.47874e-02  -4.58946e-05
   2  Mean         1.66524e+03   5.41045e-01   4.99571e-03   8.48929e-04
   3  Sigma        1.93677e+02   5.56735e-01   4.49775e-06  -1.03383e-01

#########################
import numpy

# channel = a * t + b
a = 677.2140816
var_a = 18.21537828

def ChannelToDealy( ch,  std):
	var = std**2
	sigma = ch/a
	std_sigma = sigma**2*((std/ch)**2+var_a/(a**2))
	print('[',sigma," +- ", numpy.sqrt(std_sigma), ']')


print("Sigma values of the five gaussians:")
ChannelToDealy(2.02150e+02,6.85377e-01)
ChannelToDealy(2.00755e+02,4.69555e-01) 
ChannelToDealy(1.93573e+02,5.45573e-01)
ChannelToDealy(1.94198e+02,3.64463e-01)
ChannelToDealy(1.90012e+02,1.06343e+00)
ChannelToDealy(1.93677e+02,5.56735e-01)


# Sigma values of the five gaussians:
# ('[', 0.2663648097420188, ' +- ', 0.0019405400413421009, ']')
# ('[', 0.28700082482159656, ' +- ', 0.0020152740429078319, ']')
# ('[', 0.29735943990447583, ' +- ', 0.0021707552394992506, ']')
# ('[', 0.2673748300865219, ' +- ', 0.0019051650462883118, ']')
# ('[', 0.26237198077778423, ' +- ', 0.002054248769799013, ']')

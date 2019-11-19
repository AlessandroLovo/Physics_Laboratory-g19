# TR1
# chi/ndf = 108/68
#  FCN=108.547 FROM MIGRAD    STATUS=CONVERGED      85 CALLS          86 TOTAL
#                      EDM=2.86144e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
#   EXT PARAMETER                                   STEP         FIRST   
#   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
#    1  Constant     5.07164e+03   1.29676e+01   5.46060e-02   4.37922e-05
#    2  Mean         2.18132e+03   8.42589e-01   2.46254e-03  -6.44502e-04
#    3  Sigma        1.80386e+02   6.59268e-01   3.90490e-06  -2.70571e-01

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
# chi/ndf = 201.3/78
#  FCN=201.296 FROM MIGRAD    STATUS=CONVERGED      86 CALLS          87 TOTAL
#                      EDM=2.53709e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
#   EXT PARAMETER                                   STEP         FIRST   
#   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
#    1  Constant     4.38845e+03   1.14764e+01   6.51890e-02   3.97369e-05
#    2  Mean         2.04805e+03   9.35902e-01   3.76264e-03  -6.59948e-04
#    3  Sigma        2.01376e+02   7.41918e-01   5.34731e-06  -3.44486e-01

# TR4
# chi/ndf = 115/71
#  FCN=115.305 FROM MIGRAD    STATUS=CONVERGED      83 CALLS          84 TOTAL
#                      EDM=6.64222e-08    STRATEGY= 1      ERROR MATRIX ACCURATE 
#   EXT PARAMETER                                   STEP         FIRST   
#   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
#    1  Constant     5.07397e+03   1.29030e+01   5.50244e-02  -2.95302e-05
#    2  Mean         2.18026e+03   7.29918e-01   2.47708e-03   1.56715e-04
#    3  Sigma        1.81070e+02   6.02016e-01   3.89362e-06  -9.24883e-02

# TR5
# chi/ndf = 67/58
#  FCN=67.0633 FROM MIGRAD    STATUS=CONVERGED      82 CALLS          83 TOTAL
#                      EDM=4.44423e-11    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.3 per cent
#   EXT PARAMETER                                   STEP         FIRST   
#   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
#    1  Constant     5.07405e+03   1.35467e+01  -2.58996e-02  -9.00771e-07
#    2  Mean         2.18424e+03   1.04278e+00   4.69843e-04   2.88055e-06
#    3  Sigma        1.77682e+02   8.25478e-01   1.86566e-06  -2.78659e-04





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
ChannelToDealy(1.80386e+02, 6.59268e-01)
ChannelToDealy(1.94361e+02, 6.01840e-01) 
ChannelToDealy(2.01376e+02, 7.41918e-01)
ChannelToDealy(1.81070e+02, 6.02016e-01)
ChannelToDealy(1.77682e+02, 8.25478e-01)


# Sigma values of the five gaussians:
# ('[', 0.2663648097420188, ' +- ', 0.0019405400413421009, ']')
# ('[', 0.28700082482159656, ' +- ', 0.0020152740429078319, ']')
# ('[', 0.29735943990447583, ' +- ', 0.0021707552394992506, ']')
# ('[', 0.2673748300865219, ' +- ', 0.0019051650462883118, ']')
# ('[', 0.26237198077778423, ' +- ', 0.002054248769799013, ']')
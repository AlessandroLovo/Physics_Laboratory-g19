# Analysis of LEMOs done measuring the delay of each LEMO plus 13ns (-> Base channel 13ns = (2230 +- 20) ch)
# Calibration of TAC: 
# channel = a * t + b
# ('[a,b] = ', array([  677.2140816 , -6576.89957312]))
# ('cov: ', array([[   18.21537828,  -340.21625909],
#        [ -340.21625909,  6920.47269011]]))

# LEMO_4
# (5700 +- 10)ch


# LEMO_5

# (5680 +- 40)ch

# LEMO_6

# (3025 +- 20)ch

# LEMO_7

# (3020 +- 20)ch

# LEMO_A13

# (4070 +- 30)ch

import numpy

a = 677.2140816
b = -6576.89957312
var_a = 18.21537828
var_b = 6920.47269011

def ChannelToDealy( ch,  std):
	var = std**2
	diff = ch - b
	delay = diff/a
	var_d = delay**2*((var+var_b)/(diff**2)+var_a/(a**2))
	print('[',delay-13," +- ", numpy.sqrt(var_d), ']')

ChannelToDealy(5700, 10)
ChannelToDealy(5680, 40)
ChannelToDealy(3025, 20)
ChannelToDealy(3020, 20)
ChannelToDealy(4070, 30)

# Results:
# ('[', 5.1285355793464085, ' +- ', 0.16840693955427827, ']')
# ('[', 5.099002820735205, ' +- ', 0.17773320599683978, ']')
# ('[', 1.1785291150980726, ' +- ', 0.15474660175065477, ']')
# ('[', 1.1711459254452716, ' +- ', 0.15471973810651732, ']')
# ('[', 2.7216157525334026, ' +- ', 0.16391847133503851, ']')



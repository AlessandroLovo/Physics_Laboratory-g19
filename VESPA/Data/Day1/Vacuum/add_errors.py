import sys

filename_in = sys.argv[1]
filename_out = sys.argv[2]

ifile = open(filename_in,'r')
ofile = open(filename_out,'w')

error = 1. # last digit error

while True:
    s = ifile.readline().rstrip('\n').rstrip(' ')
    if len(s) == 0:
        break
    sx,space,sy = s.partition(' ')
    sx = sx..strip(' ')
    sy = sy.strip(' ')
    x = float(sx)
    y = float(sy)
    
    i,dot,dec = sy.partition('.')
    err = error*10**(-len(dec))
    
    ofile.write(str(x)+'\t'+str(y)+'\t0\t'+str(err)+'\n')
    
ifile.close()
ofile.close()
    

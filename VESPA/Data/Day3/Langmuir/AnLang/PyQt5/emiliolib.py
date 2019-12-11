# -*- coding: utf-8 -*-
"""
Created on Sat May 09 22:15:41 2015

@author: Emilio Martines
"""

def readYoko(fileName):
    import struct
    import numpy as np
    
    try:
        f = open(fileName+'.HDR', 'r')
    except:
        try:
            f = open(fileName+'.hdr', 'r')
        except:
            print("Could not open Yokogawa header file.")
            return False
    lines = f.readlines()
    f.close()
    
    print("Reading header file")    
    if lines[0] != '//YOKOGAWA ASCII FILE FORMAT\n':
        print("Error in Yokogawa header file.")
        return False

    lines = lines[3:]
    map(str.strip, lines)
#    [line.strip() for line in lines]
    for line in lines:
        if line.strip() == "": lines.remove(line)

    properties = {}
    for line in lines[0:6]:
        col = line.split()
        properties[col[0]] = col[1]
    lines = lines[7:]
#    print properties

    gproperties = {}
    for igroup in range(int(properties["GroupNumber"])):
        for line in lines[1:3]:
            col = line.split()
            gproperties[col[0]] = col[1]
#        print gproperties
        for line in lines[3:17+int(gproperties["BlockNumber"])*2]:
#            print line
            col=line.split()
            if col[0] == "HResolution":
                hresolution = float(col[1])
            elif col[0] == "HOffset":
                hoffset = float(col[1])
            elif col[0] == "VResolution":
#                print col[1:]
                vresolution = np.array(col[1:]).astype(float)
            elif col[0] == "VOffset":
                voffset = np.array(col[1:]).astype(float)
            elif col[0] == "BlockSize":
                blocksize = int(col[1])

    ntracestot = int(properties['TraceTotalNumber'])
#    blocknumber = int(gproperties['BlockNumber'])
        
    try:
        f = open(fileName+'.WVF', 'rb')
    except:
        try:
            f = open(fileName+'.wvf', 'rb')
        except:
            print("Could not open Yokogawa waveform file.")
            return False
    
    print("Reading data from model", properties["Model"])
    t = np.linspace(hoffset, hoffset+hresolution*blocksize, num=blocksize)    
#    print(blocksize,ntracestot)
    signal = np.fromfile(f, np.int16, count =
        blocksize*ntracestot).reshape((blocksize,ntracestot), order='F')
    if struct.pack('=f', 2.3) == struct.pack('<f', 2.3):
        signal.byteswap(True)
    signal = signal.astype(float) * vresolution + voffset
   
    return signal, t


if __name__ == "__main__":
    print("sig, t = readYoko(fileName)")

    import matplotlib.pyplot as plt
    data1,t1 = readYoko('LANGM092')
    data2,t2 = readYoko('LANGM097')
    print(data1.shape)
    #plt.plot(t,data[:,0])
    plt.plot(t1,data1[:,1])
    plt.plot(t2,data2[:,1])
    plt.show()

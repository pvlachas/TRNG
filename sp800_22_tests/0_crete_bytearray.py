#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""Created by: Vlachas Pantelis, CSE-lab, ETH Zurich
"""
#!/usr/bin/env python
import numpy as np

# data = np.random.rand(2000000)
data = np.random.rand(2000000)

data[data<0.5]=0
data[data>=0.5]=1
data = [int(d) for d in data]
data = np.array(data)


print(np.sum(data==0))
print(np.sum(data==1))

s = "".join(["{:}".format(d) for d in data])

i = 0
buffer = bytearray()
while i < len(s):
    buffer.append( int(s[i:i+8], 2) )
    i += 8

file = "temp.bin"
# now write your buffer to a file
with open(file, 'bw') as f:
    f.write(buffer)


# from . import sp800_22

import sp800_22


result = sp800_22.runTestRoutines(file)

print(result)





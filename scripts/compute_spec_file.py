import math
import numpy as np
import sys

dat_name = sys.argv[1]

f_in = open(dat_name+".dat", 'r')

ok = 1
minitem=99
maxitem=0
maxtranslength=0
numtrans=0
total = 0
while (ok==1) :
	line=f_in.readline()
	items=line.split( )

	for i in items:
		if int(i)<minitem:
			minitem=int(i)
		else:
			if int(i)>maxitem:
				maxitem=int(i)
	total = total + len(items)
	if len(items)>0:

		if len(items)>maxtranslength:
			maxtranslength=len(items)

		numtrans=numtrans+1

	else:
		ok = 0


f_in.close()

print 'number of transactions: '+str(numtrans)
print 'minitem: '+str(minitem)
print 'maxitem: '+str(maxitem)
print 'maxtranslength: '+str(maxtranslength)
print 'total = '+str(total)
print 'avg '+str(float(total) / float(numtrans))

#write the spec file
alpha = 0.05
jp = 10000
newspecfilename = dat_name + "_new.spec"
f_out_spec = open(newspecfilename,'w')
f_out_spec.write(dat_name + ".dat\n")
f_out_spec.write(str(maxitem+1) + "\n")
f_out_spec.write(str(maxtranslength+1) + "\n")
f_out_spec.write(str(numtrans+1) + "\n")
f_out_spec.write(dat_name + ".labels" + "\n")
f_out_spec.write(str(alpha) + "\n")
f_out_spec.write(str(jp) + "\n")
f_out_spec.close()
print "spec file done for "+newspecfilename

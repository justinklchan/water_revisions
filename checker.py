dd='9600_test/'
for i in range(5):
	elts = open (dd+'alice/Alice-ValidBins-'+str(i)+'.txt').read().split(',')
	elts2 = open (dd+'alice/Alice-BitsAdapt_Padding-'+str(i)+'.txt').read().split(',')
	print (len(elts),len(elts2))
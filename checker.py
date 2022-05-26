import os
import sys
dd='Air_data4/'+sys.argv[1]+'/'
for i in range(5):
	elts = open (dd+'alice/Alice-ValidBins-'+str(i)+'.txt').read().split(',')
	elts2 = open (dd+'alice/Alice-BitsAdapt_Padding-'+str(i)+'.txt').read().split(',')
	elts3 = open (dd+'alice/Alice-BitsAdapt-'+str(i)+'.txt').read().split(',')
	elts4 = open (dd+'alice/Alice-FeedbackFreqs-'+str(i)+'.txt').read().split('\n')
	elts5 = open (dd+'bob/Bob-FreqEsts-'+str(i)+'.txt').read().split('\n')
	e = os.path.isfile (dd+'alice/Alice-ExactFeedbackFreqs-'+str(i)+'.txt')
	# print (elts4)
	eq1=elts4[0]==elts5[0]
	eq2=elts4[1]==elts5[1]
	print (len(elts),len(elts2),len(elts3),eq1,eq2,e)
	print (elts4)
	print (elts5)
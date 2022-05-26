import os
import sys
dd='Air_data5/'+sys.argv[1]+'/'
if sys.argv[1]=='960':
	val=50
elif sys.argv[1]=='1920':
	val=25
elif sys.argv[1]=='4800':
	val=10
elif sys.argv[1]=='9600':
	val=5
for i in range(5):
	elts = open (dd+'alice/Alice-ValidBins-'+str(i)+'.txt').read().split(',')
	elts2 = open (dd+'alice/Alice-BitsAdapt_Padding-'+str(i)+'.txt').read().split(',')
	elts3 = open (dd+'alice/Alice-BitsAdapt-'+str(i)+'.txt').read().split(',')
	elts4 = open (dd+'alice/Alice-FeedbackFreqs-'+str(i)+'.txt').read().split('\n')
	elts5 = open (dd+'bob/Bob-FreqEsts-'+str(i)+'.txt').read().split('\n')
	e = os.path.isfile (dd+'alice/Alice-ExactFeedbackFreqs-'+str(i)+'.txt')
	elts6 = open (dd+'alice/Alice-ExactFeedbackFreqs-'+str(i)+'.txt').read().split('\n')
	elts6=elts6[:len(elts6)-1]
	modvals=[1 if int(i)%val==0 else 0 for i in elts6]
	modresult=sum(modvals)==len(modvals)

	# print (elts6)
	eq1=elts4[0]==elts5[0]
	eq2=elts4[1]==elts5[1]
	print (len(elts),len(elts2),len(elts3),eq1,eq2,e,modresult)
	print (elts4)
	print (elts5)
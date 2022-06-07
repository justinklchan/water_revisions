import os
import matplotlib.pyplot as plt
import numpy as np
import shutil


def read_timstamp(name):
    text_file = open(name, "r")
    lines = text_file.readlines()
    time_list = []
    for t in lines:
        time_list.append(int(t))
    
    new_time_list = []
    for i in range(1, len(time_list)):
        new_time_list.append((time_list[i] - time_list[0])/1000)
    return new_time_list


def move_files(folder, pair_index):
    new_folder = folder + '/sync_file'
    if not os.path.exists(new_folder):                   
        os.makedirs(new_folder) 

    for i in range(0, len(pair_index)):
        Alice_index = pair_index[i][0]
        Bob_index = pair_index[i][1]
        # Alice data
        folder1 = folder + "/alice"
        shutil.copy(folder1+ '/Alice-BitsAdapt-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsAdapt-' + str(i)+'.txt')
        shutil.copy(folder1+ '/Alice-Bit_Fill_Adapt-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-Bit_Fill_Adapt-' + str(i)+'.txt')
        shutil.copy(folder1+ '/Alice-BitsAdapt_Padding-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsAdapt_Padding-' + str(i)+'.txt')
        
        
        shutil.copy(folder1+'/Alice-BitsFull_1000_4000-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsFull_1000_4000-' + str(i)+'.txt')
        shutil.copy(folder1+ '/Alice-Bit_Fill_1000_4000-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-Bit_Fill_1000_4000-' + str(i)+'.txt')        
        shutil.copy(folder1+'/Alice-BitsFull_1000_4000_Padding-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsFull_1000_4000_Padding-' + str(i)+'.txt')
        

        shutil.copy(folder1+'/Alice-BitsFull_1000_2500-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsFull_1000_2500-' + str(i)+'.txt')
        shutil.copy(folder1+ '/Alice-Bit_Fill_1000_2500-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-Bit_Fill_1000_2500-' + str(i)+'.txt')
        shutil.copy(folder1+'/Alice-BitsFull_1000_2500_Padding-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsFull_1000_2500_Padding-' + str(i)+'.txt')

        shutil.copy(folder1+'/Alice-BitsFull_1000_1500-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsFull_1000_1500-' + str(i)+'.txt')
        shutil.copy(folder1+ '/Alice-Bit_Fill_1000_1500-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-Bit_Fill_1000_1500-' + str(i)+'.txt')        
        shutil.copy(folder1+'/Alice-BitsFull_1000_1500_Padding-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-BitsFull_1000_1500_Padding-' + str(i)+'.txt')

        shutil.copy(folder1+'/Alice-DataAdapt-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-DataAdapt-' + str(i)+'.txt')
        shutil.copy(folder1+'/Alice-DataFull_1000_4000-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-DataFull_1000_4000-' + str(i)+'.txt')
        shutil.copy(folder1+'/Alice-DataFull_1000_2500-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-DataFull_1000_2500-' + str(i)+'.txt')
        shutil.copy(folder1+'/Alice-DataFull_1000_1500-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-DataFull_1000_1500-' + str(i)+'.txt')
        shutil.copy(folder1+'/Alice-FeedbackFreqs-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-FeedbackFreqs-' + str(i)+'.txt')
        shutil.copy(folder1+'/Alice-ExactFeedbackFreqs-' + str(Alice_index)+'.txt', folder+'/sync_file/Alice-ExactFeedbackFreqs-' + str(i)+'.txt')

        #Bob data
        folder2 = folder + "/bob"
        shutil.copy(folder2+'/Bob-DataRx-' + str(Bob_index)+'-0-bottom.txt', folder+'/sync_file/Bob-DataRx-' + str(i)+'-bottom.txt')
        shutil.copy(folder2+'/Bob-DataRx-' + str(Bob_index)+'-0-top.txt',folder+ '/sync_file/Bob-DataRx-' + str(i)+'-top.txt')
        shutil.copy(folder2+'/Bob-FreqEsts-' + str(Bob_index)+'.txt', folder+'/sync_file/Bob-FreqEsts-' + str(i)+'.txt')
        shutil.copy(folder2+'/Bob-SNRs-' + str(Bob_index)+'.txt',folder+ '/sync_file/Bob-SNRs-' + str(i)+'.txt')
        shutil.copy(folder2+'/Bob-Sounding-' + str(Bob_index)+'-0-top.txt', folder+'/sync_file/Bob-Sounding-' + str(i)+'-top.txt')
        shutil.copy(folder2+'/Bob-Sounding-' + str(Bob_index)+'-0-bottom.txt', folder+'/sync_file/Bob-Sounding-' + str(i)+'-bottom.txt')


folder = 'Air_data20/4800'
sounding_time = read_timstamp(folder + '/alice/Alice-Sounding-log.txt')
data_time = read_timstamp(folder + '/alice/Alice-Data-log.txt')
feedback_time = read_timstamp(folder + '/bob/Bob-Feedback-log.txt')

y_feedback = np.ones(len(feedback_time))
y_sounding = np.ones(len(sounding_time))
y_data = np.ones(len(data_time))

min_difference = 1000
min_pair = []
min_offset = 0

plt.figure(1)
plt.subplot(211)
plt.stem(sounding_time, y_sounding,  'r-x', label='sounding')
plt.stem(np.array(data_time), y_data, 'g--^', label='data')
plt.xlim([0, 130])
plt.subplot(212)
plt.stem(feedback_time, y_feedback, 'b-o',label='feedback')
plt.xlim([0, 130])
plt.legend()
interval = 5

print(feedback_time)
print(data_time)

for i in range(0, len(feedback_time)):
    if((feedback_time[i] - data_time[0]) > interval or (feedback_time[i] - data_time[0]) < -interval): continue
    print("-------------")
    data_time_new = np.array(data_time) + feedback_time[i] - data_time[0] 
    pair_index = []

    difference = 0
    for j in range(0, len(data_time)):
        delta_list = np.absolute(np.array(feedback_time) - data_time_new[j])
        match_idx = np.argmin(delta_list)     
        difference += np.min(delta_list)
        pair_index.append([j, match_idx])
    
    print(i, feedback_time[i] - data_time[0], difference)
    if difference < min_difference:
        min_difference = difference
        min_pair = pair_index
        min_offset = feedback_time[i] - data_time[0]


min_pair = []
for i in range(0, 10):
    min_pair.append([i, i])


print(min_pair)
plt.figure(2)
#plt.stem(sounding_time, y_sounding,  'r--x', label='sounding')
plt.stem(feedback_time, y_feedback, 'b-o',label='feedback')
plt.stem(np.array(data_time) + min_offset, y_data, 'g--^', label='data')
plt.legend()
plt.show()

move_files(folder, min_pair)


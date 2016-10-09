import os
import matplotlib.pyplot as plt
from matplotlib import animation 

data = "Output/"

def ReadInData(file):
	content = []
	file = open(file,'r')
	file = file.read().splitlines()
	for line in file:
		line = line.split(' ')
		row  = [float(item) for item in line if item != '']
		content.append(row)
	return content 

file_name_e	= [ data+fn for fn in os.listdir(os.getcwd()+"/"+data) if fn.endswith('.txt')]

a       = ReadInData(file_name_e[0])    # Reads a.txt
a_true  = ReadInData(file_name_e[1])    # Reads a_true.txt
r       = ReadInData(file_name_e[2])[0] # Reads r.txt



fig   = plt.figure(1)	
ax1   = fig.add_subplot(1,1,1)
scat1 = ax1.scatter(r,a[0],facecolors='none',edgecolors='k')
plt1, = ax1.plot(r,a_true[0])

def update(i,fig,ax1):
    scat1.set_offsets([[r[j],a[i][j]] for j in range(len(a[i])-1)])
    plt1.set_data(r,a_true[i])
    plt.suptitle('%i' % i)
    return ax1
    
anim = animation.FuncAnimation(fig,update,fargs=(fig,ax1),frames=len(a),interval=100)
plt.show()

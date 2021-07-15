''' Plotting Intensities vs frequencies using Vlines '''
import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np



d1=pd.read_csv( "2_AA.data" ,sep='\s+', header=0) #os.path other than writing the full path
d2=pd.read_csv( "4_AA.data" ,sep='\s+', header=0) #os.path other than writing the full path
d3=pd.read_csv( "6_AA.data" ,sep='\s+', header=0) #os.path other than writing the full path
d4=pd.read_csv( "usual_raman" ,sep='\s+', header=0) #os.path other than writing the full path


# Raman along z

plt.vlines(d1['Freq'],ymin=0, ymax=d1['Raman_zz'], colors='green', linestyles='solid',label=' 2 AA from tip', linewidth=4)
plt.vlines(d2['Freq']+20,ymin=0, ymax=d2['Raman_zz'], colors='blue', linestyles='solid',label=' 4 AA from tip', linewidth=4)
plt.vlines(d3['Freq']-20,ymin=0, ymax=d3['Raman_zz'], colors='red', linestyles='solid',label=' 6 AA from tip', linewidth=4)
plt.vlines(d4['Freq']+30,ymin=0, ymax=d4['Raman_zz'], colors='black', linestyles='solid',label=' Usual raman', linewidth=4)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid()
plt.xlabel('Frequency cm-1', fontsize=11)
plt.ylabel('Raman_z Intensity [Ang^4/amu]', fontsize= 11)
plt.legend(fontsize=10)
plt.tight_layout() #better
plt.savefig('z.png', dpi=400)
plt.show()

# Averange Raman

plt.vlines(d1['Freq'],ymin=0, ymax=d1['Total_Raman'], colors='green', linestyles='solid',label=' 2 AA from tip', linewidth=4)
plt.vlines(d2['Freq']+20,ymin=0, ymax=d2['Total_Raman'], colors='blue', linestyles='solid',label=' 4 AA from tip', linewidth=4)
plt.vlines(d3['Freq']-20,ymin=0, ymax=d3['Total_Raman'], colors='red', linestyles='solid',label=' 6 AA from tip', linewidth=4)
plt.vlines(d4['Freq']+30,ymin=0, ymax=d4['Total_Raman'], colors='black', linestyles='solid',label=' Usual raman', linewidth=4)
plt.legend(fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.grid()
plt.xlabel('Frequency cm-1', fontsize=11)
plt.ylabel('Total Raman Intensity [Ang^4/amu]', fontsize= 11)
plt.tight_layout() #better
plt.savefig('Total.png', dpi=400)
plt.show()




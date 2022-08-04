import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt


#data2 = pd.read_csv("Default Dataset.csv", names=["mu", "q"], header=0)
data = pd.read_csv("real_absorb_tau=1.csv", names=["mu", "q", "u", "N","w"], header=0)
data2 = pd.read_csv("real_absorb_tau=5.csv", names=["mu", "q", "u", "N","w"], header=0)
data3 = pd.read_csv("real_absorb_tau=10.csv", names=["mu", "q", "u", "N","w"], header=0)
data4 = pd.read_csv("real_absorb_tau=15.csv", names=["mu", "q", "u", "N","w"], header=0)


fig, ax1= plt.subplots()
#ax2 = ax1.twinx()
#ax3 = plt.subplot(sharey=ax1)

N0 = np.sum(data["N"])
N1 = np.sum(data2["N"])
N2 = np.sum(data3["N"])
N3 = np.sum(data4["N"])

print(N0,N1,N2,N3)


#plt.scatter(data["mu"], data["q"], s=16,
#          edgecolors="black", c='white', alpha=0.8, label="sample")

#p = np.sqrt((data["q"]/data["sum"])**2+(data["u"]/data["sum"])**2)
y = [0.11713, 0.08979, 0.07448, 0.06311, 0.05410, 0.04667, 0.04041, 0.03502, 0.03033, 0.02619,
     0.02252, 0.01923, 0.01627, 0.01358, 0.011123, 0.00888, 0.006818, 0.004919, 0.003155, 0.001522, 0]
x = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]  
y_ave = [0.10346, 0.082135, 0.068795, 0.058605, 0.050385, 0.04354, 0.037715, 0.032675, 0.02826, 
     0.024355, 0.020875, 0.01775, 0.014925, 0.0123515, 0.0100015, 0.007849, 0.0058685, 0.004037, 0.0023385, 0.000761]    



#for j in range(20):
ax1.scatter(data["mu"], -data["q"], s=5, edgecolors="red",
           color='red', alpha=0.8, label=r"$\tau_s=1$")
ax1.scatter(data2["mu"], -data2["q"], s=5, edgecolors="blue",
            color='blue', alpha=0.8, label=r"$\tau_s=5$")

ax1.scatter(data3["mu"], -data3["q"], s=5, edgecolors="green",
           color='green', alpha=0.8, label=r"$\tau_s=10$")
ax1.scatter(data4["mu"], -data4["q"], s=5, edgecolors="yellow",
           color='yellow', alpha=0.8, label=r"$\tau_s=15$")

#ax1.scatter(data4["mu"], -data4["u"], s=16, edgecolors="yellow",
#            color='yellow', alpha=0.8, label=r"$\tau=12$")

#ax1.scatter(data2["mu"], data2["q"], s=16, edgecolors="blue",
#            color='blue', alpha=0.8, label=r"Agol")
#ax1.scatter(data2["mu"], -data2["q"], s=16, edgecolors="blue",
#            color='blue', alpha=0.8, label=r"$\alpha=0.25,1.0$")


#ax2.scatter(data["mu"], -data["u"], s=16, edgecolors="black",
#            color='black', alpha=0.8, label="u")

#ax1.scatter(data["mu"], data["q"]/data["sum"], s=16, edgecolors="blue",
#            color='black', alpha=0.8, label="relative")

#ax1.scatter(data["mu"], data["q"]/data["sum"], s=16, edgecolors="black",
#           color='white', alpha=0.8, label="$q$")
#ax1.scatter(data["mu"], data["u"]/data["sum"], s=16, edgecolors="black",
#            color='black', alpha=0.8, label="$u$")
#ax1.scatter(data2["mu"], data2["q"]/data2["sum"], s=16, edgecolors="green",
 #           color='green', alpha=0.8, label="old")

#ax1.plot(x, y, label="chandrasekhar")




#ax2.scatter(data["mu"], abs((data["q"]/data["sum"]-y_ave)*100/y_ave), s=16, edgecolors="red",
#            color='red', alpha=0.8, label="difference")


ax1.set_xlabel("$\mu=cos(i)$")
ax1.set_ylabel("q")
ax1.grid(True)
#ax2.set_ylabel("u")

plt.title('$T=10^9$')
plt.xlim(0, 1)

handler1, label1 = ax1.get_legend_handles_labels()
#handler2, label2 = ax2.get_legend_handles_labels()
#handler3, label3 = ax3.get_legend_handles_labels()

ax1.legend(handler1, label1,
          loc="upper right", borderaxespad=0.)
#ax1.legend(handler1 + handler2, label1 + label2,
#           loc="upper right", borderaxespad=0.)
#plt.legend(("plot","histgram"))

ax1.set_ylim([-0.1, 0.12])
#ax2.set_ylim([-0.03, 0.03])
#ax2.set_ylim([0, 70])

plt.tight_layout()
plt.show()

#fig.savefig("beta_x_mu-q.png")
#fig.savefig("mu-Qi.png")
#fig.savefig("scatter2(q=-1).png")
#fig.savefig("mu-p.png")

#print(df)
#print(sum)
#print(df/sum)
#print(df_u/df)
#print(df_q/df)
#print((data["q"]-y))


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set()

history = pd.read_csv('./results/roul_ox1_sim/history.csv',
                      header=None)
distribution = pd.read_csv('./results/roul_ox1_sim/distribution.csv',
                           header=None)
bestpath = pd.read_csv('results/roul_ox1_sim/bestpath.csv',
                       header=None)

history[0] = history[0] * 100
distribution[0] = distribution[0] * 100

fig = plt.figure(figsize=[12, 6])

ax1 = fig.add_subplot(121)
ax1.set_title('fitness history')
ax1.set_xlabel('generation')
ax1.set_ylabel('fitness')
ax1.plot(history[0], history[1])
ax1.plot(distribution[0], distribution[1], marker='.', alpha=0.1, linestyle='')

ax2 = fig.add_subplot(122)
ax2.set_title('Best path')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.plot(bestpath[0], bestpath[1], marker='.')

plt.show()

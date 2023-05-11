import numpy as np
import matplotlib.pyplot as plt

data_on_node = np.loadtxt('On-Node——ping_pong .txt')
data_out_of_node = np.loadtxt('Out-of-Node——ping_pong.txt')

p_on_node = np.polyfit(data_on_node[:,0], data_on_node[:,1], 1)
p_out_of_node = np.polyfit(data_out_of_node[:,0], data_out_of_node[:,1], 1)
print("###### On-Node ######")
print(f"bandwidth: {1/p_on_node[0]/10e9:.4f} GB/s")
print(f"latency: {p_on_node[1]*10e6:.4f} us")
print("###### Out-of-Node ######")
print(f"bandwidth: {1/p_out_of_node[0]/10e9:.4f} GB/s")
print(f"latency: {p_out_of_node[1]*10e6:.4f} us")


plt.plot(data_on_node[:,0], data_on_node[:,1], label='Intra-node')
plt.plot(data_out_of_node[:,0], data_out_of_node[:,1], label='Inter-node')
plt.xlabel("Message size (bytes)")
plt.ylabel('Time (s)')
# plt.ylim(0, 2* 10e-5)
plt.legend()
plt.show()



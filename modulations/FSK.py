from matplotlib import pyplot as plt
import sys

def cpu_fsk():
    import numpy as np

    N = 8
    bit_stream = np.array([0, 0, 1, 1, 0, 1, 1, 0])

    f1 = 5
    f2 = 10
    fs = 2000
    time_step = 1/fs
    print("time_step="+str(time_step))


    # build time array
    t = [0]
    last = 0
    for i in range(fs):
        t.append(last + time_step)
        last = t[-1]
    t = np.array(t)

    # Digital 0 and 1 (for time t)
    zeros = np.zeros(len(t))
    ones = np.ones(len(t))

    # Modulated 0 and 1 (for time t at frequencies f1 and f2)
    a = np.sin(np.multiply((2 * np.pi * f1),t))
    b = np.sin(np.multiply((2 * np.pi * f2),t))

    # build the overall signal
    time = np.array([])
    digital_signal = np.array([])
    fsk_signal = np.array([])
    for i in range(len(bit_stream)):

        # Digital signal
        if bit_stream[i] == 0:
            digital_signal = np.append(digital_signal,zeros)
        if bit_stream[i] == 1:
            digital_signal = np.append(digital_signal,ones)
        
        # Modulated signal
        if bit_stream[i] == 0:
            fsk_signal = np.append(fsk_signal,a)
        if bit_stream[i] == 1:
            fsk_signal = np.append(fsk_signal,b)
       
       
        time = np.append(time,t)
        t = t + 1     
        
        #print(fsk_signal)
        #print(time)
        #sys.exit(0)

    plt.plot(time, digital_signal)
    plt.plot(time, fsk_signal)
    plt.show()


def gpu_fsk():
    import cupy as cp
    
    N = 8
    bit_stream = cp.array([0, 0, 1, 1, 0, 1, 1, 0])

    f1 = 5
    f2 = 10
    fs = 2000
    time_step = 1/fs
    print("time_step="+str(time_step))


    # build time array
    t = [0]
    last = 0
    for i in range(fs):
        t.append(last + time_step)
        last = t[-1]
    t = cp.array(t)

    # Digital 0 and 1 (for time t)
    zeros = cp.zeros(len(t))
    ones = cp.ones(len(t))

    # Modulated 0 and 1 (for time t at frequencies f1 and f2)
    a = cp.sin(cp.multiply((2 * cp.pi * f1),t))
    b = cp.sin(cp.multiply((2 * cp.pi * f2),t))

    # build the overall signal
    time = cp.array([])
    digital_signal = cp.array([])
    fsk_signal = cp.array([])
    for i in range(len(bit_stream)):

        # Digital signal
        if bit_stream[i] == 0:
            digital_signal = cp.append(digital_signal,zeros)
        if bit_stream[i] == 1:
            digital_signal = cp.append(digital_signal,ones)
        
        # Modulated signal
        if bit_stream[i] == 0:
            fsk_signal = cp.append(fsk_signal,a)
        if bit_stream[i] == 1:
            fsk_signal = cp.append(fsk_signal,b)
       
       
        time = cp.append(time,t)
        t = t + 1     
        
        #print(fsk_signal)
        #print(time)
        #sys.exit(0)

    plt.plot(time, digital_signal)
    plt.plot(time, fsk_signal)
    plt.show()



if __name__=="__main__":
    #cpu_fsk()
    gpu_fsk()










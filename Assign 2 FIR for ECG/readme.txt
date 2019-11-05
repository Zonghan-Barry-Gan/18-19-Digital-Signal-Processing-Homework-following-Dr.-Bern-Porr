the ecg_fir.py provide a set of real time FIR filter that designed for ECG data, the best filter is (0,35]HZ yet this rid of the base frequency of heart beat and change the shape of ECG curve
the ecg_count_hbt provied a real time match filter to count heart beat and out put frequency
The whole set of code is designed to realize real time analysis of a ecg device we assembled using USB-DUX. A link to video showing how to use this device would be added
the I/O of the USB-DUX running on Linux platform is provided by Dr. Bern 
(berndporr/pyusbdux, GitHub. (2019). berndporr/pyusbdux. [online] Available at: https://github.com/berndporr/pyusbdux [Accessed 5 Nov. 2019].)

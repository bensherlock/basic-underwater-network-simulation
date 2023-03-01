# Basic Underwater Acoustic Network Simulation

A basic underwater acoustic network simulation framework with example scenario. Two sensor nodes periodically send data to a gateway node, the gateway node replys with an acknowledgement, and the sensor node retrys if it doesnt receive the acknoweledgement. 

Acoustic propagation is greatly simplified to use straightline ranges for propagation delay and transmission losses. Probability of packet delivery is based on simulations of a modem waveform and receiver structure found in [Ultra-Low-Cost and Ultra-Low-Power, Miniature Acoustic Modems Using Multipath Tolerant Spread-Spectrum Techniques](https://doi.org/10.3390/electronics11091446). 



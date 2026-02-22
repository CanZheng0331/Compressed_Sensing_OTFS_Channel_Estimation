% Set system parameters
numDelayBins = 64;
numDopplerBins = 32;
carrierFreq = 4 * 10^9;
subcarrierSpacing = 15 * 10^3;
SNRdB = 0:5:30;
modType = 'QPSK'; % 'BPSK', '16QAM'
algorithmList = {'OMP'};
numExperiments = 1;

% Perform simulation and receive results
[nmse, timeInfo] = simulateOTFSChannelEstimation(numDelayBins, numDopplerBins, carrierFreq, subcarrierSpacing, SNRdB, modType, algorithmList, numExperiments);

% Print out data
disp('NMSE of each algorithm under different SNRs:');
disp(nmse);
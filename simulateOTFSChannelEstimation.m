function [nmseResults, timeResults] = simulateOTFSChannelEstimation(numDelayBins, numDopplerBins, carrierFreq, subcarrierSpacing, SNRdB, modType, algorithmList, numExperiments)
% simulateOTFSChannelEstimation: Performance simulation of the channel estimation algorithm in the OTFS system
%
% Input:
%   numDelayBins      - M
%   numDopplerBins    - N
%   carrierFreq       
%   subcarrierSpacing 
%   SNRdB             - 信噪比范围 (向量，例如 0:5:30)
%   modType           - 'BPSK', 'QPSK', '16QAM'
%   algorithmList     - {'OMP', ...}
%   numExperiments    - Monte Carlo simulation times
%
% Output:
%   nmseResults       - Normalized mean square error matrix (number of algorithms x number of SNR points)
%   timeResults       - Algorithm runtime matrix (number of algorithms x SNR points)

    %% --- 1. System parameter initialization ---                   
    fsamp = numDelayBins * subcarrierSpacing;               % Sampling Rate
    Meff = numDelayBins;                                    % 
    numSamps = numDelayBins * numDopplerBins;
    symbolDuration = (Meff / (numDelayBins * subcarrierSpacing)); % Symbol Duration Time (秒)
    
    % The modulator and symbol energy are dynamically set according to the input
    switch upper(modType)
        case 'BPSK'
            M_mod = 2; 
            modFunc = @(x) pskmod(x, 2, InputType="bit");
            Es = mean(abs(pskmod(0:1, 2)).^2);
        case 'QPSK'
            M_mod = 4; 
            modFunc = @(x) pskmod(x, 4, pi/4, InputType="bit");
            Es = mean(abs(pskmod(0:3, 4, pi/4)).^2);
        case '16QAM'
            M_mod = 16; 
            modFunc = @(x) qammod(x, 16, InputType="bit", UnitAveragePower=true);
            Es = 1; % UnitAveragePower already normalized
        otherwise
            error('The modulation method is not supported. Please select BPSK, QPSK, or 16QAM.');
    end
    M_bits = log2(M_mod);
    
    %% --- 2. Pilot and frame structure parameters ---
    centerDoppler = floor(numDopplerBins / 2);              % Doppler Center Index
    centerDelay = floor(numDelayBins / 2);                  % Delay Center Index

    pilotDelayLen = 6;                   % Pilot block length - delay dimension
    pilotDopplerLen = 10;                % Pilot block length - Doppler
    guardDelay = 5;                      % Delay domain guard interval
    guardDoppler = 9;                    % Doppler domain guard interval

    %% --- 3. Visualize pilot structure ---
    frameGrid = zeros(numDopplerBins, numDelayBins);
    % 0: Data, 1: Guard, 2: Pilot
    frameGrid(:, :) = 0; 
    frameGrid(numDopplerBins/2 - pilotDopplerLen/2 - guardDoppler + 1 : numDopplerBins/2 + pilotDopplerLen/2 + guardDoppler, ...
              numDelayBins/2 - pilotDelayLen/2 - guardDelay + 1 : numDelayBins/2 + pilotDelayLen/2 + guardDelay) = 1; 
    frameGrid(numDopplerBins/2 - pilotDopplerLen/2 + 1 : numDopplerBins/2 + pilotDopplerLen/2, ...
              numDelayBins/2 - pilotDelayLen/2 + 1 : numDelayBins/2 + pilotDelayLen/2) = 2; 

    figure(1);
    imagesc(frameGrid);
    colormap([0.6 0.9 0.6; 1 1 0.6; 1 0.4 0.2]); 
    title(sprintf('OTFS Frame Structure (%s, %dx%d)', modType, numDopplerBins, numDelayBins));
    xlabel('Delay Index');
    ylabel('Doppler Index');
    axis xy; 
    grid on;

    %% --- 4. Algorithm and simulation parameter configuration ---
    numAlgos = length(algorithmList);

    % Simplified channel parameters
    numTaps = 9;                         
    delayRange = [0, 2];                 
    dopplerRange = [-4, 3];              

    % Performance record matrix initialization
    nmseResults = zeros(numAlgos, length(SNRdB));
    timeResults = zeros(numAlgos, length(SNRdB));

    %% --- 5. Main loop ---
    for expNum = 1:numExperiments
        fprintf("Executing the %d/%d th Monte Carlo simulation...\n", expNum, numExperiments)
        
        % Generate random channel taps
        delayTaps = randi([delayRange(1), delayRange(2)], 1, numTaps);
        dopplerTaps = randi([dopplerRange(1), dopplerRange(2)], 1, numTaps);
        channelGains = (randn(1, numTaps) + 1j * randn(1, numTaps)) / sqrt(2);
        channelGains = channelGains / sqrt(sum(abs(channelGains).^2));
        
        chanParams.pathDelays      = delayTaps;
        chanParams.pathGains       = channelGains;                                     
        chanParams.pathDopplers    = dopplerTaps;
        chanParams.pathDopplerFreqs = chanParams.pathDopplers * 1/(numDopplerBins * symbolDuration);
        
        % Pilot sequence generation (Zadoff-Chu)
        nSeq = (0:pilotDelayLen*pilotDopplerLen-1)';                              
        zcSeq = exp(-1j*pi*1*nSeq.*(nSeq+1)/(pilotDelayLen*pilotDopplerLen));      
        pilotMatrix = reshape(zcSeq, pilotDopplerLen, pilotDelayLen);                 
        
        % Data symbol generation
        dataDD = zeros(numDopplerBins, numDelayBins); 
        dataBits = randi([0, 1], numDopplerBins * M_bits, numDelayBins);         
        % Dynamically calling the modulation function
        dataDD(1:numDopplerBins, :) = modFunc(dataBits); 
        
        % Mapping the guard interval and pilot signal in the DD domain
        dataDD(numDopplerBins/2 - pilotDopplerLen/2 - guardDoppler + 1 : numDopplerBins/2 + pilotDopplerLen/2 + guardDoppler, ...
               numDelayBins/2 - pilotDelayLen/2 - guardDelay + 1 : numDelayBins/2 + pilotDelayLen/2 + guardDelay) = 0;     
        dataDD(numDopplerBins/2 - pilotDopplerLen/2 + 1 : numDopplerBins/2 + pilotDopplerLen/2, ...
               numDelayBins/2 - pilotDelayLen/2 + 1 : numDelayBins/2 + pilotDelayLen/2) = pilotMatrix;             
                          
        % A temporary variable records the current loop error and time.
        tempNmse = zeros(numAlgos, length(SNRdB));
        tempTime = zeros(numAlgos, length(SNRdB));
        
        for snrIdx = 1:length(SNRdB)
            n0 = Es / (10^(SNRdB(snrIdx)/10));  % Noise
            noiseStdDev = sqrt(n0);
            %sigma_min2 = 1 * noiseStdDev; 
            
            % Real channel matrix H
            H_grid = zeros(numDopplerBins, numDelayBins);
            b_grid = zeros(numDopplerBins, numDelayBins);
            for i = 1:numTaps
                t_idx = chanParams.pathDelays(i) + centerDelay + 1;            
                f_idx = chanParams.pathDopplers(i) + centerDoppler + 1;          
                H_grid(f_idx, t_idx) = H_grid(f_idx, t_idx) + chanParams.pathGains(i);  
                b_grid(f_idx, t_idx) = 1;                              
            end
            
            % Generate the receive grid Y and the effective channel H_eff
            Y = zeros(numDopplerBins, numDelayBins); 
            Heff = zeros(numDopplerBins, numDelayBins);                     
            for k = 1:numDopplerBins
                for l = 1:numDelayBins
                    tempSum1 = 0;
                    for ki = 1:numDopplerBins
                        tempSum2 = 0;
                        for li = 1:numDelayBins                   
                            tempDelay = li - centerDelay;
                            tempDop = ki - centerDoppler;
                            tempH = H_grid(ki, li) * exp(-1j * 2 * pi * (tempDelay / (numDelayBins * subcarrierSpacing)) * (tempDop / (numDopplerBins * symbolDuration)));
                            Heff(ki, li) = tempH;
                            tempSum2 = tempSum2 + b_grid(ki, li) * tempH * dataDD(mod(k - tempDop, numDopplerBins) + 1, mod(l - tempDelay, numDelayBins) + 1);
                        end
                        tempSum1 = tempSum1 + tempSum2;
                    end
                    Y(k, l) = tempSum1;                 
                end
            end    
            % Additive noise
            noise = (randn(size(Y)) + 1j * randn(size(Y))) * sqrt(n0 / 2);
            Y = Y + noise;
        
            % Extract the received pilot signal yPilot
            yPilot = zeros(pilotDopplerLen * pilotDelayLen, 1);      
            yIndex = zeros(pilotDopplerLen * pilotDelayLen, 1);           
            
            idx = 1;
            for k_pie = -pilotDopplerLen/2 : pilotDopplerLen/2 - 1
                for l_pie = 0 : pilotDelayLen - 1
                    kReal = k_pie + centerDoppler + 1;                         
                    lReal = l_pie + centerDelay + 1;
                    tempIdx = l_pie * pilotDopplerLen + k_pie + pilotDopplerLen/2 + 1;
                    yIndex(idx) = tempIdx;
                    yPilot(tempIdx) = Y(kReal, lReal);
                    idx = idx + 1;
                end
            end        
            
            % Extract the actual generated equivalent channel hTrue   
            hTrue = zeros((2 * guardDoppler) * (guardDelay + 1), 1); 
            hIndex = zeros((2 * guardDoppler) * (guardDelay + 1), 1);        
            idx = 1;
            for kk = -guardDoppler : guardDoppler - 1
                for ll = 0 : guardDelay
                    kkReal = kk + centerDoppler + 1;                 
                    llReal = ll + centerDelay + 1;
                    tempHIdx  = (2 * ll + 1) * guardDoppler + kk + 1;
                    hIndex(idx) = tempHIdx;
                    hTrue(tempHIdx) = Heff(kkReal, llReal);
                    idx = idx + 1;
                end
            end        
            
            % Generate the sensing matrix Psi 
            Psi = zeros(pilotDopplerLen * pilotDelayLen, (2 * guardDoppler) * (guardDelay + 1));                
            fx = 1;
            for k_pie = -pilotDopplerLen/2 : pilotDopplerLen/2 - 1
                for l_pie = 0 : pilotDelayLen - 1
                    kReal = k_pie + centerDoppler + 1;                        
                    lReal = l_pie + centerDelay + 1;
                    fy = 1;
                    for kk = -guardDoppler : guardDoppler - 1
                        for ll = 0 : guardDelay
                            Psi(yIndex(fx), hIndex(fy)) = dataDD(kReal - kk, lReal - ll);
                            fy = fy + 1;
                        end
                    end
                    fx = fx + 1;
                end
            end      
            
            Psi_pinv = pinv(Psi);
            
            %% --- Algorithm Call Interface ---
            algoParams.taps = numTaps;
            algoParams.A_pinv = Psi_pinv;
            algoParams.hTrue = hTrue;
            
            truePower = sum(abs(hTrue).^2);
            
            for algIdx = 1:numAlgos
                algoName = algorithmList{algIdx};
                funcName = sprintf('Algorithm_%s', algoName);
                
                tStart = tic;
                % 核心调用
                hEst = feval(funcName, Psi, yPilot, algoParams);
                runTime = toc(tStart);
                
                tempTime(algIdx, snrIdx) = runTime;
                tempNmse(algIdx, snrIdx) = sum(abs(hTrue - hEst).^2) / truePower;
            end
        end
        
        nmseResults = nmseResults + tempNmse;
        timeResults = timeResults + tempTime;
    end

    %% --- 6. Data averaging and visualization ---
    nmseResults = nmseResults / numExperiments;
    timeResults = timeResults / numExperiments;

    nmse_dB = 10 * log10(nmseResults);

    lineStyles = {'>-', 'd-', '<-', 'o-', 's-'};

    figure(2);
    hold on;
    for algIdx = 1:numAlgos
        plot(SNRdB, nmse_dB(algIdx, :), lineStyles{mod(algIdx-1, length(lineStyles))+1}, 'LineWidth', 2);
    end
    title(sprintf('NMSE (dB) vs SNR (%s)', modType));
    xlabel('SNR (dB)');
    ylabel('NMSE (dB)');
    legend(algorithmList, 'Location', 'best');
    grid on;
    hold off;

    figure(3);
    hold on;
    for algIdx = 1:numAlgos
        plot(SNRdB, timeResults(algIdx, :), lineStyles{mod(algIdx-1, length(lineStyles))+1}, 'LineWidth', 2);
    end
    title('Running Time vs SNR');
    xlabel('SNR (dB)');
    ylabel('Running Time (s)');
    legend(algorithmList, 'Location', 'best');
    grid on;
    hold off;
end
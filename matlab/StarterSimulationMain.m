% StarterSimulationMain   A very simple simulation to get started.
%    All in one matlab file. 
%    Providing node locations and configuring the acoustic properties of
%    the channel we can derive the propagation delay and packet success
%    probabilities between two nodes. Knowing these measn we can simulate
%    different scenarios in a relatively straightforward way.
%
%   Benjamin Sherlock, Newcastle University, 2023


% Clean the workspace
clear all;
close all;


% Configuration
% ========================================================================

% Acoustic Properties
% Attenuation Alpha (dB/km at Carrier Frequency)
AttenuationAlpha = 6.0;

% Speed of sound (m/s)
SpeedOfSound = 1520.0;

% Ambient noise spectral density in dB re 1uPa^2/Hz (See Wenz Curves at Carrier Frequency)
AmbientNoiseSpectralDensity = 50.0;

% Modem Properties
% Modem Source Level (dB re 1uPa @ 1m)
ModemSourceLevel = 168.0;

% Modem Bandwidth (Hz)
ModemBandwidth = 8000.0;


% Node position: [ [x, y, depth]; ... ]
NodePositions = [ ...
    [0.0, 0.0, 10.0]; ... 
    [1800.0, 0.0, 10.0];
    [0.0, 2000.0, 10.0]];

% Calculate the Node Count
NodeCount = length(NodePositions(:,1));

% Node label: ['Name'; 'Name'; ... ]
NodeLabels = ['A'; 'B'; 'C'];

% Node addresses: [ 1; 2; 3; ... ]
NodeAddresses = [1; 2; 3];


% Generate Lookup Tables
% ========================================================================
NodeRanges = CalculateNodeRanges(NodePositions);
NodePropagationDelays = NodeRanges ./ SpeedOfSound;
NodeDeliveryProbabilities = CalculateNodeDeliveryProbabilities(NodePositions, ...
    AttenuationAlpha, AmbientNoiseSpectralDensity, ModemBandwidth, ModemSourceLevel);

% Run
% ========================================================================

% Print out information
PrintNodeDeliveryProbabilities(NodeDeliveryProbabilities, NodeLabels);

% Display the scene
DisplayNodes(NodePositions, NodeLabels);

% Run a simulation
RunSimulationA(NodePropagationDelays, NodeDeliveryProbabilities, NodeLabels, NodeAddresses, 1);


% Simulation Functions
% ========================================================================

function RunSimulationA(NodePropagationDelays, NodeDeliveryProbabilities, NodeLabels, NodeAddresses, IterationCount)
    % RunSimulationA - This simulation runs a scenario where each sensor
    % node periodically transmits data to the gateway node (1). 

    RunDuration = 1.0*60.0*60.0; % seconds (1hours)
    StepFrequency = 16e3; % 16kHz update rate
    StepPeriod = 1 / StepFrequency; % seconds between each step.
    RunStepCount = RunDuration/StepPeriod;


    SensorTransmitPeriodTime = 1.0*60.0*60.0; % seconds (1 hour)
    SensorTransmitPeriodStepCount = floor(SensorTransmitPeriodTime / StepPeriod); 

    NodeCount = length(NodeDeliveryProbabilities(:,1));
    NodeTransmitTimes = [0:NodeCount-1] * (SensorTransmitPeriodTime / NodeCount); % seconds
    NodeTransmitStepCounts = floor(NodeTransmitTimes ./ StepPeriod);
    
    ModemBitRate = 100.0; % Arbitrary data rate of 100bps 
    
    % Each Channels is unidirectional (From, To) and has a queue of
    % Packets currently travelling through the water.
    Channels = cell(NodeCount, NodeCount); 
    for FromIdx = 1:NodeCount
        for ToIdx = 1:NodeCount
            Channels{FromIdx, ToIdx} = {}; % Empty Queue of Transmitted Packets
        end
    end
    ChannelsCount = NodeCount*NodeCount;


    % Create the queues for transmitted packets
    TransmittedPackets = cell(NodeCount);
    for NodeIdx = 1:NodeCount
        TransmittedPackets{NodeIdx} = {}; % Empty Queue of Transmitted Packets
    end

    % Create the queues for received packets
    ReceivedPackets = cell(NodeCount);
    for NodeIdx = 1:NodeCount
        ReceivedPackets{NodeIdx} = {}; % Empty Queue of Received Packets
    end


    for IterationIdx = 1:IterationCount
        
        % Optimsations to speed up the simulation
        ChannelsEmpty = true;
        TransmittedPacketsEmpty = true;
        ReceivedPacketsEmpty = true;

        ProcessStartTime = tic;
        for StepIdx = 1:RunStepCount
            CurrentTime = (StepIdx-1) * StepPeriod;

            % Update printout
            if mod((StepIdx-1), StepFrequency*600) == 0
                ProcessTime = toc(ProcessStartTime);
                fprintf('%09.3fs :Tick: Elapsed Time Since Last Tick=%0.3fs \n', CurrentTime, ProcessTime);
                ProcessStartTime = tic;
            end

            % User Scenario: Nodes create transmissions at certain times
            for NodeIdx = 1:NodeCount
                
                %if mod(StepIdx-1, NodeTransmitStepCounts(NodeIdx)) == 0
                if StepIdx-1 == NodeTransmitStepCounts(NodeIdx)
                    % Transmit Sensor Data

                    SourceAddress = NodeAddresses(NodeIdx);
                    DestinationAddress = NodeAddresses(1);
                    PayloadBytes = ['H', 'e', 'l', 'l', 'o'];                    
                    TransmitStartTime = CurrentTime;
                    TransmitDuration = (length(PayloadBytes) * 8) / ModemBitRate;

                    % Construct a Packet struct
                    TxPacket = Packet(SourceAddress, DestinationAddress, PayloadBytes, TransmitStartTime, TransmitDuration);
                    
                    % The modem is omni directional so all outgoing 
                    % channels from this node will potentially deliver the
                    % packet.
                    
                    TransmittedPackets{NodeIdx}{length(TransmittedPackets{NodeIdx})+1} = TxPacket;
                    TransmittedPacketsEmpty = false;

                    fprintf(['%09.3fs :Tx: ' 'Node ' NodeLabels(NodeIdx) ' has transmitted a packet. '], CurrentTime);
                    PrintPacket(TxPacket);
                    fprintf('\n');
                end
            end

            if ~TransmittedPacketsEmpty
                % Outgoing Transmissions from each of the modems into the channels
                for FromIdx = 1:NodeCount
                    while ~isempty(TransmittedPackets{FromIdx})
                        % The From Node has a packet to send
                        TxPacket = TransmittedPackets{FromIdx}{1};
                        % Delete from queue
                        TransmittedPackets{FromIdx}(1) = [];
    
                        for ToIdx = 1:NodeCount
                            % Use the probability of packet success and randmise
                            % whether this packet will arrive at the destination.
                            ProbabilityOfPacketDelivery = NodeDeliveryProbabilities(FromIdx, ToIdx);
                            if rand(1) < ProbabilityOfPacketDelivery
                                % Only if the random number (0.0->1.0) is less than
                                % the probability of delivery do we send the
                                % packet and append to the Channel queue.
                                Channels{FromIdx, ToIdx}{length(Channels{FromIdx, ToIdx})+1} = TxPacket;
    
                                ChannelsEmpty = false;
                            end
                        end
                    end
                    TransmittedPacketsEmpty = true;
                end
            end


            % Only process the propagation if there are packets in the
            % channels. Thsi helps to speed up the simulation.
            if ~ChannelsEmpty
                ChannelsEmpty = true;  % Set the flag but we will check each individual channel anyway.

                % Propagation of packets through the channels to the receiving modems
                for FromIdx = 1:NodeCount
                    for ToIdx = 1:NodeCount
    
                        % If there are packets in the channel, check to see if
                        % they have now arrived at the destination.
                        if ~isempty(Channels{FromIdx, ToIdx})
                            ChannelsEmpty = false;  % Not empty so clear the flag

                            if (CurrentTime >= (Channels{FromIdx, ToIdx}{1}.TransmitStartTime ... 
                                + Channels{FromIdx, ToIdx}{1}.TransmitDuration ...
                                + NodePropagationDelays(FromIdx, ToIdx) ) )   
    
                                % Get packet at the head of the queue
                                RxPacket = Channels{FromIdx, ToIdx}{1};
                                % Delete from queue
                                Channels{FromIdx, ToIdx}(1) = [];
        
                                % Put the packet into received packet queue for the
                                % destination node to process in the next stage.
                                ReceivedPackets{ToIdx}{length(ReceivedPackets{ToIdx})+1} = RxPacket;

                                ReceivedPacketsEmpty = false;
                            end
                        end
                    end
                end

            end


            if ~ReceivedPacketsEmpty
                % 3. Nodes act on received packets
                for NodeIdx = 1:NodeCount
                    while ~isempty(ReceivedPackets{NodeIdx})
                        % This node has just received some packets
        
                        % Get the Packet at the head of the queue
                        RxPacket = ReceivedPackets{NodeIdx}{1};
                        % Delete from queue
                        ReceivedPackets{NodeIdx}(1) = [];
    
                        fprintf(['%09.3fs :Rx: ' 'Node ' NodeLabels(NodeIdx) ' has received a packet. '], CurrentTime);
                        PrintPacket(RxPacket);
                        fprintf('\n');
    
                    end
                    ReceivedPacketsEmpty = true;
                end
            end
        
        end

    end
end




% Calculation Functions
% ========================================================================

function PropagationDelay = CalculatePropagationDelay(Range, SpeedOfSound)
    % CalculatePropagationDelay - Use straight line range and speed of
    % sound to calculate the propagation delay of an acoustic signal.
    PropagationDelay = Range / SpeedOfSound;
end

function TransmissionLosses = CalculateTransmissionLosses(Range, AttenuationAlpha)
    % CalculateTransmissionLosses - Use the range and attenutation
    % coefficient to calculate the dB transmission losses based on 
    % spherical spreading.
    TransmissionLosses = (20.0 * log10(Range)) + (AttenuationAlpha * Range * 0.001);
end

function ReceivedSnr = CalculateReceivedSnr(Range, AttenuationAlpha, ...
    AmbientNoiseSpectralDensity, ModemBandwidth, ModemSourceLevel)
    % CalculateReceivedSnr - Uses the source level and transmission losses
    % combined with the ambient noise levels to calculate an effective SNR
    % at the receiver.

    % Calculate the Transmission losses
    TransmissionLosses = CalculateTransmissionLosses(Range, AttenuationAlpha);

    % Calculate the Noise losses
    NoiseLosses = AmbientNoiseSpectralDensity + (10.0 * log10(ModemBandwidth));

    % Combine the Transmission and Noise losses
    TotalLosses = TransmissionLosses + NoiseLosses;

    % Calculate the effective received SNR (dB)
    ReceivedSnr = ModemSourceLevel - TotalLosses;
end

function ProbabilityOfDelivery = CalculateProbabilityOfDelivery(SNR)
    % CalculateProbabilityOfDelivery - For a given Receiver SNR calculate
    % the probability that a packet will be successfully delivered.
    % This example uses a simple look up table that was generated from
    % simulations of a modem receiver structure. 

    PacketSnrValues = -9:9;
    %PacketErrorRates = [1.000, 1.000, 1.000, 1.000, 0.998, ...
    %                    0.897, 0.388, 0.080, 0.013, 0.002, ...
    %                    0.000, 0.000, 0.000, 0.000, 0.000, ...
    %                    0.000, 0.000, 0.000, 0.000];
    PacketErrorRates = [1.000, 1.000, 1.000, 1.000, 1.000, ...
                        1.000, 0.997, 0.952, 0.804, 0.553, ...
                        0.245, 0.089, 0.020, 0.010, 0.003, ...
                        0.001, 0.000, 0.000, 0.000];



    if SNR < PacketSnrValues(1)
        % SNR is less than our lowest value so failed delivery
        ProbabilityOfDelivery = 0.0;
    elseif SNR >= PacketSnrValues(end)
        % SNR is greater than our highest value so successful delivery
        ProbabilityOfDelivery = 1.0;
    else
        for idx = 1:length(PacketSnrValues)
            if SNR <= PacketSnrValues(idx)
                % Convert Probability of Error to Probability of Success
                ProbabilityOfDelivery = 1.0 - PacketErrorRates(idx);
                break;
            end
        end
    end
end

function NodeRanges = CalculateNodeRanges(NodePositions)
    % CalculateNodeRanges - Calculate a Range between each pair of nodes
    NodeCount = length(NodePositions(:,1));

    NodeRanges = zeros(NodeCount, NodeCount);

    for FromIdx = 1:NodeCount
        for ToIdx = 1:NodeCount
            Range = sqrt( ...
                    (NodePositions(ToIdx, 1) - NodePositions(FromIdx, 1)).^2 ...
                  + (NodePositions(ToIdx, 2) - NodePositions(FromIdx, 2)).^2 ...
                  + (NodePositions(ToIdx, 3) - NodePositions(FromIdx, 3)).^2 ...
                );

            NodeRanges(FromIdx, ToIdx) = Range;
        end
    end

end

function NodeDeliveryProbabilities = CalculateNodeDeliveryProbabilities(NodePositions, ...
    AttenuationAlpha, AmbientNoiseSpectralDensity, ModemBandwidth, ModemSourceLevel)
    % CalculateNodeDeliveryProbabilities - Calculate the probability of
    % packet delivery between every pair of nodes for the given positions
    % and acoustic conditions.
    NodeCount = length(NodePositions(:,1));

    NodeDeliveryProbabilities = zeros(NodeCount, NodeCount);
    
    for FromIdx = 1:NodeCount
        for ToIdx = 1:NodeCount
            ProbabilityOfDelivery = 0.0;
            
            if FromIdx ~= ToIdx  % Cannot transmit to itself                

                Range = sqrt( ...
                        (NodePositions(ToIdx, 1) - NodePositions(FromIdx, 1)).^2 ...
                      + (NodePositions(ToIdx, 2) - NodePositions(FromIdx, 2)).^2 ...
                      + (NodePositions(ToIdx, 3) - NodePositions(FromIdx, 3)).^2 ...
                    );
    
                ReceivedSnr = CalculateReceivedSnr(Range, AttenuationAlpha, ...
                    AmbientNoiseSpectralDensity, ModemBandwidth, ModemSourceLevel);
    
                ProbabilityOfDelivery = CalculateProbabilityOfDelivery(ReceivedSnr);
            end

            NodeDeliveryProbabilities(FromIdx, ToIdx) = ProbabilityOfDelivery;
        end
    end


end

% Printout Functions
% ========================================================================
function PrintNodeDeliveryProbabilities(NodeDeliveryProbabilities, NodeLabels)
    NodeCount = length(NodeDeliveryProbabilities(:,1));
    
    fprintf('Probability of Packet Delivery Success\n')
    fprintf('\tTo:  \t')
    for ToIdx = 1:NodeCount
        fprintf([NodeLabels(ToIdx) '\t\t']);
    end
    fprintf('\n');

    for FromIdx = 1:NodeCount
        fprintf('From:\t');
        fprintf([NodeLabels(FromIdx) '\t']);
        for ToIdx = 1:NodeCount
            fprintf('%0.3f\t', NodeDeliveryProbabilities(FromIdx, ToIdx));
        end
        fprintf('\n');
    end

end

function PrintPacket(Packet)
    % PrintPacket - Print out the Packet contents
    fprintf('Packet:');

    fprintf(' SourceAddress=%d', Packet.SourceAddress);
    fprintf(' DestinationAddress=%d', Packet.DestinationAddress);
    fprintf(' PayloadBytes=[');
    for idx = 1:length(Packet.PayloadBytes)

        fprintf('%c', Packet.PayloadBytes(idx));

        if idx < length(Packet.PayloadBytes)
            fprintf(',');
        end
    end
    fprintf(']');

    fprintf(' TransmitStartTime=%0.3f', Packet.TransmitStartTime);
    fprintf(' TransmitDuration=%0.3f', Packet.TransmitDuration);
end


% Display Functions
% ========================================================================

function DisplayNodes(NodePositions, NodeLabels)
    % DisplayNodes - Show the locations of the nodes in 2D/3D space.

    LimitsMargin = 10.0;
    Xlimits = [min(NodePositions(:,1))-LimitsMargin, max(NodePositions(:,1))+LimitsMargin];
    Ylimits = [min(NodePositions(:,2))-LimitsMargin, max(NodePositions(:,2))+LimitsMargin];
    Zlimits = [0.0, max(NodePositions(:,3))+LimitsMargin];

    % 3D view
    fig = figure;
    scatter3(NodePositions(:,1), NodePositions(:,2), NodePositions(:,3), 'filled');
    hold on;
    text(NodePositions(:,1), NodePositions(:,2), NodePositions(:,3), NodeLabels(:));
    hold off;

    xlim(Xlimits);
    ylim(Ylimits);
    zlim(Zlimits);

    title('Node Positions');

    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Depth (m)');

    % Reverse direction of Z axis so that increasing depth goes downwards.
    fig.CurrentAxes.ZDir = 'Reverse';  

end



% Container Structs
% ========================================================================

function Packet = Packet(SourceAddress, DestinationAddress, PayloadBytes, TransmitStartTime, TransmitDuration)
    % Packet - Constructs a struct of a Packet containing the provided
    % information. Access the fields on the returned struct such as
    % Packet.SourceAddress .

    Packet = struct(...
        'SourceAddress', SourceAddress, ...
        'DestinationAddress', DestinationAddress, ...
        'PayloadBytes', PayloadBytes, ...
        'TransmitStartTime', TransmitStartTime, ...
        'TransmitDuration', TransmitDuration ...
        );
end

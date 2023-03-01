% StarterSimulationMain   A very simple simulation to get started.
%    All in one matlab file. 
%    Providing node locations and configuring the acoustic properties of
%    the channel we can derive the propagation delay and packet success
%    probabilities between two nodes. Knowing these means we can simulate
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
%NodePositions = [ ...
%        [0.0, 0.0, 10.0]; ... 
%        [900.0, 500.0, 10.0]; ...
%        [1100.0, -900.0, 10.0]; ...
%        [1900.0, -50.0, 10.0]; ...
%        [2400.0, -2000.0, 10.0] ...
%    ];

NodePositions = [ ...
        [0.0, 0.0, 10.0]; ... 
        [1900.0, 800.0, 30.0]; ...
        [1400.0, -1000.0, 30.0]; ...
    ];

% Calculate the Node Count
NodeCount = length(NodePositions(:,1));

% Node label: ['Name'; 'Name'; ... ]
%NodeLabels = ['A'; 'B'; 'C'; 'D'; 'E'];
NodeLabels = ['A'; 'B'; 'C'];

% Node addresses: [ 1; 2; 3; ... ]
%NodeAddresses = [1; 2; 3; 4; 5];
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
PrintNodeRanges(NodeRanges, NodeLabels);
PrintNodeDeliveryProbabilities(NodeDeliveryProbabilities, NodeLabels);

% Display the scene
DisplayNodes(NodePositions, NodeLabels, NodeDeliveryProbabilities);
drawnow;

% Run a simulation
RunSimulationA(NodePositions, NodePropagationDelays, NodeDeliveryProbabilities, NodeLabels, NodeAddresses, SpeedOfSound, 1);


% Simulation Functions
% ========================================================================

function RunSimulationA(NodePositions, NodePropagationDelays, NodeDeliveryProbabilities, NodeLabels, NodeAddresses, SpeedOfSound, IterationCount)
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

    % Log to record the state of each node at each time step.
    NodeStateLogs = zeros(NodeIdx, RunStepCount);
    NodeStateListening = 0;
    NodeStateReceiving = 1;
    NodeStateReceivedSuccessfully = 2;
    NodeStateTransmitting = 3;


    % Scenario - Acks
    NodeAwaitingAcks = zeros(NodeCount, 1);
    NodePacketsToResend = cell(NodeCount);
    NodeRetriesCount = zeros(NodeCount, 1);
    TimeBetweenRetries = 5.0; % (s)
    MaxRetries = 3;

        
    % Optimsations to speed up the simulation
    ChannelsEmpty = true;
    TransmittedPacketsEmpty = true;
    ReceivedPacketsEmpty = true;

    % Visualisation
    AnimationFigHandle = CreateAcousticPropagationFigure(NodePositions, NodeLabels);
    AnimationUpdated = false;

    % Simulator Mainloop
    ProcessStartTime = tic;
    for StepIdx = 1:RunStepCount
        CurrentTime = (StepIdx-1) * StepPeriod;

        % Update printout
        if mod((StepIdx-1), StepFrequency*600) == 0
            ProcessTime = toc(ProcessStartTime);
            fprintf('%09.3fs :Tick: Elapsed Time Since Last Tick=%0.3fs \n', CurrentTime, ProcessTime);
            ProcessStartTime = tic;
        end

        % Visualisation printout
        if mod((StepIdx-1), StepFrequency*0.1) == 0
            if ~ChannelsEmpty
                UpdateAcousticPropagationFigure(AnimationFigHandle, NodePositions, NodeLabels, Channels, SpeedOfSound, CurrentTime)
                drawnow;
                AnimationUpdated = false;
            elseif ~AnimationUpdated
                UpdateAcousticPropagationFigure(AnimationFigHandle, NodePositions, NodeLabels, Channels, SpeedOfSound, CurrentTime)
                drawnow;
                AnimationUpdated = true;
            end

        end

        % User Code: Act on Timestamp
        % Nodes create transmissions at certain times
        for NodeIdx = 1:NodeCount
            
            if NodeIdx > 1 % Sensor Nodes Only
                %if mod(StepIdx-1, NodeTransmitStepCounts(NodeIdx)) == 0
                if StepIdx-1 == NodeTransmitStepCounts(NodeIdx)
                    % Transmit New Sensor Data
    
                    SourceAddress = NodeAddresses(NodeIdx);
                    DestinationAddress = NodeAddresses(1);
                    PayloadBytes = ['S', 'e', 'n', 's', 'o', 'r', ' ', 'D', 'a', 't', 'a'];
                    TransmitStartTime = CurrentTime;
                    TransmitDuration = (length(PayloadBytes) * 8) / ModemBitRate;
    
                    % Construct a Packet struct
                    TxPacket = Packet(SourceAddress, DestinationAddress, PayloadBytes, TransmitStartTime, TransmitDuration);
                    
                    % Put the Transmit packet in the queue for the simulator
                    % core to handle.
                    TransmittedPackets{NodeIdx}{length(TransmittedPackets{NodeIdx})+1} = TxPacket;
                    TransmittedPacketsEmpty = false;
    
                    % Store information for handling Acks/Retries
                    NodeAwaitingAcks(NodeIdx) = 1;
                    NodePacketsToResend{NodeIdx} = TxPacket;
                    NodeRetriesCount(NodeIdx) = 0;
       
                elseif NodeAwaitingAcks(NodeIdx)
                    % Check for Timeout
                    if CurrentTime >= NodePacketsToResend{NodeIdx}.TransmitStartTime + TimeBetweenRetries
                        % Check if we have any retries left for this packet
                        if NodeRetriesCount(NodeIdx) < MaxRetries
                            % Resend
                            TxPacket = NodePacketsToResend{NodeIdx};
                            TxPacket.TransmitStartTime = CurrentTime; % Update the Transmit Time Information
    
                            % Put the Transmit packet in the queue for the simulator
                            % core to handle.
                            TransmittedPackets{NodeIdx}{length(TransmittedPackets{NodeIdx})+1} = TxPacket;
                            TransmittedPacketsEmpty = false;
    
                            % Increment the retry count
                            NodeRetriesCount(NodeIdx) = NodeRetriesCount(NodeIdx) + 1;
                            NodePacketsToResend{NodeIdx} = TxPacket;
                        else
                            % We have sent it multiple times already.
                            NodeAwaitingAcks(NodeIdx) = 0;
                            NodePacketsToResend{NodeIdx} = [];
                            NodeRetriesCount(NodeIdx) = 0;
                        end
                    end
    
                end
            end
        end


        % Simulator Core:
        if ~TransmittedPacketsEmpty
            % Outgoing Transmissions from each of the modems into the channels
            for FromIdx = 1:NodeCount
                while ~isempty(TransmittedPackets{FromIdx})
                    % The From Node has a packet to send
                    % Get from the head of the queue
                    TxPacket = TransmittedPackets{FromIdx}{1};
                    % Delete from queue
                    TransmittedPackets{FromIdx}(1) = [];

                    % Send the packet to all outgoing channels from
                    % this node
                    for ToIdx = 1:NodeCount
                        if ToIdx ~= FromIdx  % Don't transmit to itself
                            % Send the packet and append to the Channel queue.
                            Channels{FromIdx, ToIdx}{length(Channels{FromIdx, ToIdx})+1} = TxPacket;
                            ChannelsEmpty = false;                           
                        end
                    end

                    % Recording the transmitting state                    
                    NodeStateLogs(FromIdx, StepIdx:StepIdx+floor(TxPacket.TransmitDuration*StepFrequency)) = NodeStateTransmitting;

                    fprintf(['%09.3fs :Tx: ' 'Node ' NodeLabels(FromIdx) ' has transmitted a packet. '], CurrentTime);
                    PrintPacket(TxPacket);
                    fprintf('\n');
                end
                TransmittedPacketsEmpty = true;
            end
        end


        % Simulator Core: 
        % Only process the propagation if there are packets in the
        % channels. This helps to speed up the simulation.
        if ~ChannelsEmpty
            ChannelsEmpty = true;  % Set the flag but we will check each individual channel anyway.

            % Propagation of packets through the channels to the receiving modems
            for FromIdx = 1:NodeCount
                for ToIdx = 1:NodeCount

                    % If there are packets in the channel, check to see if
                    % they have now arrived at the destination.
                    if ~isempty(Channels{FromIdx, ToIdx})
                        ChannelsEmpty = false;  % Not empty so clear the flag

                        % CurrentTime >= TransmitStartTime + TransmitDuration + PropagationDelay
                        if (CurrentTime >= (Channels{FromIdx, ToIdx}{1}.TransmitStartTime ... 
                            + Channels{FromIdx, ToIdx}{1}.TransmitDuration ...
                            + NodePropagationDelays(FromIdx, ToIdx) ) )   

                            % Get packet at the head of the queue
                            RxPacket = Channels{FromIdx, ToIdx}{1};
                            % Delete from queue
                            Channels{FromIdx, ToIdx}(1) = [];

                            % Recording the receiving state
                            NodeStateLogs(ToIdx, StepIdx-floor(RxPacket.TransmitDuration*StepFrequency):StepIdx) = NodeStateReceiving;
    
                            % Use the probability of packet success and randmise
                            % whether this packet will arrive at the destination.
                            ProbabilityOfPacketDelivery = NodeDeliveryProbabilities(FromIdx, ToIdx);
                            if rand(1) < ProbabilityOfPacketDelivery
                                % Only if the random number (0.0->1.0) is less than
                                % the probability of delivery did the
                                % packet arrive.

                                % Put the packet into received packet queue for the
                                % destination node to process in the next stage.
                                ReceivedPackets{ToIdx}{length(ReceivedPackets{ToIdx})+1} = RxPacket;                                
                                ReceivedPacketsEmpty = false;

                                fprintf(['%09.3fs :Rx: ' 'Node ' NodeLabels(ToIdx) ' has received a packet. '], CurrentTime);
                                PrintPacket(RxPacket);
                                fprintf('\n');

                                % Recording the received successfully state
                                NodeStateLogs(ToIdx, StepIdx) = NodeStateReceivedSuccessfully;
                            end

                            
                        end
                    end
                end
            end
        end


        % User Code: Act on Received Packets
        if ~ReceivedPacketsEmpty
            
            for NodeIdx = 1:NodeCount
                while ~isempty(ReceivedPackets{NodeIdx})
                    % This node has just received some packets
    
                    % Get the Packet at the head of the queue
                    RxPacket = ReceivedPackets{NodeIdx}{1};
                    % Delete from queue
                    ReceivedPackets{NodeIdx}(1) = [];

                    % Process the packet here. 
                    if RxPacket.DestinationAddress == NodeAddresses(NodeIdx)
                        % This packet is addressed to the current node
                        if NodeIdx == 1 % This is the Gateway Node
                            % Sensor Data Received
                            fprintf('%09.3fs :Gateway: Has Received Sensor Data From %d\n', CurrentTime, RxPacket.SourceAddress);
                            
                            % Send an Acknowledgment
                            SourceAddress = NodeAddresses(NodeIdx);
                            DestinationAddress = RxPacket.SourceAddress;
                            PayloadBytes = ['A', 'c', 'k'];                    
                            TransmitStartTime = CurrentTime;
                            TransmitDuration = (length(PayloadBytes) * 8) / ModemBitRate;
            
                            % Construct a Packet struct
                            TxPacket = Packet(SourceAddress, DestinationAddress, PayloadBytes, TransmitStartTime, TransmitDuration);
                            
                            % Put the Transmit packet in the queue for the simulator
                            % core to handle.
                            TransmittedPackets{NodeIdx}{length(TransmittedPackets{NodeIdx})+1} = TxPacket;
                            TransmittedPacketsEmpty = false;

                            fprintf('%09.3fs :Gateway: Sending Acknowledgment To %d\n', CurrentTime, TxPacket.DestinationAddress);

                        else % This is a Sensor Node
                            % Acknowledgement Received
                            if RxPacket.SourceAddress == NodeAddresses(1)
                                fprintf('%09.3fs :Sensor: Has Received Acknowledgment From %d\n', CurrentTime, RxPacket.SourceAddress);

                                % Clear the retry information
                                NodeAwaitingAcks(NodeIdx) = 0;
                                NodePacketsToResend{NodeIdx} = [];
                                NodeRetriesCount(NodeIdx) = 0;
                            end
                        end
                    end

                end
                ReceivedPacketsEmpty = true;
            end
        end
    
    end



    DisplayNodeStateLogs(NodeStateLogs, NodeLabels, StepPeriod);

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
function PrintNodeRanges(NodeRanges, NodeLabels)
    NodeCount = length(NodeRanges(:,1));
    
    fprintf('Ranges (m)\n')
    fprintf('\tTo:  \t')
    for ToIdx = 1:NodeCount
        fprintf([NodeLabels(ToIdx) '\t\t\t']);
    end
    fprintf('\n');

    for FromIdx = 1:NodeCount
        fprintf('From:\t');
        fprintf([NodeLabels(FromIdx) '\t']);
        for ToIdx = 1:NodeCount
            fprintf('%09.3f\t', NodeRanges(FromIdx, ToIdx));
        end
        fprintf('\n');
    end

end

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

function DisplayNodes(NodePositions, NodeLabels, NodeDeliveryProbabilities)
    % DisplayNodes - Show the locations of the nodes in 2D/3D space.

    LimitsMargin = 50.0;
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

    if nargin > 2
        % Show Channel Probabilities
        hold on;
        NodeCount = length(NodePositions(:,1));
        for FromIdx = 1:NodeCount
            for ToIdx = 1:NodeCount
                Label = sprintf('%0.3f', NodeDeliveryProbabilities(FromIdx, ToIdx));
                DrawArrow(NodePositions(FromIdx, :), NodePositions(ToIdx, :));
            end
        end
        hold off;
    end

end

function DrawArrow(NodePositionA, NodePositionB)
    % DrawArrow - Add an arrow on the plot from A to B. 
    quiver3(NodePositionA(1), NodePositionA(2), NodePositionA(3), ...
        NodePositionB(1)-NodePositionA(1), NodePositionB(2)-NodePositionA(2), NodePositionB(3)-NodePositionA(3),...
        'filled', 'LineWidth', 1.0, 'AutoScale', 'off');
end

function DisplayNodeStateLogs(NodeStateLogs, NodeLabels, StepPeriod)
    % DisplayNodeStateLogs - Display a timeline showing the state of each
    % node

    NodeCount = length(NodeStateLogs(:,1));
    
    NumRows = NodeCount;
    NumCols = 1;
    PlotCount = 1;
    figure('NumberTitle','off', 'Name','Node State Logs', 'Position',[50 50 1200 800])
    sgtitle('Node State Logs');
    axes = [];

    t = StepPeriod*(0:length(NodeStateLogs(1,:))-1);
    
    for NodeIdx = 1:NodeCount
        ax = subplot(NumRows, NumCols, PlotCount);   
        axes = [axes ax];
        %plot(t(NodeStateLogs(NodeIdx,:) > 0), NodeStateLogs(NodeIdx, NodeStateLogs(NodeIdx,:) > 0 ), 'LineStyle','none', 'Marker','.', 'Color','r');
       
        stem(t(NodeStateLogs(NodeIdx,:) == 1), NodeStateLogs(NodeIdx, NodeStateLogs(NodeIdx,:) == 1), 'Color','#CCCCCC', 'Marker','none');
        hold on;        
        stem(t(NodeStateLogs(NodeIdx,:) == 3), NodeStateLogs(NodeIdx, NodeStateLogs(NodeIdx,:) == 3), 'Color','#999999', 'Marker','none');
        stem(t(NodeStateLogs(NodeIdx,:) == 2), NodeStateLogs(NodeIdx, NodeStateLogs(NodeIdx,:) == 2), 'Color','#222222', 'LineWidth',1.5, 'Marker','x', 'MarkerSize',10);
        hold off;
        title(['Node ' NodeLabels(NodeIdx)]);
        xlabel('Time (s)');
        ylabel('Node State');
        ylim([-0.1 3.1]);

        yticks([0, 1, 2, 3]);
        yticklabels({'Listening', 'Receiving', 'Received OK', 'Transmitting'});
        grid on;
           
        
        PlotCount = PlotCount + 1;
    end

    linkaxes(axes,'x');
    axes(1).XLim = [-10 max(t)+10];

end


function FigHandle = CreateAcousticPropagationFigure(NodePositions, NodeLabels)
    % CreateAcousticPropagationFigure

    FigHandle = figure;
    % Set size
    FigHandle.Position = [100 100 800 600];

    LimitsMargin = 500.0;
    Xlimits = [min(NodePositions(:,1))-LimitsMargin, max(NodePositions(:,1))+LimitsMargin];
    Ylimits = [min(NodePositions(:,2))-LimitsMargin, max(NodePositions(:,2))+LimitsMargin];
    Zlimits = [0.0, max(NodePositions(:,3))+LimitsMargin];

    % 3D view
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
    FigHandle.CurrentAxes.ZDir = 'Reverse';  

    %view(0, 90);
    view(2);

    axis square;

end

function UpdateAcousticPropagationFigure(FigHandle, NodePositions, NodeLabels, Channels, SpeedOfSound, CurrentTime)
    % UpdateAcousticPropagationFigure - Update the existing figure with current position of acoustic
    % propagations.

    if isempty(FigHandle)
        FigHandle = figure; 
    else
        FigHandle = figure(FigHandle); 

        % If we're updating an existing figure then only delete the
        % axes lines. Leaving the rest intact.
        axes = findobj(FigHandle, 'type', 'axes', '-not', 'tag', 'legend', '-not', 'tag', 'Colorbar');
        if ~isempty(axes)                
            delete(axes.Children);
        end
    end

    LimitsMargin = 500.0;
    Xlimits = [min(NodePositions(:,1))-LimitsMargin, max(NodePositions(:,1))+LimitsMargin];
    Ylimits = [min(NodePositions(:,2))-LimitsMargin, max(NodePositions(:,2))+LimitsMargin];
    Zlimits = [0.0, max(NodePositions(:,3))+LimitsMargin];

    % 3D view
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
    FigHandle.CurrentAxes.ZDir = 'Reverse';  

    hold on;
    NodeCount = length(NodePositions(:,1));
    for FromIdx = 1:NodeCount
        % Centre = NodePositions(FromIdx,1), NodePositions(FromIdx,2), NodePositions(FromIdx,3)
        TxX = NodePositions(FromIdx,1);
        TxY = NodePositions(FromIdx,2);
        TxZ = NodePositions(FromIdx,3);

        for ToIdx = 1:NodeCount
            if ~isempty(Channels{FromIdx, ToIdx})
                % Packets in transit so display the acoustic waves.
                for PacketIdx = 1:length(Channels{FromIdx, ToIdx})
                    Packet = Channels{FromIdx, ToIdx}{PacketIdx};
                    OuterDiameter = (CurrentTime - Packet.TransmitStartTime) * SpeedOfSound;
                    InnerDiameter =  (CurrentTime - Packet.TransmitStartTime - Packet.TransmitDuration) * SpeedOfSound;
                    if InnerDiameter < 0
                        % Still being transmitted
                        InnerDiameter = 0;
                    end

                    % Plot Circles
                    if OuterDiameter > 0
                        color = '#F7DC6F';
                        t = linspace(0, 2*pi, 32);
                        OuterCircle = polyshape(TxX + OuterDiameter*cos(t), TxY + OuterDiameter*sin(t), 'Simplify',false);
                        Wave = OuterCircle;
                        if InnerDiameter > 0
                            InnerCircle = polyshape(TxX + InnerDiameter*cos(t), TxY + InnerDiameter*sin(t), "Simplify",false);
                            Wave = subtract(OuterCircle, InnerCircle);
                        end
                        plot(Wave, "FaceColor",color, "EdgeColor",color);
                    end
                  

                end
            end
        end
    end
    hold off;

    %view(0, 90);
    view(2);

    axis square;
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

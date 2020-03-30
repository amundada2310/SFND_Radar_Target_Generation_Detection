clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant

%given values:
%fr = 77*10^9;
Rmax = 200;
Rr = 1;
Vmax = 100;
c = 3*10^8;

% defined values:
Rtarget = 100;
Vtarget = 50;

 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

Bsweep = c/(2*Rr); %Bandwidth for each chirp at given resolution
Tchirp  = 5.5*2*Rmax/c; %Chirp time for each chirp
slope = Bsweep/Tchirp; %Slope of the chirp signal

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = Rtarget+(Vtarget*t(i));% distance = speed * time
    td(i) = 2*r_t(i)/c; %twice the distnace to and fro
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*((fc*t(i))+((slope*t(i)*t(i))/2)));
    Rx(i) = cos(2*pi*((fc*(t(i)-td(i)))+((slope*(t(i)-td(i))*(t(i)-td(i)))/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Mix_fft = fft(Mix,Nr);
 % *%TODO* :
% Take the absolute value of FFT output and %normalize
Mix_fft = abs(Mix_fft/Nr);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
Mix_fft = Mix_fft(1:(Nr/2)+1);

%plotting the range
figure ('Name','Range from First FFT')
%subplot(2,1,1)

 % *%TODO* :
 % plot FFT output 

plot(Mix_fft); 
axis ([0 200 0 1]);
xlabel('Range'); 
ylabel('Amplitude'); 



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
ylabel('Range'); 
xlabel('Velocity'); 
zlabel('Amplitude');
title(' 2D FFT Output');

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;%20%14%12%12%10
Td = 8;%20%14%12%10%8

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5;%8%5
Gd = 5;%6%5
% *%TODO* :
% offset the threshold by SNR value in dB
Offset = 6;%6
% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);

signal_cfar = zeros(Nr/2,Nd);
% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
 
trainingcells = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1) - (2*Gr+1)*(2*Gd+1);%total training cells

for i = 1:(Nr/2 - (2*Gr+2*Tr))%outerloop of range
    for j = 1:(Nd - (2*Gd+2*Td))%inner loops of doppler
        
        %sig1-sig2 is done to only calculate training cells noise signal
        %(Tr and Td)- note we just considering the leading training cells
        %For every iteration sum the signal level within all the training cells.
         %To sum convert the value from logarithmic to linear using db2pow function.
        sig1 = sum(db2pow(RDM(i:i+2*Tr+2*Gr, j:j+2*Td+2*Gd)),'all');
        sig2 = sum(db2pow(RDM(i+Tr:i+Tr+2*Gr, j+Td:j+Td+2*Gd)),'all');    
        noise_level = sig1 - sig2;
        
        %Average the summed values for all of the training cells used for this -
        %convert back to the logarithmic using pow2db funtion.
        threshold = pow2db(noise_level/trainingcells);     
        
        %Further add the offset to it to determine the threshold.
        threshold = threshold + Offset;
        
        %Measuring the signal within the CUT
        signal = RDM(i+Tr+Gr, j+Td+Gd);
        
        % Filter the signal above the threshold
        % Next, compare the signal under CUT with this threshold. 
        %If the CUT level > threshold assign it a value of 1,
        % *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
        if (db2pow(signal) > db2pow(threshold))
            signal = 1;
        else 
            signal = 0;
        end
        
         if(signal ~= 0 & signal ~= 1)
             signal = 0;
         end
         
        signal_cfar(i+Tr+Gr,j+Td+Gd) = signal; 
        
    end
end

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;
xlabel('Doppler'); 
ylabel('Range');
zlabel('Amplitude');
title('CA-CFAR 2D Output');

figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;
xlabel('Doppler'); 
ylabel('Range');
zlabel('Amplitude');
title('CA-CFAR Range and Amplitude Output View');
view(90,0);
 
figure,surf(doppler_axis,range_axis,signal_cfar);
colorbar;
xlabel('Doppler'); 
ylabel('Range');
zlabel('Amplitude');
title('CA-CFAR Velocity and Amplitude Output View');
view(0,0); 
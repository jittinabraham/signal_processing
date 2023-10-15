clear all 
close all
A=1;
 fx=1000;
 ft=25.25*fx;
 At=10;
 F=500*fx;
 SNRdB=20
 n=0:1000;
 for N=1:length(n)
    Xn(N)=A*cos(2*pi*fx*n(N)/F);
    Xc(N)=At*cos(2*pi*ft*n(N)/F);
 end
Xt=transmitter_wave(Xn,Xc);
Xt_AWGN=noise_generator(SNRdB,Xt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Receiver SIde%%%%%%%%%%%%%%%%%%%%%%%

 
b=bandpass_filter(Xt_AWGN);

%%rectifier
rect_out = rectifier(b)

%%lowpass%%%%
td=1/(2*pi*(fx+ft));
output_signal=lowpass_filter(rect_out,F,td);

subplot(2, 3, 1);
plot(n, Xn);
title('input signal');

% Plot the second graph (1 row, 5 columns, second plot)
subplot(2, 3, 2);
plot(n, Xt);
title('transmitted signal');

% Plot the third graph (1 row, 5 columns, third plot)
subplot(2, 3, 3);
plot(n, Xt_AWGN);
title('Signal with AWGN noise');

% Plot the fourth graph (1 row, 5 columns, fourth plot)
subplot(2, 3, 4);
plot(n, b);
title('bandpass filterd');

% Plot the fifth graph (1 row, 5 columns, fifth plot)
subplot(2, 3, 5);
plot(n, rect_out);
title('Rectified signal');

subplot(2,3,6);
plot(output_signal);
title('output signal')







function rect_out = rectifier(b) 
        for n=1:length(b)
             if b(n)>0
             rect_out(n)=b(n);
             end
             if b(n)<0
             rect_out(n)=abs(b(n));
             end
        end
end






function f = bandpass_filter(e)
    % Initialize the state variables
    f = zeros(size(e)); % Initialize the output signal
    f_1 = 0; % Initialize f(n-1)
    f_2 = 0; % Initialize f(n-2)
    e_2 = 0; % Initialize e(n-2)

    % Filter coefficients
    a = [0.0245216, -0.0245216];
    b = [1.88873, -0.950957];

    % Apply the bandpass filter
    for n = 1:length(e)
        % Update f(n-2) and e(n-2)
        f_2 = f_1;
        e_2 = e(n);
        
        % Calculate f(n) using the filter equation
        f(n) = a(1) * e(n) - a(2) * e_2 + b(1) * f_1 - b(2) * f_2;
        
        % Update f(n-1) for the next iteration
        
    end
end

function output_signal=lowpass_filter(rect_out,F,td)
    output_signal=zeros(1,length(rect_out));
    output_signal=(1/(1+td*F))*rect_out(1);
    for i=2:1:length(rect_out)-1
      output_signal(i)=(1/(1+td*F))*rect_out(i)+ (1/(1+1/(td*F)))*output_signal(i-1);
    end
end

function Xt_AWGN = noise_generator(SNRdB,Xt)
signal_power = var(Xt)
% Calculate the power of the noise
noise_power = signal_power / (10^(SNRdB/10));

% Generate AWGN
noise = sqrt(noise_power) * randn(size(Xt));

% Add AWGN to the function output
Xt_AWGN = Xt + noise;

end

function Xt = transmitter_wave(Xn,Xc)
    Xt=(1+Xn).*Xc %modulated_signal
end


% Plot the first graph (1 row, 5 columns, first plot)



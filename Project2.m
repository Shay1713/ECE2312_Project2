
%% KEY

% creating audio recording parameters
Fs = 44100; 
nBits = 16; 
nChannels = 2; 
ID = -1;       % default audio input device 

load handel.mat % plotting parameters
window = hamming(512);
N_overlap = 256;
N_fft = 1024;

% plotting phrase 1 spectrogram
% [y1, Fs] = audioread('phrase1.wav');
% [S, F, T, P] = spectrogram(y1(:,1), window, N_overlap, N_fft, Fs, 'yaxis');
% figure;
% surf(T, F, 10*log10(P), 'edgecolor', 'none');
% axis tight;
% view(0,90);
% colormap(jet);
% set(gca, 'clim', [-80,-20]);
% ylim([0 8000]); xlabel('Time (seconds)'); % ylabel('Frequency (Hz)'); title(' Spectogram Phrase 1');



%% QUESTION 1

% creating audio recording parameters
Fs = 44100; % sample frequency
nBits = 16; 
nChannels = 2; 
ID = -1;       % default audio input device 
load handel.mat % plotting parameters
window = hamming(512);
N_overlap = 256;
N_fft = 1024;

% CITATION: https://www.mathworks.com/matlabcentral/answers/10575-how-to-check-length-of-a-wav-file-through-gui
[y,Fs] = audioread('phrase1.wav'); % reading the audio file into a variable
time = length(y)./Fs; % time in seconds
disp(time); % 6
disp(Fs); % 44100

% CITATION: https://www.mathworks.com/matlabcentral/answers/157727-generating-a-6-second-sine-wave 
w0=5000; % signal parameters
ts=1/Fs;
t=0:ts:time;
x=sin(2*pi*w0*t);
filename1 = 'team[8]-sinetone.wav'; % creating a .wav file for the sine tone
audiowrite(filename1,x,Fs); % writing and reading the sine tone into the .wav file
[x,Fs] = audioread(filename1);

[S, F, T, P] = spectrogram(x(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the sine tone in a spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-sinetone.wav'); % labels

% sound(x, Fs); % 6 second long beeeeeeeeeeppppppppp



%% QUESTION 2


% CITATION: https://www.mathworks.com/help/signal/ref/chirp.html 
%           https://www.mathworks.com/matlabcentral/answers/390583-how-can-i-generate-audio-chirp-signal

y = chirp(t,0,6,8000); % making the chrip with required parameters
%sound(y, Fs); % 6 second long chrip WARNING LOUD

filename2 = 'team[8]-chirp.wav'; % creating a .wav file to store the chirp
audiowrite(filename2,y,Fs); % writing and reading the file
[y,Fs] = audioread(filename2);

[S, F, T, P] = spectrogram(y(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-chirp.wav '); % labels


%% QUESTION 3

t1=0:ts:0.35; % first sound in sequence
w01=2200;
x1=sin(2*pi*w01*t1);
%sound(x1,Fs);
filenameBeep1 = 'beep1.wav';
audiowrite(filenameBeep1,x1,Fs);
[x1,Fs] = audioread(filenameBeep1);

t2=0:ts:0.7; % second sound in sequence
w02=2300;
x2=sin(2*pi*w02*t2);
%sound(x2,Fs);
filenameBeep2 = 'beep2.wav';
audiowrite(filenameBeep2,x2,Fs);
[x2,Fs] = audioread(filenameBeep2);

t3=0:ts:1.1; % third sound in sequence
w03=2000;
x3=sin(2*pi*w03*t3);
%sound(x3,Fs);
filenameBeep3 = 'beep3.wav';
audiowrite(filenameBeep3,x3,Fs);
[x3,Fs] = audioread(filenameBeep3);

t4=0:ts:0.7; % fourth sound in sequence
w04=1500;
x4=sin(2*pi*w04*t4);
%sound(x4,Fs);
filenameBeep4 = 'beep4.wav';
audiowrite(filenameBeep4,x4,Fs);
[x4,Fs] = audioread(filenameBeep4);

t5=0:ts:2.1; % fifth sound in sequence
w05=1800;
x5=sin(2*pi*w05*t5);
%sound(x5,Fs);
filenameBeep5 = 'beep5.wav';
audiowrite(filenameBeep5,x5,Fs);
[x5,Fs] = audioread(filenameBeep5);

% CITATION: https://www.mathworks.com/matlabcentral/answers/476508-concatenate-two-audio-files
combinedWav = [x1(1:Fs:numel(x1)); x2;x3;x4;x5]; % conocating all the sounds together
filenameBeepTotal = 'team[8]-cetk.wav'; % creating a .wav file to store the sound
audiowrite(filenameBeepTotal,combinedWav, Fs); % writing and reading the file
[combinedWav, Fs] = audioread(filenameBeepTotal);

[S, F, T, P] = spectrogram(combinedWav(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-cetk.wav'); % labels




%% QUESTION 4

[xP1, Fs] = audioread('phrase1.wav'); % calling the .wav file from project 1

timeP1 = length(xP1)./Fs; % 6 % double checking parameters did not change
disp(timeP1);
timex = length(x)./Fs; % 6.000
disp(timex);
disp(Fs); % 44100

szX2 = [x; x]; % making the sizes of the sine and phrase files  the same
% disp(size(szX2));  %   529200    1
singleX = reshape(szX2,[],2);
% disp(size(singleX));  %   264600     2

speechChirp = [xP1; singleX]; % conocating the phrase and sine signals
% disp(size(speechChirp)); 

% sound(speechChirp, Fs);

filenameSC = 'team[8]-speechchirp.wav'; % creating a .wav file to store sound
audiowrite(filenameSC, speechChirp, Fs); % wriiting and readnig the file
[speechChirp, Fs] = audioread(filenameSC);

[S, F, T, P] = spectrogram(speechChirp(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-speechchirp.wav'); % labels



%_______________________________________________________________

[xP12, Fs2] = audioread('phrase1.wav'); % calling the .wav file from project 1

timeP12 = length(xP12)./Fs; % 6 % double checking parameters did not change
disp(timeP12);
timex2 = length(x)./Fs; % 6.000
disp(timex2);
disp(Fs); % 44100

szX22 = [x; x]; % making the sizes of the sine and phrase files  the same
%disp(size(szX22)); 
singleX2 = reshape(szX22,[],2);
%disp(size(singleX2)); 
singleX2 = singleX2(1:end-1,:);
speechChirp2 = [xP12, singleX2]; % conocating the phrase and sine signals

%sound(speechChirp2, Fs);

filenameSC2 = 'team[8]-speechchirp2.wav'; % creating a .wav file to store sound
audiowrite(filenameSC2, speechChirp2, Fs); % wriiting and readnig the file
[speechChirp2, Fs] = audioread(filenameSC2);

[S, F, T, P] = spectrogram(speechChirp2(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-speechchirp2.wav'); % labels




%% QUESTION 5

% CITATION: https://www.mathworks.com/help/signal/ref/firls.html
% https://www.mathworks.com/help/signal/ref/lowpass.html


sigFilter = lowpass(speechChirp,4000,Fs); % lowpass filter of speech sine with cutoff freq of 4000Hz
% lowpass(speechChirp,4000,Fs); % plotting filtered signal
%sound(sigFilter, Fs); % playing sound

filenameFilter = 'team[8]-filteredspeechsine.wav'; % creating a .wav file to store sound
audiowrite(filenameFilter, sigFilter, Fs); % wriiting and reading the file
[sigFilter, Fs] = audioread(filenameFilter);

[S, F, T, P] = spectrogram(sigFilter(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-filteredspeechsine.wav'); % labels




% ________________________________________________________________

sigFilter2 = lowpass(speechChirp2,4000,Fs); % lowpass filter of speech sine with cutoff freq of 4000Hz
%lowpass(speechChirp2,4000,Fs); % plotting filtered signal
%sound(sigFilter, Fs); % playing sound

filenameFilter2 = 'team[8]-filteredspeechsine2.wav'; % creating a .wav file to store sound
audiowrite(filenameFilter2, sigFilter2, Fs); % wriiting and reading the file
[sigFilter2, Fs] = audioread(filenameFilter2);

[S, F, T, P] = spectrogram(sigFilter2(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-filteredspeechsine2.wav'); % labels




%% QUESTION 6


[phraseWav,Fs] = audioread('phrase1.wav'); % calling the phrase file again
zeroWav = zeros(size(phraseWav)); % making an array of zeros to make up for the file being shorter
conWav = [phraseWav, zeroWav]; % convolving with the zero array

phraseWav0 = reshape(conWav,[],2); % making the arrays the same dimension
sigFilter0 = sigFilter(1:end-1,:);
% disp(size(phraseWav0)); 
% disp(size(sigFilter0)); 
sigStereo = [phraseWav0 , sigFilter0]; % convolving to make a stero signal

% sound(sigStereo, Fs);

filenameStereo = 'team[8]-stereospeechsine.wav'; % creating a .wav file to store sound
audiowrite(filenameStereo, sigStereo, Fs); % wriiting and readnig the file
[sigStereo, Fs] = audioread(filenameStereo);

[S, F, T, P] = spectrogram(sigStereo(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-stereospeechsine.wav'); % labels



% _________________________________________________

[phraseWav,Fs] = audioread('phrase1.wav'); % calling the phrase file again

sigStereo2 = [phraseWav , sigFilter2]; % convolving to make a stero signal

% sound(sigStereo2, Fs);

filenameStereo2 = 'team[8]-stereospeechsine2.wav'; % creating a .wav file to store sound
audiowrite(filenameStereo2, sigStereo2, Fs); % wriiting and readnig the file
[sigStereo2, Fs] = audioread(filenameStereo2);

[S, F, T, P] = spectrogram(sigStereo2(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
figure;
surf(T, F, 10*log10(P), 'edgecolor', 'none');
axis tight;
view(0,90);
colormap(jet);
set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-stereospeechsine2.wav'); % labels



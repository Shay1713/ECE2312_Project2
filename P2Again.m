
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
% disp(time); % 6
% disp(Fs); % 44100

% CITATION: https://www.mathworks.com/matlabcentral/answers/157727-generating-a-6-second-sine-wave 
w0=5000; % signal parameters
ts=1/Fs;
t=0:ts:time;
% x = sin(2*pi*w0*t);

filename1 = writeReadFile('team[8]-sinetone.wav', sin(2*pi*w0*t), 44100);
% sound(filename1, Fs); % 6 second long beeeeeeeeeeppppppppp

% [S, F, T, P] = spectrogram(filename1(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the sine tone in a spectrogram
% figure;
% surf(T, F, 10*log10(P), 'edgecolor', 'none');
% axis tight;
% view(0,90);
% colormap(jet);
% set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
% xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-sinetone.wav'); % labels




%% QUESTION 2


% CITATION: https://www.mathworks.com/help/signal/ref/chirp.html 
%           https://www.mathworks.com/matlabcentral/answers/390583-how-can-i-generate-audio-chirp-signal


filename2 = writeReadFile('team[8]-sinetone.wav', chirp(t,0,6,8000), 44100);
%sound(filename2, Fs); % 6 second long chrip WARNING LOUD

% [S, F, T, P] = spectrogram(filename2(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
% figure;
% surf(T, F, 10*log10(P), 'edgecolor', 'none');
% axis tight;
% view(0,90);
% colormap(jet);
% set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
% xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-chirp.wav '); % labels


%% QUESTION 3

% % CITATION: https://www.mathworks.com/matlabcentral/answers/476508-concatenate-two-audio-files

sound1 = sounds(0:ts:0.35, 2200, 'beep1.wav', 44100); % individual sounds
sound2 = sounds(0:ts:0.7, 2300, 'beep2.wav', 44100);
sound3 = sounds(0:ts:1.1, 2000, 'beep3.wav', 44100);
sound4 = sounds(0:ts:0.7, 1500, 'beep4.wav', 44100);
sound5 = sounds(0:ts:2.1, 1800, 'beep5.wav', 44100);
combinedWav = [sound1(1:Fs:numel(sound1)); sound2;sound3;sound4;sound5]; % combining the sounds

filename3 = writeReadFile('team[8]-cetk.wav', combinedWav, 44100);
%sound(filename3, Fs); % cool tones

% [S, F, T, P] = spectrogram(filename3(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
% figure;
% surf(T, F, 10*log10(P), 'edgecolor', 'none');
% axis tight;
% view(0,90);
% colormap(jet);
% set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
% xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-cetk.wav'); % labels




%% QUESTION 4

[xP1, Fs] = audioread('phrase1.wav'); % calling the .wav file from project 1

timeP1 = length(xP1)./Fs; % 6 % double checking parameters did not change
disp(timeP1);
timex = length(filename1)./Fs; % 6.000
disp(timex);
disp(Fs); % 44100

szX2 = [filename1; filename1]; % making the sizes of the sine and phrase files  the same
% disp(size(szX2));  %   529200    1
singleX = reshape(szX2,[],2);
% disp(size(singleX));  %   264600     2

speechChirp = [xP1; singleX]; % conocating the phrase and sine signals
% disp(size(speechChirp)); 

filename4 = writeReadFile('team[8]-speechchirp.wav', speechChirp, 44100);
% sound(speechChirp, Fs);

% [S, F, T, P] = spectrogram(filename4(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
% figure;
% surf(T, F, 10*log10(P), 'edgecolor', 'none');
% axis tight;
% view(0,90);
% colormap(jet);
% set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
% xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-speechchirp.wav'); % labels


%_______________________________________________________________

[xP12, Fs2] = audioread('phrase1.wav'); % calling the .wav file from project 1

timeP12 = length(xP12)./Fs; % 6 % double checking parameters did not change
disp(timeP12);
timex2 = length(filename1)./Fs; % 6.000
% disp(timex2); % disp(Fs); % 44100

szX22 = [filename1; filename1]; % making the sizes of the sine and phrase files  the same
singleX2 = reshape(szX22,[],2);
%disp(size(szX22));  %disp(size(singleX2)); 
singleX2 = singleX2(1:end-1,:);
speechChirp2 = [xP12, singleX2]; % conocating the phrase and sine signals

filename5 = writeReadFile('team[8]-speechchirp2.wav', speechChirp2, 44100);
%sound(speechChirp2, Fs);

% [S, F, T, P] = spectrogram(filename5(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
% figure;
% surf(T, F, 10*log10(P), 'edgecolor', 'none');
% axis tight;
% view(0,90);
% colormap(jet);
% set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
% xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-speechchirp2.wav'); % labels





%% QUESTION 5

% CITATION: https://www.mathworks.com/help/signal/ref/firls.html
% https://www.mathworks.com/help/signal/ref/lowpass.html

sigFilter = lowpass(speechChirp,4000,Fs); % lowpass filter of speech sine with cutoff freq of 4000Hz
% lowpass(speechChirp,4000,Fs); % plotting filtered signal

filename6 = writeReadFile('team[8]-filteredspeechsine.wav', sigFilter, 44100);
%sound(sigFilter, Fs); % playing sound

%[S, F, T, P] = spectrogram(filename6(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
%figure;
%surf(T, F, 10*log10(P), 'edgecolor', 'none');
%axis tight;
%view(0,90);
%colormap(jet);
%set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
%xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-filteredspeechsine.wav'); % labels




% ________________________________________________________________

sigFilter2 = lowpass(speechChirp2,4000,Fs); % lowpass filter of speech sine with cutoff freq of 4000Hz
%lowpass(speechChirp2,4000,Fs); % plotting filtered signal

filename7 = writeReadFile('team[8]-filteredspeechsine.wav', sigFilter2, 44100);
%sound(sigFilter, Fs); % playing sound

%[S, F, T, P] = spectrogram(filename7(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
%figure;
%surf(T, F, 10*log10(P), 'edgecolor', 'none');
%axis tight;
%view(0,90);
%colormap(jet);
%set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
%xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-filteredspeechsine2.wav'); % labels


%%% firls below


% freqFs = []; % creating a Fs frequency array
% for i = 1:length(speechChirp2)
%     freqFs = [freqFs, Fs];
%     
% end
% 
% ampSC = []; % creating an amplitude array
% for i = 1:length(speechChirp2)
%     ampSC = [ampSC, speechChirp2(i)];
%     
% end
% 
% %disp(freqFs);
% %disp(ampSC);
% 
% %lpFilter = firls(4000,freqFs, ampSC);
% 
% %fvtool(lpFilter,'impulse');
% 
% 
% 
%% QUESTION 6


[phraseWav,Fs] = audioread('phrase1.wav'); % calling the phrase file again
zeroWav = zeros(size(phraseWav)); % making an array of zeros to make up for the file being shorter
conWav = [phraseWav, zeroWav]; % convolving with the zero array

phraseWav0 = reshape(conWav,[],2); % making the arrays the same dimension
sigFilter0 = sigFilter(1:end-1,:);
sigStereo = [phraseWav0 , sigFilter0]; % convolving to make a stero signal
% disp(size(phraseWav0)); % disp(size(sigFilter0)); 

filename8 = writeReadFile('team[8]-stereospeechsine.wav', sigStereo, 44100);
% sound(sigStereo, Fs);


%[S, F, T, P] = spectrogram(filename8(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
%figure;
%surf(T, F, 10*log10(P), 'edgecolor', 'none');
%axis tight;
%view(0,90);
%colormap(jet);
%set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
%xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-stereospeechsine.wav'); % labels



% _________________________________________________

[phraseWav,Fs] = audioread('phrase1.wav'); % calling the phrase file again
sigStereo2 = [phraseWav , sigFilter2]; % convolving to make a stero signal

filename9 = writeReadFile('team[8]-stereospeechsine2.wav', sigStereo2, 44100);
% sound(sigStereo2, Fs);

%[S, F, T, P] = spectrogram(filename9(:,1), window, N_overlap, N_fft, Fs, 'yaxis'); % plotting the spectrogram
%figure;
%surf(T, F, 10*log10(P), 'edgecolor', 'none');
%axis tight;
%view(0,90);
%colormap(jet);
%set(gca, 'clim', [-80,-20]); ylim([0 8000]); 
%xlabel('Time (seconds)'); ylabel('Frequency (Hz)'); title('Spectogram team[8]-stereospeechsine2.wav'); % labels




%% FUNCTION STORAGE

function [x1, Fs] = sounds(t1, w01, fileName, Fs) % QUESTION 3 SOUND FUNCTION
    x1=sin(2*pi*w01*t1);
    audiowrite(fileName,x1,Fs);
    [x1,Fs] = audioread(fileName);
end

function [x, Fs] = writeReadFile(fileNameDo, x, Fs) % READ AND WRITE AUDIO FILES 
    filename1 = fileNameDo; % creating a .wav file for the sine tone
    audiowrite(filename1,x,Fs); % writing and reading the sine tone into the .wav file
    [x,Fs] = audioread(filename1);
end





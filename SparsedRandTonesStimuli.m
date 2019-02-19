clear;clc;

%Variables:
Duty_cycle = 3/4;
JittVar_octaves = 1/12;
JittersTime_ms = 250;
JittersTime_s = JittersTime_ms/1000;
PlayingTime_s = JittersTime_s*Duty_cycle;
PureToneTotal_s = 4;
RepInPureTone = PureToneTotal_s/JittersTime_s;
nFreqs = 8;
SamplingRate_Hz = 44100;
ImageAcqTime_s = 1;
SilencePureTone_s = 3;
SoundPureTone = PlayingTime_s*SamplingRate_Hz;
PureToneTime_s = PureToneTotal_s + ImageAcqTime_s + SilencePureTone_s;
nRepsEachFreq = 10;
TotalTime_s = (((PureToneTime_s)*(nFreqs+1))*(nRepsEachFreq));
startfreq = 200;
endfreq = 8000;
EndPureTone = (PureToneTotal_s*SamplingRate_Hz);
EndAcqSilence = ((PureToneTotal_s+ImageAcqTime_s)*SamplingRate_Hz);

%Creates a list of octaves
OctaveList=(log2(startfreq):(log2(endfreq)-log2(startfreq))/(nFreqs-1):log2(endfreq))';

%Converts the list of octaves into frequencies
freqlist=pow2(OctaveList);

%Creates a random list of nFreqs * nRepsEachFreq
ListRandFreqs = RandOrderRowVector((nFreqs),(nRepsEachFreq));

%From the time 0 to (PlayingTime_s*SamplingRate_Hz), play the frequency with a jitter of
%random number between |JittVar_octaves|. Then silence from the
%(PlayingTime_s*SamplingRate_Hz) number and until the
%(JittersTime_s*SamplingRate_Hz) number is reached.

%Chooses the frequency to be played in each 4s tone block
RandFreqList = zeros(1,length(ListRandFreqs));
for ii = 0:(nFreqs)
    for hh = 1:length(ListRandFreqs)
        if ii == ListRandFreqs(hh) && ii == 0
        RandFreqList(hh) = 0;
        elseif ii == ListRandFreqs(hh) && ii ~= 0
        RandFreqList(hh) = freqlist(ii);    
        else
            continue
        end
    end
end


%%%Sets the characteristics of the whole audio%%%
%The blocks of 8s are randomly ordered in a 12 min audio.


for repSound = 1:length(RandFreqList) %Repeats the ToneBlock of 8s during the whole audio.

    %%%Sets the characteristics of each complete tone block (8s)%%%

    ToneBlockLength = (PureToneTotal_s+ImageAcqTime_s+SilencePureTone_s)*SamplingRate_Hz;

    ToneBlock_8s = zeros(((ToneBlockLength)),6);
    %Defines the matrix so that it doesn't change size each loop.

    ToneBlock_8s (:, 1) = (0:(1/SamplingRate_Hz):(PureToneTime_s-(1/SamplingRate_Hz)))';
    %Defines the 1st column of the matrix and this is equal to the 'sampling 
    %rate counter / time' that is going to be used to produce de wave.


    %%%Sets the characteristics of each 4s pure tone block%%%

    PureToneBlock_4s = zeros(EndPureTone,5);
    %Defines the matrix so that it doesn't change size each loop.

    PureToneBlock_4s (:, 1) = (0:(1/SamplingRate_Hz):PureToneTotal_s-(1/SamplingRate_Hz))';
    %Defines the 1st column of the matrix and this is equal to the 'sampling 
    %rate counter / time' that is going to be used to produce de wave.


    for kk = 0:(RepInPureTone-1) %Repeats the 250 tone block for the given number of repetitions.
        ToneBlock_250ms = zeros((JittersTime_s*SamplingRate_Hz),1);
        ToneBlock_250ms (:,2) = 1;

        %Create a cos gated ramp so the amplitude increases slowly until it is
        %one in each Pure Tone Block
        ampramp_time = 20; % ramp time, in ms
        ramp = 0:1/SamplingRate_Hz:(ampramp_time/1000);
        gate = cos((-pi:pi/length(ramp):0));
        gate = (gate+1)/2;

        upramp=gate';
        downramp=(fliplr(gate))'; 

        for mm = 1:(JittersTime_s*SamplingRate_Hz) %Sets the characteristics of each 250 tone block
           if mm <= SoundPureTone
               if mm <= length(upramp)
                   ToneBlock_250ms (mm,1) = 1;
                   ToneBlock_250ms (mm,2) = upramp (mm); %apply upramp
               elseif mm > length(upramp) && mm < ((SoundPureTone)-(length(downramp)))
                   ToneBlock_250ms (mm,1) = 1;
                   ToneBlock_250ms (mm,2) = 1;
               else %mm > length(upramp) && mm >= ((JittersTime_s*Duty_cycle*SamplingRate_Hz)-(length(downramp)+1))
                   ToneBlock_250ms (mm,1) = 1;
                   ToneBlock_250ms (mm,2) = downramp (ceil(length(downramp)-(SoundPureTone-(mm)))); %apply downramp
               end
           else
               ToneBlock_250ms (mm, 1) = 0;
               ToneBlock_250ms (mm, 2) = 0;
           end 
        end

        ptb_4s_timing = ((kk*JittersTime_s*SamplingRate_Hz)+1):((kk+1)*JittersTime_s*SamplingRate_Hz);

        %Defines the 2nd column of the matrix that specifies if the sound
        %is played (1) or if it's silence (0)
        PureToneBlock_4s(ptb_4s_timing,2)=ToneBlock_250ms(:,1);

        %Defines the 3rd column of the matrix that specifies the gating of each
        %block of sound.
        PureToneBlock_4s(ptb_4s_timing,3)=ToneBlock_250ms(:,2);

        %Defines the 4th column of the matrix that gives the played 
        %frequency and converts the zeros to silence.
        if RandFreqList(repSound) == 0
            PureToneBlock_4s(ptb_4s_timing,4) = 0;
            PureToneBlock_4s(:,2) = 0;
        else
            PureToneBlock_4s(ptb_4s_timing,4)= RandFreqList(repSound)+ pow2(RealRand(JittVar_octaves,1)); 
        end

        %Defines the 5th column of the matrix that makes the wave with the
        %formula 'sine wave = A(sin(2*pi*freq*time))' where time is the second 
        %column, and freq is the 3th column;
        PureToneBlock_4s(ptb_4s_timing,5)= PureToneBlock_4s(ptb_4s_timing,2).*PureToneBlock_4s(ptb_4s_timing,3).*(sin(2*pi...
            .*PureToneBlock_4s(ptb_4s_timing,1).*(PureToneBlock_4s(ptb_4s_timing,4))));

    end
    
    %Sets the characteristics of the 8s blocks the tone (4s) after the silence (4s)
    %Sets the silence before the block of tones
    ToneBlock_8s(1:EndPureTone,2:5) = 0;
        
    %Sets the tones after the silence
    ToneBlock_8s((EndPureTone+1):ToneBlockLength,2:5) = PureToneBlock_4s(:,2:5);
    
        
    STS_Blocks{repSound} = ToneBlock_8s; %#ok<SAGROW>
    MAT_file{repSound} = ToneBlock_8s(:,5); %#ok<SAGROW>
end

Sound_STS = cell2mat(STS_Blocks');

Sound_STS (:,1) = (0:(1/SamplingRate_Hz):(TotalTime_s)-(1/SamplingRate_Hz))';


%Add the noise to the whole audio
    
%%%%Creates the column with the noise
%Converts from decibels to amplitude
dB_tone = 75;
dB_noise = dB_tone - 40;
Amp_tone = Amp_to_dB(dB_tone);
Amp_noise = Amp_to_dB(dB_noise);
Proportion = (Amp_noise / Amp_tone); %As the amplitude of the tone is 1, this calculates the proportion of the noise in respect of 1.
noise = -Proportion + (Proportion+Proportion)*rand(1,length(Sound_STS))'; %Multiplies per the calculated proportion (40 dB below the tone)
Reduction_amp = round((1-(max(noise))),3);

%ADD UPRAMP and DOWNRAMP to the noise
Gated_noise = noise;
Gated_noise(1:length(upramp)) = Gated_noise(1:length(upramp)).*upramp;
Gated_noise((end-length(downramp)+1):end) = Gated_noise((end-length(downramp)+1):end).*downramp;

%Reduces the amplitude of the tone column to let the noise be summed to it.
Sound_STS (:,6) = Sound_STS (:,5) .* Reduction_amp;

%Adds the column that is equal to the noise
Sound_STS (:,7) =  Gated_noise .* Reduction_amp;

%%%Creates a MAT file with the noise
noiseVector = -Proportion + (Proportion+Proportion)*rand(1,(SamplingRate_Hz*((12*60)+16))); %12*60+16 = number of seconds to play the noise
Gated_noiseVector = noiseVector';
Gated_noiseVector(1:length(upramp)) = Gated_noiseVector(1:length(upramp)).*upramp;
Gated_noiseVector((end-length(downramp)+1):end) = Gated_noiseVector((end-length(downramp)+1):end).*downramp;
save('Noise_run7.mat','Gated_noiseVector');

%Adds a column that is the sum of the tone and the noise
Sound_STS (:,8) = Sound_STS (:,6) + Sound_STS (:,7);

%%Creates the MAT file of the tones
TOTAL_MAT_file = cell2mat(MAT_file);
TOTAL_MAT_file_reduced = TOTAL_MAT_file .* Reduction_amp; %Needs to be reduced to be the same as the .wav
TOTAL_MAT_file_reduced = TOTAL_MAT_file_reduced';
save('SparsedRandTonesStimuli_run7.mat','TOTAL_MAT_file_reduced');

%Creates the files .txt to be printed
PureToneBlock_250ms_list = 2:JittersTime_s*SamplingRate_Hz:length(Sound_STS); %Begins in 2 because #1 = 0
ToneBlock_8s_list = ((PureToneTotal_s*SamplingRate_Hz)+1):(PureToneTime_s*SamplingRate_Hz):length(Sound_STS);
JitteredFrequencies_list = [Sound_STS(PureToneBlock_250ms_list,1) Sound_STS(PureToneBlock_250ms_list,4)];
Frequencies_list = [Sound_STS(ToneBlock_8s_list,1) Sound_STS(ToneBlock_8s_list,4)];

audiowrite(('TOTALSparsedTonesStimuli_tones_run7.wav'),(Sound_STS (:,6))',SamplingRate_Hz);
audiowrite(('TOTALSparsedTonesStimuli_noise_run7.wav'),(Sound_STS (:,7))',SamplingRate_Hz);
audiowrite(('TOTALSparsedTonesStimuli_complete_run7.wav'),(Sound_STS (:,8))',SamplingRate_Hz);
audiowrite(('TOTALSparsedTonesStimuli_MEG_TTL_run7.wav'),(Sound_STS (:,2))',SamplingRate_Hz);
dlmwrite('Frequencies_list_run7.txt',Frequencies_list,'delimiter','\t', 'precision', 9);
dlmwrite('JitteredFrequencies_list_run7.txt',JitteredFrequencies_list,'delimiter','\t', 'precision', 9);
fprintf('Everything done. I am ready for more.\n');

%Finally, SOUND_STS is a 9 column matrix in which the columns represent:
    %1. The Timing (1/sampling rate) with a duration of 24 minutes.
    %2. The TTL for MEG. Whenever there is sound in the file.
    %3. The gated TTL applying a gate in each 250 ms block.
    %4. The frequency that is played (with 1 semitone variations each 250
    %ms, and a change of frequency each 36 s).
    %5. The wave of the tone with an amplitude of 1.
    %6. The wave of the tone with an amplitude of .9 .
    %7. The signal of the noise.
    %8. The final audio with the noise added to it.
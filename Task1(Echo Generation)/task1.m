#variables

# x:input signal in time domain(original audio)
# X:input signal in frequency domain
# h:impulse response
# H:frequency response
# y:output audio in time domain(audio with echo)
# Y:output audio in frequency domain
# x2:audio after echo removal in time domain 


#read audio file
[x,fs]=audioread("audio1.wav");
original_audio_length = length(x)

#plot spectrum of signal
figure(1)
plot(x);
title("Input Audio spectrum in time domain")
ylabel("x[n]")
xlabel("n")

#calculate h[n]
h=zeros(fs,1);
h(cast(ceil(0.25*fs),"int32"))=0.7;
h(cast(ceil(0.5*fs),"int32"))=0.8;
h(cast(ceil(0.75*fs),"int32"))=0.9;
h(fs)=1;

#plot h[n]
figure(2)
plot(h);
title("Impulse response")
ylabel("h[n]")
xlabel("n")

#calculate echo function using convolution
y=conv(x,h);
figure(3)
plot(y)
title("Output Audio spectrum in time domain")
ylabel("y[n]")
xlabel("n")

#write output audio file after apply echo
audiowrite("EchoAudio.wav",y,fs);

#Echo removal part
L=max([length(y);length(h)]);
#make sure that length is not less than y by padding the difference with zeros
H=fft([h;zeros(L-length(h),1)]);  
Y=fft([y;zeros(L-length(y),1)]);  

#find x2(decoded audio) using inverse of forier transform with original length
x2 = real(ifft(Y ./ H))

#removes very small values due to error in precision of calculations which causes redundant memory and time
x2 = x2(abs(x2) > 10^-15);

#plotting decoded audio after echo removal 
figure(4)
plot(x2,"r");
title("Decoded Audio spectrum in time domain after echo removal")
ylabel("x2[n]")
xlabel("n")

#write decoded file
audiowrite("decodedAudio.wav",x2,fs)

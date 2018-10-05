%return logorithmically spaced PSD curves for the Peterson Models given a
%logorithmic period spacing of dlP and sampling rate fs
%lpd1 and lpd2 are the log Peterson parameterization periods for the low and high
%noise models, respectively
%Reference is Peterson, 1993
function [LNMA,HNMA,lpd1,lpd2]=peterson_acc(dlP,fs)
%Peterson low-noise acceleration parameterization values  (Peterson, 1993)
%periods
P1=[1.01/(fs/2) .1 .17 .40 .80 1.24 2.40 4.30 5.0 6.0 10.0 12.0 15.6 21.9 31.6 45 70 101 154 328 600];
A1=[-162.36 -162.36 -166.7 -170 -166.4 -168.6 -159.98 -141.10 -71.36 -97.26 -132.18 -205.27 -37.65 -114.37 ...
    -160.58 -187.50 -216.47 -185.00 -168.34 -217.43];
B1=[5.64 5.64 0.00 -8.30 28.90 52.48 29.81 0.00 -99.77 -66.49 -31.57 36.16 ...
    -104.33 -47.10 -16.28 0.00 15.70 0.00 -7.61 11.90];

%Peterson high-noise acceleration parameterization values
P2=[1.01/(fs/2) .1 .22 .32 .80 3.80 4.60 6.30 7.90 15.40 20 354.8];
A2=[-108.73 -108.73 -150.34 -122.31 -116.85 -108.48 -74.66 0.66 -93.37 73.54 -151.52];
B2=[-17.23 -17.23 -80.5 -23.87 32.51 18.08 -32.95 -127.18 -22.43 -162.98 10.01];

%use logorithmic frequency spacing to generate Peterson Curves
dlP1=dlP;
lP1=log10(P1);
k=1;
for i=1:length(P1)-1

    for lp1=lP1(i):dlP1:lP1(i+1)-dlP1

        lpd1(k)=lp1;

        %Acceleration PSD
        LNMA(k)=A1(i)+B1(i)*lp1;
        %Velocity PSD
        %LNMV(k)=A1(i)+B1(i)*lp1+20*(lp1-log10(2*pi));
        k=k+1;

    end
end

dlP2=dlP;
lP2=log10(P2);
k=1;
for i=1:length(P2)-1

    for lp2=lP2(i):dlP2:lP2(i+1)-dlP2

        lpd2(k)=lp2;

        %Acceleration PSD
        HNMA(k)=A2(i)+B2(i)*lp2;
        %Velocity PSD
        %HNMV(k)=A2(i)+B2(i)*lp2+20*(lp2-log10(2*pi));
        k=k+1;

    end
end
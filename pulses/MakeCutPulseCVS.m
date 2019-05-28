function [ energyRatio ] = MakeCutPulseCVS( pulseFWHM, cutFWHM, residual )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gauss = @(a, b, c, x) a*exp(-(x-b).^2/(2*(c/2.35)^2));
gerf = @(a, b, c, x) a+(1-a)*0.5*(1+erf((x-b)/sqrt((2*(c/2.35)^2))));

t = -pulseFWHM*4:0.05:pulseFWHM*4;

pulse = gauss(1, 0, pulseFWHM, t);
pulseArea = trapz(t, pulse);


cut = gerf(residual, 0, cutFWHM, -t);
pulse(t>0) = 1;%so that the cut persists
cutArea = trapz(t, (pulse.*cut));
%cut = 1;

pulse = pulse.*cut;

figure;
plot(t, pulse); hold on
firstVal = find(pulse > 1e-10, 1, 'first');
lastVal = find(pulse > 1e-10, 1, 'last');

t = t(firstVal:lastVal);
pulse = pulse(firstVal:lastVal);
plot(t, pulse); 



pulseMat = [t', (pulse)'];

fileName = sprintf('pulse%gps_cut%gps.csv', pulseFWHM, cutFWHM)
dlmwrite(fileName, pulseMat, 'delimiter', ',', 'precision', '%0.6f');


energyRatio = cutArea ./ pulseArea;

%plot(t,  pulse*cutArea./energyRatio)
end


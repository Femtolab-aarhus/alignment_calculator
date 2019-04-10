function [ energyRatio ] = MakeCutPulseWithKickCVS( pulseFWHM, cutFWHM, residual, kickHeight, kickPos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gauss = @(a, b, c, x) a*exp(-(x-b).^2/(2*(c/2.35)^2));
gerf = @(a, b, c, x) a+(1-a)*0.5*(1+erf((x-b)/sqrt((2*(c/2.35)^2))));

%if kickPos + 2  < pulseFWHM*2
    t = -pulseFWHM*2:0.05:pulseFWHM*2-0.0;
%else
 %   t = -pulseFWHM*2:0.05:kickPos + 5;
%end
    



pulse = gauss(1, 0, pulseFWHM, t);
pulseArea = trapz(t, pulse);


cut = gerf(residual, 0, cutFWHM, -t);
pulse(t>0) = 1;%so that the cut persists

%now add the kick
cutPulse = pulse.*cut + gauss(kickHeight, kickPos, 0.5, t);
cutArea = trapz(t, (cutPulse));


figure;
plot(t, cutPulse);

pulseMat = [t', (cutPulse)'];

fileName = sprintf('pulse%gps_cut%gps_residual%g_kick%gat%g.csv', pulseFWHM, cutFWHM, residual, kickHeight, kickPos)
dlmwrite(fileName, pulseMat, 'delimiter', ',', 'precision', '%0.5f');


energyRatio = cutArea ./ pulseArea;

plot(t,  cutPulse*cutArea./energyRatio)
end


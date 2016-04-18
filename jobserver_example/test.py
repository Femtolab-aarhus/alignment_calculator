#!/usr/bin/python3

import sqlite3
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy
import propagation
import config

I_max = 9.6e16;
K,M,KMsign = 0,0,1;
FWHM = 330e-15;

molecule = config.molecule("OCS","conf/molecules.conf");

t = 3*FWHM*numpy.array([-1,1]);
transfer_matrix120 = propagation.get_transfer_matrix(I_max,FWHM,t,120,K,M,KMsign,molecule);
transfer_matrix160 = propagation.get_transfer_matrix(I_max,FWHM,t,160,K,M,KMsign,molecule);

a = numpy.zeros_like(transfer_matrix160);
a[0:121,0:121] = transfer_matrix120;

transfer_matrix = transfer_matrix160-a;

print(numpy.where(transfer_matrix!=0));

plt.matshow(numpy.real(transfer_matrix));
plt.title('Real part');
plt.colorbar();
plt.matshow(numpy.imag(transfer_matrix));
plt.title('Imaginary part');
plt.colorbar();
plt.matshow(numpy.abs(transfer_matrix160)**2-numpy.abs(a)**2);
plt.title('Norm squared');
plt.colorbar();
plt.matshow((numpy.angle(transfer_matrix160)-numpy.angle(a))*180/numpy.pi);
plt.title('Angle');
plt.colorbar();


plt.show();

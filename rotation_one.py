# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

for i in range(1,31):
	nombres='AURIGA_ALEX/dat_auriga/py/stars_main_H%s_ids_127.mat' %(i)
	print nombres
	diccio= sio.loadmat('%s' %(nombres))  #cargar datos a python
	catgo= diccio['a'] #usar solo los datos numericos
	# 0)X  ; 1)Y  ; 2)Z  ; 3)vx, 4)vy, 5)vz, 6)mass, 7)Ids
	catgo[:,0:3] = 1000.*catgo[:,0:3]/0.6777 #pasar de mpc a kpc
	#viene normalizado por h=0.6777 
	catgo[:,6] = catgo[:,6]*(1.e10/0.6777) #masas solares
	catgo[:,7] = catgo[:,7]/1.e9 #giga años

	lb_time = 0 #segun el snapshot (ahora estamos a z=0)

	print catgo.shape

	r_loc = np.sqrt((posX**2) + (posY**2) + (posZ**2))
	
	#selecion de estrellas jovenes menores a 3gyr para obtener disco frio
	age_sel = np.where((catgo[:,7] <= 3 + lb_time)&(catgo[:,7] > 0.))[0]

	print len(age_sel)

	h_0 = 0.6777 
	omega_lambda= 0.693
	omega_m=0.307
	fac_esc=1 #factor de escala (cambiar segun snapchot)
	


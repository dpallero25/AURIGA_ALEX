# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

log = np.zeros([1,5])

for i in range(1,31):
    nombres='../dat_auriga/py/stars_main_H%s_ids_127.mat' %(i)
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

    #selecion de estrellas jovenes menores a 3gyr para obtener disco frio
    age_sel = np.where((catgo[:,7] <= 3 + lb_time)&(catgo[:,7] > 0.))[0]

    print len(age_sel)

    h_0 = 0.6777 
    omega_lambda= 0.693
    omega_m=0.307
    fac_esc=1 #factor de escala (cambiar segun snapchot)

    w = omega_m / (omega_lambda*(fac_esc**3))
    hubble_t = h_0 * np.sqrt(omega_lambda * np.sqrt(1+w))

	#distancia de objetos al centro
    r_loc = np.sqrt(np.square(catgo[:,0]) + np.square(catgo[:,1]) + np.square(catgo[:,2]))

	#seleccion a menos de 10 kpc
    inloc = np.where(r_loc <= 10.0)[0]


	#obtener velocidades en unidades fisicas, dentro del radio
    for j in range(3,6):
        catgo[:,j] = (hubble_t * catgo[:,j-3] + fac_esc * catgo[:,j] - np.mean(catgo[inloc,j]))

	#velocidad promedio	
    mean_v	= np.mean(catgo[inloc,3:6])
    print('mean velocity: %s' %(mean_v))

    #establecer limites para centrar muestra
    Zlim=6.
    Rmax=7.
    Rmin=0.
    Rlim=5
	#rotacion iterativa 4 veces 
    for k in range(4):
        print k                
        edad = catgo[age_sel,0:6]
        R = np.sqrt(np.square(edad[:,0]) + np.square(edad[:,1]))

        ind=np.where( (R <= Rmax)&(R >= Rmin))[0]

        atemp = edad[ind,:]

        #obtener la componente para momentum
        #producto cruz para obtener momentum (lx=p_y*v_z + p_z*v_y)
        l1=atemp[:,1]*atemp[:,5] - atemp[:,2]*atemp[:,4]
        l2=atemp[:,2]*atemp[:,3] - atemp[:,0]*atemp[:,5]
        l3=atemp[:,0]*atemp[:,4] - atemp[:,1]*atemp[:,3]

        lt1=np.sum(l1)
        lt2=np.sum(l2)
        lt3=np.sum(l3)
        lt4=np.sqrt(np.square(lt1)+np.square(lt2)+np.square(lt3))

        lt=[lt1,lt2,lt3,lt4]

        print 'original orientation'
        print 'lx/lt', lt[0]/lt[3]
        print 'ly/lt', lt[1]/lt[3]
        print 'lz/lt', lt[2]/lt[3]
        print '#################################################################'
        
        #Definir angulo de rotacion
        #situacion para cada cuadrante
        
        theta = np.rad2deg(np.arccos(lt[2]/lt[3]))
        
        if  lt[0] > 0. and lt[1] > 0 :
            phi = -np.rad2deg(np.arctan(lt[0]/lt[1]))
            print 'case 1'

        elif lt[0] > 0. and lt[1] < 0. :
            phi = -(180. + np.rad2deg(np.arctan(lt[0]/lt[1])))
            print 'case 2'

        elif ( lt[0] < 0.) and (lt[1] < 0.)  :
            phi = -(180. + np.rad2deg(np.arctan(lt[0]/lt[1])))
            print 'case 3' 	

        elif ( lt[0] < 0. and lt[1] > 0.  ):
            phi = -(360. + np.rad2deg(np.arctan(lt[0]/lt[1])))
            print 'case 4'
            
        
        # ROTACION EN 3 EJES #
        print 'begin rotation'
        print 'please wait a moment'
        print ' *_* '
        
        #primera rotacion#
        datrot = np.copy(catgo)
        datrot[:,0] = datrot[:,0]*np.cos(np.deg2rad(phi)) + datrot[:,1]*np.sin(np.deg2rad(phi))
        datrot[:,1] = -datrot[:,0]*np.sin(np.deg2rad(phi)) + datrot[:,1]*np.cos(np.deg2rad(phi))
        datrot[:,2] = datrot[:,2]
        datrot[:,3] = datrot[:,3]*np.cos(np.deg2rad(phi)) + datrot[:,4]*np.sin(np.deg2rad(phi))
        datrot[:,4] = -datrot[:,3]*np.sin(np.deg2rad(phi)) + datrot[:,4]*np.cos(np.deg2rad(phi))
        datrot[:,5] = datrot[:,5]
        
        print ' O.X '
        #segunda rotacion#
        datrot[:,0] = datrot[:,0] - 0
        datrot[:,1] = np.sin(np.deg2rad(theta))*datrot[:,2] + np.cos(np.deg2rad(theta))*datrot[:,1]
        datrot[:,2] = datrot[:,2]*np.cos(np.deg2rad(theta)) - datrot[:,1]*np.sin(np.deg2rad(theta))
        datrot[:,3] = datrot[:,3] - 0
        datrot[:,4] = datrot[:,5]*np.sin(np.deg2rad(theta)) + datrot[:,4]*np.cos(np.deg2rad(theta))
        datrot[:,5] = datrot[:,5]*np.cos(np.deg2rad(theta)) - datrot[:,4]*np.sin(np.deg2rad(theta))
        
        print ' X.0 '
        
        #tercera rotacion#
        phi = -phi

        datrot[:,0] = datrot[:,0]*np.cos(np.deg2rad(phi)) + datrot[:,1]*np.sin(np.deg2rad(phi))
        datrot[:,1] = -datrot[:,0]*np.sin(np.deg2rad(phi)) + datrot[:,1]*np.cos(np.deg2rad(phi))
        datrot[:,2] = datrot[:,2] - 0
        datrot[:,3] = datrot[:,3]*np.cos(np.deg2rad(phi)) + datrot[:,4]*np.sin(np.deg2rad(phi))
        datrot[:,4] = -datrot[:,3]*np.sin(np.deg2rad(phi)) + datrot[:,4]*np.cos(np.deg2rad(phi))
        datrot[:,5] = datrot[:,5] - 0
        
        print 'X.X'
        
        ## guardar coordenadas rotadas ##
        
        af=np.copy(datrot)
        
        if i == 0:
            Zlim = 4
            Rlim = 5
        
        else: Zlim = Zlim - 1
        
        Rlim = Rlim - 1
        
        #volver a extraer momentum con las mismas restricciones para rotar#
        l1 = af[ind,1]*af[ind,5] - af[ind,2]*af[ind,4]
        l2 = af[ind,2]*af[ind,3] - af[ind,0]*af[ind,5]
        l3 = af[ind,0]*af[ind,4] - af[ind,1]*af[ind,3]
        lt1 = np.sum(l1)
        lt2 = np.sum(l2)
        lt3 = np.sum(l3)
        lt4 = np.sqrt(np.square(lt1) + np.square(lt2) + np.square(lt3))
        lt = [lt1, lt2, lt3, lt4]
        print 'final'
        print 'lx/lt : ', lt[0]/lt[3]
        print 'ly/lt : ', lt[1]/lt[3]
        print 'lz/lt : ', lt[2]/lt[3]
    
    print '################End#of#Rotation#######################'
        
        
    #devolver los momentum completos pero ya rotados#
    l1 = af[:,1]*af[:,5] - af[:,2] * af[:,4]
    l2 = af[:,2]*af[:,3] - af[:,0] * af[:,5]
    l3 = af[:,0]*af[:,4] - af[:,1] * af[:,3]
    lt1 = np.sum(l1)
    lt2 = np.sum(l2)
    lt3 = np.sum(l3)
    lt4 = np.sqrt(np.square(lt1) + np.square(lt2) + np.square(lt3))
    lt = [lt1, lt2, lt3, lt4]
    
    print 'final'
    print 'lx/lt : ', lt[0]/lt[3]
    print 'ly/lt : ', lt[1]/lt[3]
    print 'lz/lt : ', lt[2]/lt[3]
    print '##############'
    #af[:,6] = mass
    #af[:,7] = edad
	#af[:,8] = pot
    #af[:,8] = ids
    #af[:,9] = feh
    
    dat = np.asarray([i, mean_v,lt[0]/lt[3],lt[1]/lt[3],lt[2]/lt[3]])
    log = np.vstack([log, dat])
    print np.shape(af)
    
    np.savetxt('../dat_au_alx/%s_rotated.dat' %(i), af)
    
log = np.delete(log,0,0)
np.savetxt('../dat_au_alx/log_rotacion' , log)
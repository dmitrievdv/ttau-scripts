import readmodel as rm
import numpy as np
import os

def cartesian(*X):
	if(len(X) == 1):
		return X[0]
	res = []
	X = list(X)
	# print(X)
	recursive = False
	for x in X[0]:
		Y = X[1:]
		# print(Y)
		Y = cartesian(*Y)
		for y in Y:
			# print('fory: ',y)
			# print(y.__class__)
			if(y.__class__ == tuple):
				res.append((x,)+y)
			else:
				res.append((x,y))
	return res

def Mdot_to_string(Mdot):
	Mdot10 = int(-np.log10(Mdot)*10)
	return str(Mdot10)
	# if Mdot10%10 == 0:
	# 	return str(Mdot10//10)
	# else:
	# 	return str(Mdot10)


model_dir = 'models/data/'
os.makedirs(model_dir, exist_ok = True)

logMdots = [-8.0, -8.5, -9.0]
logMdots = np.array(logMdots)
Mdots = 10.0**logMdots
# print(Mdots)
# quit()
Tmaxs = [7000]
rmis = [2.2, 4.2]
drm = [0.8, 1.8]
mags = cartesian(rmis, drm)
# mags = list(zip(rmis, drm))

models = cartesian(Mdots,Tmaxs,mags)

# print(models)

# quit()

Mstar = 0.8
Rstar = 2
Tstar = 4000
hotspot = False
Thot = 5000
dhot = 0.1

star_name = 'hart94'
model_names = []

for Mdot,Tmax,rmi,drm in models:
	rmo = rmi+drm
	# rmi = rmo-dr
	model_name = (star_name+'_'
		         +Mdot_to_string(Mdot)+'_'
				 +str(int(Tmax/100))+'_'
				 +str(int(rmi*10)))+'-'+str(int(rmo*10))
	if (hotspot): 
		model_name = model_name + ('_hot'
				 +str(int(Thot))+'l')
				 
	model_names.append(model_name)
	model_data = ('# Parameters for model '+model_name+'\n'+
				  str(Mstar)+' # Mstar\n'+
				  str(Tstar)+' # Tstar\n'+
				  str(Rstar)+' # Rstar\n'+
				  str(Mdot)+' # Mdot\n'+
				  str(Tmax)+' # Tmag\n'+
				  'dipole # field type\n'+
				  '# borders\n'+
				  str(rmi)+' # first\n'+
				  str(rmo)+' # second\n')
	if (hotspot):
		model_data = model_data + ('# hot spot:\n'+
				  str(dhot)+'# Dhot\n'+
				  str(Thot)+'# Thot')
	data_file = open(model_dir+model_name+'_data.dat', 'w')
	data_file.write(model_data)
	data_file.close()


print(model_names)
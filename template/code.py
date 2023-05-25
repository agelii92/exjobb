

import os
import random
import csv
import numpy as np
from functions import make_SASE
import sobol

def random_photon_energy():
    if photon_energy_variable:
        return str(int(random.choice(X_ray_Ebins)))
    else:
        return photon_energy

def crash_line(x):
    k = 1700/0.2e6
    m = -500
    if x < 0.4e6:
        return k*x+m
    else:
        return k*0.4e6+m
    
#pre-defined discrete values of the photon energy
X_ray_Ebins = np.logspace(-1,4,500)
idx = [n for n,i in enumerate(X_ray_Ebins) if (i >= min_photon_energy and i <= max_photon_energy)]
X_ray_Ebins = list(X_ray_Ebins[idx])    

#Update simulation file to the start/stop values
directory = f'runs_{start}-{end}'
os.system(f'mkdir {directory}')
os.system(f'sed -i \'s/<MODIFY_START>/' + str(start) + f'/g\' ./get_database.sh')
os.system(f'sed -i \'s/<MODIFY_END>/' + str(end) + f'/g\' ./get_database.sh')
os.system(f'cp davinci.sh ./{directory}/davinci_runs_{start}-{end}.sh')
os.system(f'sed -i \'s/<START>/' + str(start) + f'/g\' ./{directory}/davinci_runs_{start}-{end}.sh')
os.system(f'sed -i \'s/<END>/' + str(end) + f'/g\' ./{directory}/davinci_runs_{start}-{end}.sh')

num_variables = sum([photon_energy_variable, pulse_length_variable, fluence_variable, density_variable, water_percent_variable, hydrogen_variable, carbon_variable, nitrogen_variable, oxygen_variable])

run_no = list(range(start,end+1))

for run in run_no:
#generate parameters: the rescaling algorithm is to scale the [0,1] samples to correct magnitude
        
    if fluence_variable:
        fluence = random.uniform(min_fluence,max_fluence)

    if photon_energy_variable:
        photon_energy = float(random_photon_energy())
        while crash_line(fluence) > photon_energy:
            if fluence_variable:
                fluence = random.uniform(min_fluence,max_fluence)
            photon_energy = float(random_photon_energy())

    def randomize_atom(minval,maxval):
        binary = np.random.choice([1,1,1,1,1,1,1,1,1,0])
        if binary == 0:
            ret = 0
        if binary == 1:
            ret = random.uniform(minval,maxval)
        return ret

    if density_variable:
        density = random.uniform(min_density,max_density)
        
    if pulse_length_variable:
        pulse_length = float(run)
    
    loop = 1
    while loop:

        if hydrogen_variable:
            hydrogen_factor = randomize_atom(min_hydrogen,max_hydrogen)

        if carbon_variable:
            carbon_factor = randomize_atom(min_carbon,max_carbon)

        if nitrogen_variable:
            nitrogen_factor = randomize_atom(min_nitrogen,max_nitrogen)

        if oxygen_variable:
            oxygen_factor = randomize_atom(min_oxygen,max_oxygen)

        if sulphur_variable:
            sulphur_factor = randomize_atom(min_sulphur,max_sulphur)  

        if 0 in [hydrogen_factor,carbon_factor,nitrogen_factor,oxygen_factor,sulphur_factor]:
            loop = 0
        
    #Normalize the atomic factors to sum up to 100%
    atoms = [hydrogen_factor,oxygen_factor,nitrogen_factor,carbon_factor,sulphur_factor]
    normalizing_factor = sum(atoms)

    hydrogen = hydrogen_factor/normalizing_factor
    carbon = carbon_factor/normalizing_factor
    nitrogen = nitrogen_factor/normalizing_factor
    oxygen = oxygen_factor/normalizing_factor
    sulphur = sulphur_factor/normalizing_factor

    #Update and copy simulation files
    os.system(f'cp template.gen ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_FLUENCE>/{fluence}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_BEAM>/{photon_energy}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_RHO>/{density}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_PULSE_LENGTH>/{pulse_length}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_HYDROGEN>/{hydrogen}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_CARBON>/{carbon}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_NITROGEN>/{nitrogen}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_OXYGEN>/{oxygen}/g\' ./{directory}/run_{run}.gen')
    os.system(f'sed -i \'s/<MODIFY_SULPHUR>/{sulphur}/g\' ./{directory}/run_{run}.gen')


    if SASE:
        #generate SASE pulse
        make_SASE(run,pulse_length,directory)
        os.system(f'sed -i \'s/c <SASE>//g\' ./{directory}/run_{run}.gen')
        os.system(f'sed -i \'s/<MODIFY_XFILENAME>/run_{run}/g\' ./{directory}/run_{run}.gen')

        I_laser = (fluence*1.e+7)/(pulse_length*1.e-15)
        data = open(f'./runs_{start}-{end}/run_{run}.xfile', 'r')
        profile = data.read()
        data.close()
        t = []
        power = []
        for i in profile.split("\n"):
            if len(i) != 0:
                if i[0] == 't':
                    j=i.split(' ')
                    t.append(float(j[1]))
                    power.append(I_laser*float(j[2]))
                    t_fs = np.array(t)*1e-15
                    I_pulse = np.array(power)/1e7 #divide by 1e7 to convert erg->J
                    pulse_fluence = np.trapz(I_pulse,t_fs)
                    
    if gaussian:
        os.system(f'sed -i \'s/c <NO_GAUSSIAN>//g\' ./{directory}/run_{run}.gen')

    if SASE:
        keys = ['fluence','photon_energy','density','pulse_length','water_percent','hydrogen_factor','carbon_factor','nitrogen_factor','oxygen_factor','sulphur_factor','pulse_fluence']
        values = [fluence,photon_energy,density,pulse_length,water_percent,hydrogen,carbon,nitrogen,oxygen,sulphur,pulse_fluence]
    else:
        keys = ['fluence','photon_energy','density','pulse_length','water_percent','hydrogen_factor','carbon_factor','nitrogen_factor','oxygen_factor','sulphur_factor']
        values = [fluence,photon_energy,density,pulse_length,water_percent,hydrogen,carbon,nitrogen,oxygen,sulphur]

    #Generate a csv to store the parameter values
    dict = [{key:value for (key,value) in zip(keys,values)}]
    fields = keys
    with open(f'./{directory}/run_{run}_parameters.csv','w',newline='') as file:
        writer = csv.DictWriter(file,fieldnames=fields)
        writer.writeheader()
        writer.writerows(dict)

os.system(f'cp simulation_settings.txt ./runs_{start}-{end}')

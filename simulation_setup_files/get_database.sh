#!/home/harrya/.conda/envs/proj/bin/python
#SBATCH --job-name=vdb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p regular
#SBATCH --time=0:40:00


#Import packages
import os
import h5py
import numpy as np
import pandas as pd
import csv
pulse_shift = 30
HPEV = 4.135667516E-15  # Planks constant [eV s]
ERG2J = 1E-7

#['TIME', 'SPENERGY', 'ENERGY']
#[['TEVvsTIME', 'TIVvsTIME', 'ZBARvsTIME', 'ERADvsTIME', 'EDENSvsTIME', 'IONDENSITYvsTIME'], ['INTENSITYvsSPENERGY'], ['INTENSITYvsENERGY']]

add_SASE_pulse = 1

lexikon = {'TEVvsTIME':'electron_temperature',
           'TIVvsTIME':'ion_temperature',
           'ZBARvsTIME':'total_zbar',
           'ZBAR1vsTIME':'hydrogen_zbar',
           'ZBAR2vsTIME':'carbon_zbar',
           'ZBAR3vsTIME':'nitrogen_zbar',
           'ZBAR4vsTIME':'oxygen_zbar',
           'ZBAR5vsTIME':'sulphur_zbar',
           'ERADvsTIME':'radiation_energy_density',
           'EDENSvsTIME':'electron_density',
           'IONDENSITYvsTIME':'ion_density',
           'TIME':'time',
           'SPENERGY':'spectrum_energy',
           'ENERGY':'avg_energy',
           'INTENSITYvsSPENERGY':'spectrum',
           'INTENSITYvsENERGY':'time_average_spectrum',
           'YISOFRAC1vsTIME' : 'yisofrac_1',
           'YISOFRAC2vsTIME' : 'yisofrac_2',
           'YISOFRAC3vsTIME' : 'yisofrac_3',
           'YISOFRAC4vsTIME' : 'yisofrac_4',
           'YISOFRAC5vsTIME' : 'yisofrac_5',
           'YISOFRAC5vsTIME+' : 'yisofrac_5'
            }
special_cases = ['YISOFRAC1vsTIME','YISOFRAC2vsTIME','YISOFRAC3vsTIME','YISOFRAC4vsTIME','YISOFRAC5vsTIME','YISOFRAC5vsTIME+']
special_labels = [lexikon[i] for i in special_cases]
#helper function when extracting data from .plt        
def get_line(line):
    #print(line[1:])
    temp = list(line[1:].replace('  ', ',').replace('E', 'e').split(","))
    #print(temp)
    
    #temp = list(line[1:].replace('  ', ',').replace('E', 'e').split(","))
    timestep = temp[0]
    #print(temp[1:])
    if 'NaN' in temp:
        datavalue = 1
        nanval = 1
    else:
        temp = [np.float_(i) for i in temp]
        datavalue = np.mean(temp[1:])
        nanval = 0
    return [datavalue,timestep],nanval

#helper function when extracting yiso data from .plt        
def get_special_line(line):
    temp = list(np.float_(line[1:].replace('  ', ',').replace('E', 'e').split(",")))
    #print(temp)
    datavalue = temp[1:]
    return datavalue

#the yisofrac is printed in the plt file columnwise and must be treated separately
def get_yisofrac_data(start,end,run):
    special_lex = {}
    if os.path.isfile(f'/home/harrya/Master_project/simulations/data/generator_3/runs_{start}-{end}/run_{run}.plt'):
        lines = open(f'/home/harrya/Master_project/simulations/data/generator_3/runs_{start}-{end}/run_{run}.plt')
        data = []
        for i in lines:
            if i[0] == '#':
                if len(data) != 0:
                    special_lex[header] = np.array(data)
                content = ''.join(c for c in i if (c.isalpha() or c.isnumeric())).split('vs') #get names of dep/indep var
                header = f'{content[0]}vs{content[1]}'
                while header in special_lex:
                    header = header + '+'
                data = []
            if i[0] == '$':
                this_line = ''.join(c for c in i if (c.isalpha() or c.isnumeric())).split('vs') #get names of dep/indep var
            if header in special_cases and (i[0] not in ['#','$']):
                data.append(get_special_line(i))
        #append the last entry:
        if 'header' in locals():
            special_lex[header] = np.array(data)

            yisofrac5 = []
            for i in list(range(special_lex['YISOFRAC5vsTIME'].shape[1])):
                yisofrac5.append(special_lex['YISOFRAC5vsTIME'][:,i])
            for i in list(range(special_lex['YISOFRAC5vsTIME+'].shape[1])):
                yisofrac5.append(special_lex['YISOFRAC5vsTIME+'][:,i])
            yisofrac5 = np.array(yisofrac5).transpose()

            del special_lex['YISOFRAC5vsTIME']
            del special_lex['YISOFRAC5vsTIME+']
            special_lex['YISOFRAC5vsTIME'] = yisofrac5
        return special_lex
    
def make_datafiles(start,end,add_time_matrix):
#extract the time data from edits - note that every property is the average over the different sample zones
#define a by-hand-modified dictionary to map the "standard .plt" names to a more interpretable name
#Get data from .plt files
    for run in list(range(start,end+1)):
        #print(f'run: {run}')
        if os.path.isfile(f'./runs_{start}-{end}/run_{run}.plt'):
        #if True:
            if not os.path.getsize(f'./runs_{start}-{end}/run_{run}.plt') == 0:
            #if True:
                if not os.path.isfile(f'./database/files/data_{run}.hdf5'):
                    print(f'run: {run}')
                    lines = open(f'./runs_{start}-{end}/run_{run}.plt')
                    dependent_labels = []
                    independent_labels = []
                    dependent_data = []
                    independent_data = []
                    temp_dep = []
                    temp_indep = []
                    for i in lines:
                        if i[0] == '#':
                            if len(temp_dep) != 0:
                                dependent_data[independent_labels.index(current_independent)].append(temp_dep)
                                temp_dep = []
                            if len(temp_indep) != 0:
                                independent_data.append(temp_indep)
                                temp_indep = []
                            #REMOVE THIS LINE!
                            i = i.replace('VS','vs')
                            content = ''.join(c for c in i if (c.isalpha() or c.isnumeric())).split('vs') #get names of dep/indep var
                            current_independent = content[1]
                            current_dependent = f'{content[0]}vs{content[1]}'
                            if current_independent not in independent_labels:
                                independent_labels.append(current_independent)                                     #add new indep var
                                dependent_labels.append([])
                                dependent_data.append([])
                                add_independent = 1
                            else:
                                add_independent = 0

                            dependent_labels[independent_labels.index(current_independent)].append(current_dependent)     #add new dep var 

                        elif not i[0] == '$':
    #                        print(f'run {run}, {current_dependent}')
                            numbers,nanval = get_line(i)

                            temp_dep.append(numbers[0])
                            if add_independent:
                                temp_indep.append(numbers[1])

                    if len(temp_dep) != 0:
                        dependent_data[independent_labels.index(current_independent)].append(temp_dep)
                    if len(temp_indep) != 0:
                        independent_data.append(temp_indep)

                #change names
                    independent_labels = [lexikon[i] for i in independent_labels]
                    dependent_labels = [[lexikon[i] for i in j] for j in dependent_labels]

                    if not os.path.isfile(f'./database/files/data_{run}.hdf5'):
                        #print('The hdf5 of this run does not exist')
                        if len(independent_labels) != 0:
                            #print('there is data to write')
                            if not nanval:
                                #print('the ')
                                hf = h5py.File(f'./database/files/data_{run}.hdf5','w')
                                hf.close()
                                hf = h5py.File(f'./database/files/data_{run}.hdf5','r+')
                                hf.create_group(f'/data_{run}')
                                for idx in list(range(0,len(independent_labels))):
                                    hf.create_group(f'/data_{run}/{independent_labels[idx]}')
                                    hf[f'/data_{run}/{independent_labels[idx]}'].create_dataset(name=independent_labels[idx],data=independent_data[idx])
                                    for sub_idx in list(range(0,len(dependent_labels[idx]))):
                                        if dependent_labels[idx][sub_idx] not in special_labels:
                                            #print(np.array(dependent_data[idx][sub_idx]).shape)
                                            hf[f'/data_{run}/{independent_labels[idx]}'].create_dataset(name=dependent_labels[idx][sub_idx],data=dependent_data[idx][sub_idx])

                                #add attributes
                                parameters = pd.read_csv(f'./runs_{start}-{end}/run_{run}_parameters.csv')
                                attribute_names = list(parameters.columns)
                                attributes = [float(parameters[i]) for i in attribute_names]
                                for idx in list(range(0,len(attributes))):
                                    hf[f'/data_{run}'].attrs[attribute_names[idx]] = attributes[idx]

                        #add pulse from .xfile
                                if os.path.isfile(f'./runs_{start}-{end}/run_{run}.xfile'):
                                    fluence = float(hf[f'/data_{run}'].attrs['fluence'])
                                    pulse_length = float(hf[f'/data_{run}'].attrs['pulse_length'])
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
                                    #calculate the pulse_fluence [J/m2] by integrating the pulse
                        #            t_fs = np.array(t)*1e-15
                                    I_pulse = np.array(power)
                        #            pulse_fluence = np.trapz(I_pulse,t_fs)
                        #            hf[f'/data_{run}'].attrs['pulse_fluence'] = pulse_fluence
                                    hf.create_group(f'/data_{run}/pulse')
                                    hf[f'/data_{run}/pulse'].create_dataset('time',data=np.array(t))
                                    hf[f'/data_{run}/pulse'].create_dataset('intensity',data=np.array(I_pulse))

                                #Add metadata (here: units)
                                hf[f'/data_{run}'].attrs['units'] = """fluence  photon_energy  pulse_length  pulse_fluence,
                        ---------------------------------------------------,
                        J/cm2        g/cm3          fs           J/cm2"""

                                if add_time_matrix:
                        #make sure run was not terminated (only way I can work around the weird terminated cases where ntimes exsit, but there are less than ntimes timesteps)
                                    h5filename = f'./runs_{start}-{end}/run_{run}.d00'
                                    h5file = h5py.File(h5filename,'r')
                                    keys = list(h5file.keys())
                                    if 'ntimes' in keys:
                                        if "time_"+str(np.array(h5file['/ntimes'])-1) in list(h5file.keys()):
                                        # gathering spectra and simulation time
                                            ntimes = np.array(h5file['/ntimes'])
                                            evsp = np.array(h5file['/evsp'])
                                            jsp = np.zeros((evsp.size, ntimes))
                                            for t in range(0, ntimes-1):
                                                jsp[:, t] = np.mean(np.array(h5file["/time_"+str(t)+"/jsp"]), axis=0)
                                            jsp[jsp == 0.0] = 1E-40  # log doesn't like zeros
                                            spectrum_time_matrix = jsp
                                            hf.create_group(f'/data_{run}/gothere')
                                            hf[f'/data_{run}/time'].create_dataset('spectrum_time_matrix',data=spectrum_time_matrix)

                                #Add yisofrac:
                                yisofrac_lex = get_yisofrac_data(start,end,run)
                                for yiso_label in list(yisofrac_lex.keys()):
                                    hf[f'/data_{run}/time'].create_dataset(lexikon[yiso_label],data=yisofrac_lex[yiso_label])

                                hf.close()

start = 0
end = 1
add_time_matrix = 0
if not os.path.isdir('./database'):
    os.system('mkdir database')
    os.system('mkdir database/files')
os.system('cp create_virtual_database.py ./database')
os.system('cp functions.py ./database')


make_datafiles(start,end,add_time_matrix)

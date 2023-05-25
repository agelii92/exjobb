#Import packages
import os
import h5py
from matplotlib import pyplot as plt
from matplotlib import colors
import numpy as np
import pandas as pd
import csv
HPEV = 4.135667516E-15  # Planks constant [eV s]
ERG2J = 1E-7
import random
from scipy import stats


#return data specified by name 'dataset' of hdf5 file named 'h5filename'
def h5read(h5filename, dataset):
    with h5py.File(h5filename, "r") as infile:
        return np.array(infile[dataset])

    
#extract data from .d00 file and return as array objects
def get_simulation_data(run,start,end):
    run = str(run)
    h5filename = f'./runs_{start}-{end}/run_{run}.d00'
    # gathering spectra and simulation time
    ntimes = h5read(h5filename, "/ntimes")      #number of time step
    stime = h5read(h5filename, '/times')        #time values, one value for each time step
    evsp = h5read(h5filename, "/evsp")          #photon energy values, kind of the x-axis when plotting spectrum
    eav = h5read(h5filename, "/eav")            #values of the "ebins 1e-1 1e4 500" energy values on a exp scale (?) used for making calculation of continuum radiation more efficient?
    #jnu = np.zeros((eav.size, ntimes))          #continuum radiation intensity
    jsp = np.zeros((evsp.size, ntimes))         #spectral radiation intensity

    for t in range(0, ntimes-1):
    #    jnu[:, t] = np.mean(h5read(h5filename, "/time_"+str(t)+"/jnu"), axis=0) #fill up the time step values of emitted continuum radiation. returns array of size 500: each value is the average emitted intensity for a specific energy in eav, average taken over the 4 zones
        jsp[:, t] = np.mean(h5read(h5filename, "/time_"+str(t)+"/jsp"), axis=0) ##fill up the time step values of emitted spectral radiation. returns array of size 10000: each value is the average emitted intensity for a specific energy in evsp, average taken over the 4 zones

    # ergs/cm^2/sec/Hz/str -> W/cm^2  #Change units accordingly
    E1 = 2000  # [eV]
    E2 = E1 * 1.001  # [eV]
    mult = HPEV / (E2 - E1)
    #jnu = (4*np.pi * ERG2J * jnu)/mult
    jsp = (4*np.pi * ERG2J * jsp)/mult
    fstime = stime * 1E15  # s -> fs #change units accordingly

    #jnu[jnu == 0.0] = 1E-40  # log doesn't like zeros       #prevent singular behaviour
    jsp[jsp == 0.0] = 1E-40  # log doesn't like zeros'
    int_jsp = np.trapz(jsp, stime, axis=1)                  #integrate (not just sum!) over the time steps to obtain a accumulative spectrum (again, both spectral and continuum part)
    #int_jnu = np.trapz(jnu, stime, axis=1)
    parameters = pd.read_csv(f'./runs_{start}-{end}/run_{run}_parameters.csv')

    attribute_names = list(parameters.columns)
    attributes = [float(parameters[i]) for i in attribute_names]

    return ntimes, fstime, evsp, int_jsp, attribute_names, attributes


#rearange the data from get_simulation_data into instances of .hdf5 files
def add_to_database(database, run,start,end):
    run=str(run)
    ntimes, fstime, evsp, int_jsp, attribute_names, attributes = get_simulation_data(run,start,end)

    data_group = database.create_group('/data_'+run)
    for i in list(range(0,len(attributes))):
        data_group.attrs[attribute_names[i]] = attributes[i]
 
    energy_spectrum_data = data_group.create_dataset('emission_spectrum',data=int_jsp)
    energy_spectrum_data.attrs['unit'] = 'W/cm^2' 
    energy_bins = data_group.create_dataset('energy_bins',data=evsp)
    energy_bins.attrs['unit'] = 'eV'
    energy_spectrum_data.attrs['readme'] = 'Cretin simulated emission spectra after interraction with X-ray.'
        
    #energy_spectrum_time_development = data_group.create_dataset('energy_spectrum_time_development',data=jsp)
    #energy_spectrum_time_development.attrs['number_of_timesteps'] = ntimes
    #energy_spectrum_time_development.attrs['time'] =fstime
    #energy_spectrum_time_development.attrs['unit'] = 'fs'
    #energy_spectrum_time_development.attrs['readme'] = 'The average time development of the energy distribution of emission. Shown here as a matrix of size (n_timesteps x n_energy_bins). Average is with respect to the different measuring zones within the sample.'

    #last security check:
    delete = False
    for i in list(np.array(database['/data_'+run]['emission_spectrum'])):
        if np.isnan(i):
            print('found NaN')
            delete = True
    if delete:
        del database['/data_'+run]
        print('deleted')

    database.close
#create the .hdf5 files    
def make_datafiles():
    start = input('start run')
    end = input('end run')
    if not os.path.isdir('./database'):
        os.system('mkdir database')
        os.system('mkdir database/files')
    os.system('cp create_virtual_database.py ./database')
    os.system('cp functions.py ./database')
    runs = list(range(int(start),int(end)+1))
    for run in runs:
        if os.path.isfile(f'./runs_{start}-{end}/run_{run}.d00'):
            if 'ntimes' in list(h5py.File(f'./runs_{start}-{end}/run_{run}.d00', "r").keys()): #this prevents trying to add simulations that have been aborted: such simulations will not have 'ntimes' property stored.
                if 'time_'+str(np.array(h5py.File(f'./runs_{start}-{end}/run_{run}.d00', "r")['/ntimes'])-1) in list(h5py.File(f'./runs_{start}-{end}/run_{run}.d00', "r").keys()): #think about it... badly written but ony way I find to stop error
                    if not os.path.isfile(f'./database/files/data_{run}'):
                        hf=h5py.File(f'./database/files/data_{run}.hdf5', 'w')
                        add_to_database(hf,run,start,end)
                        hf.close()
                        
#create references of each hdf5 file to a virtual database                        
def make_virtual_database():
    hf = h5py.File('VDB.hdf5','w')
    hf.close()
    hf=h5py.File('VDB.hdf5','r+')
    for run in list(range(0,2000)):
        source_filename = f'./files/data_{run}.hdf5'
        if os.path.isfile(source_filename):    
            src = h5py.File(source_filename,'r')
            hf.create_group(f'/data_{run}')
            for attribute in list(src[f'/data_{run}'].attrs):
                hf[f'/data_{run}'].attrs[attribute] = src[f'/data_{run}'].attrs[attribute] #should be OK up to here
            for groups in list(src[f'/data_{run}'].keys()):
                hf.create_group(f'/data_{run}/{groups}')
                for sets in list(src[f'/data_{run}/{groups}'].keys()):
                    datatype = src[f'/data_{run}/{groups}/{sets}'].dtype
                    datashape = src[f'/data_{run}/{groups}/{sets}'].shape
                    layout = h5py.VirtualLayout(shape=datashape, dtype=datatype) # Virtual target is a representation of the output dataset
                    vsource = h5py.VirtualSource(source_filename, f'/data_{run}/{groups}/{sets}', shape=datashape)    #   Sources represent the input datasets
                    layout[:] = vsource                                         # Map the inputs into the virtual dataset

                    hf[f'/data_{run}/{groups}'].create_virtual_dataset(sets, layout)
    hf.close()

    
def make_SASE(run,pulse_length,directory):
    #make pulse
    L = random.uniform(5,50)
    maxval = random.uniform(5e9,5e10)
    L=pulse_length
    maxval = 5e10
    pulse_shift = 30
    start = -L/2
    end = L/2
    time = np.linspace(start, end, int((end-start)/0.04))
    y = stats.norm.pdf(time, 0, L/6)
    y = np.abs([y[i]*np.random.uniform(0, 2)*np.random.choice([random.uniform(0.7,1.3),random.uniform(0.1,0.3),random.uniform(0.1,0.3),random.uniform(0.1,0.3)]) for i in list(range(0,len(time)))])
    y = y/max(y)
    
    #make .dat
    with open(f'./{directory}/run_{run}.dat', 'w') as f:
        f.write('# x = Time (fs); y = Power (W)')
        f.write('\n')
        for i in list(range(0,len(time))):
            f.write(f'{time[i]} {y[i]}')
            f.write('\n')
    #read and edit info from .dat
    data = open(f'./{directory}/run_{run}.dat', 'r')
    profile = data.read()
    data.close()
    t = []
    power = []
    for line in profile.split("\n"):
        if line.startswith("#") or len(line) == 0:
            continue
        split_line = line.split()
        t.append(float(split_line[0])+pulse_shift)
        power.append(float(split_line[1]))
            
    #make .xfile
    xfile = open(f"./{directory}/run_{run}.xfile", 'w')
    source = f"source jbndry 1 E1 E2 value history 1 MULT\n"
    pulse = f"history 1 ILASER 1e-15\n"
    xfile.writelines([source, pulse])
    for index in range(1, len(time)):
        tvline = f"tv {t[index-1]:.4f} {power[index-1]:.4f}\n"
        xfile.write(tvline)
    xfile.close()
    return


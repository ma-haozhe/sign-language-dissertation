################
## Imports
################
import mne
import os
import glob
import numpy as np
import pandas as pd
from scipy.io import savemat
import matplotlib.pyplot as plt

################
## Parameters
################
# General
#root_dir = './inputs/EEG_raw/'
#output_base = './outputs/64Hz/EEG/'
#root_dir = './eegProject/datasets/SLdata/eeg/'
#output_base = './eegProject/datasets/SLdata'
#June 19, 24: let's switch back to use EEG_raw folder, more generalized
root_dir = './EEG_raw/EEG/'
output_base = './outputs_new/64Hz/EEG'



csv_sequence = 'stimuli_sequence.csv'
csv_behave = 'answers_to_questions.csv'
csv_pre = 'preprocessing_pipeline.csv'
plot = False
#subjects_to_process = ['698908', '752086']  # ID of the subject
#subjects_to_process = ['164123', '319885'] - these two are done
subjects_to_process = ['698908', 
                       '992123', '888521', '966706', 
                       '910810', '780931', '437658', 
                       '176622']

# Preprocessing
# Notch filtering
notch_applied = False
freq_notch = 50

# Bandpass filtering 
bpf_applied = True
bandpass = '1-30Hz'
freq_low   = 1
freq_high  = 30
ftype = 'butter'
order = 3

# Spherical interpolation
int_applied = False
interpolation = 'spline'

# Rereferencing using average of mastoids electrodes
reref_applied = True
reref_type = 'Mastoids'  #Mastoids/Average

# Downsampling
down_applied = True
downfreq = 64

############################################################################
## Loop over .bdf recordings
############################################################################
# Find BioSemi files in root_dir and itrate over them
files = glob.glob(os.path.join(root_dir, '**', '*.bdf'), recursive=True)
# files = [file for file in files if subjects_to_process.any() in file]
for idx, file in enumerate(files):
     # Extract the subject ID from the file name
    subject_id = file.split('/')[-1].split('.')[0]
    if subject_id not in subjects_to_process:
        continue
    # if file.split('/')[-1].split('.')[0] not in subjects_to_process:
    #     continue
    print(idx, file.split('/')[-1])

    # Select which file to open
    # idx_file_to_open = input('Enter idx of subject to open: ')
    # file_to_open = files[int(idx_file_to_open)]
    file_to_open = file

    # Create output folder
    name_subject = file_to_open.split('/')[-1].split('.')[0]
    output_dir = os.path.join(output_base, bandpass, reref_type, name_subject)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # Create file that keeps track of the preprocessing
    df_pre = pd.DataFrame()
    print("file to open: ", file_to_open)
    ##################################################################

    # # Load and resave stimuli sequence
    # path = file_to_open.split('/')
    # print("path: ", path)
    # path[-1] = csv_sequence
    # print("path for csv_sequence: ", path[-1])
    # path = '/'.join(path)
    # print("path for csv_sequence-2: ", path)
    # df_stimuli = pd.read_csv(path)  
    # print(df_stimuli)
    # df_stimuli.to_csv(os.path.join(output_dir, csv_sequence), index=False)

    # # Load and resave behavioural questions
    # path = file_to_open.split('/')
    # path[-1] = csv_behave
    # path = '/'.join(path)
    # if os.path.exists(path):
    #     df_behave = pd.read_csv(path)  
    #     df_behave.to_csv(os.path.join(output_dir, csv_behave), index=False)

    ##################################################################
    # Construct the path for the subject-specific folder
    subject_folder_pattern = f'{subject_id}_*'
    print("subject_folder_pattern: ", subject_folder_pattern)
    print("root_dir: ", root_dir)
    subject_folder = glob.glob(os.path.join(root_dir, subject_folder_pattern))
    print("subject_folder: ", subject_folder)
    if subject_folder:
        subject_folder = subject_folder[0]

        # Load and resave stimuli sequence
        stim_path = os.path.join(subject_folder, csv_sequence)
        if os.path.exists(stim_path):
            df_stimuli = pd.read_csv(stim_path)  
            print("we do have df_stimuli: ", df_stimuli)
            df_stimuli.to_csv(os.path.join(output_dir, csv_sequence), index=False)
        else:
            print(f"Stimuli sequence file not found for subject {subject_id}")

        # Load and resave behavioural questions
        behave_path = os.path.join(subject_folder, csv_behave)
        if os.path.exists(behave_path):
            df_behave = pd.read_csv(behave_path)  
            df_behave.to_csv(os.path.join(output_dir, csv_behave), index=False)
        else:
            print(f"Behavioural questions file not found for subject {subject_id}")
    else:
        print(f"Subject folder not found for subject {subject_id}")
    ############################################################################
    ## Load EEG data
    ############################################################################
    # Load EEG data
    raw = mne.io.read_raw_bdf(file_to_open, eog=["LO1", "LO2", "IO1", "IO2"], stim_channel='auto', 
                            misc=['SO1', 'SO2', 'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp'], 
                            infer_types=False, preload=False, verbose=None)
    raw.info.set_montage('standard_1020', match_case=False)
    print(raw)
    raw.load_data()

    # Check metadata
    n_time_samps = raw.n_times
    time_secs = raw.times
    ch_names = raw.ch_names
    n_chan = len(ch_names) 
    print('the (cropped) sample data object has {} time samples and {} channels.'
        ''.format(n_time_samps, n_chan))
    print('The last time sample is at {} seconds.'.format(time_secs[-1]))
    print('The first few channel names are {}.'.format(', '.join(ch_names[:3])))
    print('bad channels:', raw.info['bads'])  # chs marked "bad" during acquisition
    print(raw.info['sfreq'], 'Hz')            # sampling frequency
    print(raw.info['description'], '\n')      # miscellaneous acquisition info
    print(raw.info)
    if plot:
        raw.plot(start=100, duration=10)

    ############################################################################
    ## Find events in trigger channel
    ############################################################################
    # Find starting sample of events in the trigger channel
    N_start_events = mne.find_events(raw)
    print(N_start_events)

    # Select only relevant events
    N_start_events = N_start_events[N_start_events[:, 2] == 65281]

    # Get starting samples of events
    N_start_events = N_start_events[:, 0]

    # Get corresponding time
    T_start_events = N_start_events / raw.info['sfreq']

    # Check there are as many triggers as the listened stimuli
    assert len(N_start_events) == len(df_stimuli)

    ############################################################################
    ## Segment and then process EEG data
    ############################################################################
    
    t_max = 0
    for index, row_stimuli in df_stimuli.iterrows():

        ###############################
        ## Crop data
        ###############################
        # Get the name of the listened stimulus
        filename = row_stimuli['File'].split('.')[0].split('/')[-1]
        
        # Get the onset and duration in seconds and in samples
        onset_seconds = T_start_events[index] 
        assert onset_seconds > t_max

        duration_seconds = row_stimuli['Time']
        assert duration_seconds > 0

        t_min = onset_seconds
        t_max = onset_seconds + duration_seconds
        eeg = raw.copy().crop(tmin=t_min, tmax=t_max)

        ###############################
        ## Preprocessing
        ###############################
        ## -------------
        ## Select channels
        ## -------------
        eeg_channels = ch_names[0:66]
        eeg = eeg.pick_channels(eeg_channels)
        if plot:
            eeg.plot(start=100, duration=10, n_channels=len(raw.ch_names))

        ## -------------
        ## Notch filtering
        ## -------------
        df_pre['notch_applied'] = [notch_applied]
        if notch_applied:
            eeg = eeg.notch_filter(freqs=freq_notch)
            df_pre['notch'] = [freq_notch]
            if plot:
                eeg.plot()

        ## -------------
        ## BPFiltering
        ## -------------
        df_pre['bpf_applied'] = [bpf_applied]
        if bpf_applied:
            iir_params = dict(order=order, ftype=ftype)
            filter_params = mne.filter.create_filter(eeg.get_data(), eeg.info['sfreq'], 
                                                    l_freq=freq_low, h_freq=freq_high, 
                                                    method='iir', iir_params=iir_params)

            if plot:
                flim = (1., eeg.info['sfreq'] / 2.)  # frequencies
                dlim = (-0.001, 0.001)  # delays
                kwargs = dict(flim=flim, dlim=dlim)
                mne.viz.plot_filter(filter_params, eeg.info['sfreq'], compensate=True, **kwargs)
                # plt.savefig(os.path.join(output_dir, 'bpf_ffilt_shape.png'))

            eeg = eeg.filter(l_freq=freq_low, h_freq=freq_high, method='iir', iir_params=iir_params)
            df_pre['bandpass'] = [iir_params]
            df_pre['HPF'] = [freq_low]
            df_pre['LPF'] = [freq_high]
            if plot:
                eeg.plot()

        ## -------------
        ## Intrpolation
        ## -------------
        df_pre['int_applied'] = [int_applied]
        if int_applied: 
            eeg = eeg.interpolate_bads(reset_bads=False)  #, method=interpolation

            # Get the indices and names of the interpolated channels
            interp_inds = eeg.info['bads']
            interp_names = [eeg.info['ch_names'][i] for i in interp_inds]

            # Print the number and names of the interpolated channels
            print(f'{len(interp_inds)} channels interpolated: {interp_names}')

            df_pre['interpolation'] = [interpolation]
            df_pre['interp_inds'] = [interp_inds]
            df_pre['interp_names'] = [interp_names]

            if plot:
                eeg.plot()

        ## -------------
        ## Rereferencing
        ## -------------
        df_pre['reref_applied'] = [reref_applied]
        if reref_applied:
            # Set electrodes for rereferencing
            if reref_type == 'Mastoids':
                reref_channels = ['M1', 'M2']   
            else:
                reref_channels = 'average'           

            # Actually r-referencing signals
            eeg = eeg.set_eeg_reference(ref_channels=reref_channels)
            df_pre['reref_type'] = [reref_type]
            df_pre['reref_channels'] = [reref_channels]
            if plot:
                eeg.plot()

        ## -------------
        ## Resampling
        ## -------------
        df_pre['down_applied'] = [down_applied]
        if down_applied:
            eeg = eeg.resample(sfreq=downfreq)
            df_pre['downfreq'] = [downfreq]
            print(eeg.info)
            if plot:
                eeg.plot()

        ## -------------
        ## Save preprocessing stages
        ## -------------
        df_pre.to_csv(os.path.join(output_dir, csv_pre), index=False)

        # Get only data matrix
        eeg = eeg.get_data()
        
        # Save file to mupy in the subject's folder
        print('Saving EEG responses to ', filename, eeg.shape)
        savemat(os.path.join(output_dir, filename + '.mat'), {'trial_data': eeg[0:64, :], 
                                                              'trial_mastoids': eeg[64:, :]})
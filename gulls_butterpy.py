import pandas as pd
import numpy as np
import butterpy as bp

import time

from scipy import stats
import glob
from joblib import Parallel, delayed

from astropy.modeling.models import BlackBody
from astropy import units as u

class GullsButterpyInterface:
    def __init__(self, input_lightcurve_path, output_directory,ndays=5000,chi2_threshold=5000, var_fraction = 0.25):
        self.parameters_array = None # initialize but fill later

        self.renormalized_stellar_flux = None
        self.input_lightcurve_path = input_lightcurve_path
        self.output_directory = output_directory
        self.ndays = ndays # number of days to simulate over
        self.days_offset = None
        #columns = ['Simulation_time', 'measured_relative_flux', 'measured_relative_flux_error', 'true_relative_flux',
        #           'true_relative_flux_error', 'observatory_code', 'saturation_flag', 'best_single_lens_fit',
        #           'parallax_shift_t', 'parallax_shift_u',
        #           'BJD', 'source_x', 'source_y', 'lens1_x', 'lens1_y', 'lens2_x', 'lens2_y', 'satX', 'satY', 'satZ']

        self.gulls_dataframe = pd.read_csv(input_lightcurve_path, sep='\s+', comment='#',na_filter=False,low_memory=False)
        with open(input_lightcurve_path, 'r') as f:
            fs1_string = f.readline()
            fs2_string = f.readline()
            f.readline()
            f.readline()
            source1data_string = f.readline()
            source2data_string = f.readline()
        fs1_string = fs1_string.split(' ')
        fs2_string = fs2_string.split(' ')

        fs1_list = []
        fs2_list = []
        for i in range(1, len(fs1_string) - 1):
            fs1_list.append(float(fs1_string[i]))
            fs2_list.append(float(fs2_string[i]))
        self.fs1_list = fs1_list
        self.fs2_list = fs2_list
        source1data_string = source1data_string.split(' ')
        self.source1_teff = float(source1data_string[11])

        source2data_string = source2data_string.split(' ')
        self.source2_teff = float(source2data_string[11])

        if fs2_list[0] == 0.:
            #print("Event is single source")
            self.N_sources = 1
        else:
            #print("Event is binary source")
            self.N_sources = 2


        self.chi2_threshold = chi2_threshold
        self.realization_done = False  # initialize to none so we can check if the realization was created
        self.multiband_done = False


        #check to see if we will add variabilty to this  with a random check against var_fraction
        if stats.uniform.rvs(0, 1, size=1)[0] <= var_fraction:
            self.use_var_columns = True
        else:
            self.use_var_columns = False

    def draw_parameters(self):
        parameters_array = []
        parameters_array.append(stats.loguniform.rvs(0.01,50, size=1)[0]) #activity level log_uniform 1-10
        parameters_array.append(stats.loguniform.rvs(1, 40, size=1)[0]) #cycle period log uniform 1-40 yr
        parameters_array.append(stats.loguniform.rvs(0.1, parameters_array[1], size=1)[0])  # cycle overlap log_uniform 0.1 - T_cycle
        parameters_array.append(stats.uniform.rvs(0, 40, size=1)[0]) # minlat uniform 0-40
        parameters_array.append(stats.uniform.rvs(parameters_array[3]+5, 80, size=1)[0])  # maxlat
        parameters_array.append(stats.uniform.rvs(0.1, 180, size=1)[0]) # period in days
        sin2i = stats.uniform.rvs(0, 1, size=1)[0]
        i = np.arcsin(np.sqrt(sin2i))*180/np.pi#0-90 uniform in sin2i
        parameters_array.append(i)
        parameters_array.append(stats.loguniform.rvs(1,10, size=1)[0]) # spot liftime
        choice_index = np.random.choice([0,1,2],size=1,p=[0.5,0.25,0.25])
        if choice_index ==0:
            parameters_array.append(stats.loguniform.rvs(0.1, 1,size = 1)[0])
        elif choice_index == 1:
            parameters_array.append(0)
        else:
            parameters_array.append(-stats.loguniform.rvs(0.1, 1, size=1)[0])
        parameters_array.append(self.source1_teff)
        parameters_array.append(self.source1_teff-stats.uniform.rvs(100., 1000., size=1)[0])
        self.parameters_array = parameters_array
        return parameters_array

    def create_butterpy_realization(self, days_offset=None):
        t_min = self.gulls_dataframe['BJD'].min()
        t_max = self.gulls_dataframe['BJD'].max()
        delta_t = t_max - t_min
        times_array = self.gulls_dataframe['BJD'] - t_min
        if days_offset:
            self.days_offset = days_offset
        else:
            self.days_offset = np.random.randint(500, self.ndays - np.ceil(
                delta_t))  # place GBDTS survey at some random time after start of LC
        self.star = bp.Surface()
        self.star.emerge_regions(
            ndays=self.ndays,
            butterfly=True,
            activity_level=self.parameters_array[0],
            cycle_period=self.parameters_array[1],
            cycle_overlap=self.parameters_array[2],
            max_lat=self.parameters_array[4],
            min_lat=self.parameters_array[3],
        );
        self.spot_lc = self.star.evolve_spots(
            times_array + self.days_offset,
            # leave as is
            alpha_med=1e-4,
            threshold=0.1,
            period=self.parameters_array[5],
            inclination=self.parameters_array[6],
            tau_evol=self.parameters_array[7],  # leave as is
            shear=self.parameters_array[8]  # leave as is
        );
        # get mean survey flux and renormalize to that
        self.mean_flux = np.mean(self.spot_lc.flux)
        self.renormalized_flux = self.spot_lc.flux / self.mean_flux
        self.time = self.spot_lc.time
        self.realization_done = True
        return None
    def create_full_butterpy_lightcurve(self):
        time = np.arange(0, self.ndays, 1)
        spot_lc_full = self.star.evolve_spots(
            time,
            alpha_med=1e-4,  # leave as is
            period=self.parameters_array[5],
            inclination=self.parameters_array[6],
            tau_evol=10,  # leave as is
            shear=0.25  # leave as is
        );
        self.renormalized_flux_full = spot_lc_full.flux / self.mean_flux
        self.time_full = spot_lc_full.time

    def apply_blackbody_to_spot_lc(self, band_waves=[1.464, 0.869, 2.125]):
        '''
        Make the variable lightcurves band dependent.
        Based on a code Robby Wilson shared with me.
        '''
        if not self.realization_done:
            raise ValueError('you must first call create_butterpy_realization')
        self.spot_fraction = 1 - self.spot_lc.flux # treat diff from 1 fraction of surface covered in spots
        spot_blackbody = BlackBody(temperature=self.parameters_array[10]* u.Kelvin)
        star_blackbody = BlackBody(temperature=self.parameters_array[9]* u.Kelvin)
        cond_list = []
        for i in range(len(self.fs1_list)):
            cond_list.append(self.gulls_dataframe['observatory_code'] == i)
        flux_fraction = (spot_blackbody(np.array(band_waves) * u.um).value / star_blackbody(np.array(band_waves) * u.um).value)
        flux_fraction_array = np.select(cond_list, flux_fraction, default=-1)

        # get spot fluxes
        spot_fluxes = self.spot_fraction * flux_fraction_array
        #get difference from mean
        self.renormalized_stellar_flux = np.zeros(self.spot_lc.flux.shape)
        for i in range(len(self.fs1_list)):
            filter_indices = np.where(self.gulls_dataframe['observatory_code'] == i)[0]
            sp_flx = spot_fluxes[filter_indices]
            var_lc = self.spot_lc.flux[filter_indices]
            band_flux = (sp_flx + var_lc) / np.mean(sp_flx + var_lc)
            self.renormalized_stellar_flux[filter_indices] = band_flux

        self.multiband_done = True


    def apply_renormalized_flux_to_event(self):
        if not self.multiband_done:
            raise ValueError('you must first call apply_blackbody_to_spot_lc')
        cond_list = []

        for i in range(len(self.fs1_list)):
            cond_list.append(self.gulls_dataframe['observatory_code'] == i)
        fs1_array = np.select(cond_list, self.fs1_list, default=-1)

        self.gulls_dataframe['variability_flux'] = self.renormalized_stellar_flux
        self.gulls_dataframe['variability_flux_bol'] = self.spot_lc.flux

        mu_true = (self.gulls_dataframe['true_relative_flux'] - 1 + fs1_array) / fs1_array


        self.gulls_dataframe['true_relative_flux_var'] = fs1_array * self.renormalized_stellar_flux * mu_true + (1 - fs1_array)
        self.gulls_dataframe['true_relative_flux_var_error'] = self.renormalized_stellar_flux * self.gulls_dataframe[
            'true_relative_flux_error']

        fractional_shift = (self.gulls_dataframe['measured_relative_flux'] - self.gulls_dataframe['true_relative_flux']) / self.gulls_dataframe['true_relative_flux']
        self.gulls_dataframe['measured_relative_flux_var'] = (1 + fractional_shift) * self.gulls_dataframe['true_relative_flux_var']

        fractional_error = self.gulls_dataframe['measured_relative_flux_error'] / self.gulls_dataframe['measured_relative_flux']
        self.gulls_dataframe['measured_relative_flux_var_error'] = fractional_error * self.gulls_dataframe['measured_relative_flux_var']

        return None


    def apply_renormalized_flux_to_event_2source(self):
        if not self.multiband_done:
            raise ValueError('you must first call apply_blackbody_to_spot_lc')
        cond_list = []
        for i in range(len(self.fs1_list)):
            cond_list.append(self.gulls_dataframe['observatory_code'] == i)
        fs1_array = np.select(cond_list, self.fs1_list, default=-1)
        fs2_array = np.select(cond_list, self.fs2_list, default=-1)


        source1_flux = self.gulls_dataframe['source1_relative_flux']
        source2_flux = self.gulls_dataframe['source2_relative_flux']
        self.gulls_dataframe['variability_flux'] = self.renormalized_stellar_flux
        self.gulls_dataframe['variability_flux_bol'] = self.renormalized_flux


        self.gulls_dataframe['true_relative_flux_var'] = source1_flux * self.renormalized_stellar_flux + source2_flux + (1 - fs1_array - fs2_array)
        self.gulls_dataframe['true_relative_flux_var_error'] = self.renormalized_stellar_flux * self.gulls_dataframe['true_relative_flux_error']

        fractional_shift = (self.gulls_dataframe['measured_relative_flux'] - self.gulls_dataframe['true_relative_flux']) / self.gulls_dataframe['true_relative_flux']
        self.gulls_dataframe['measured_relative_flux_var'] = (1 + fractional_shift) * self.gulls_dataframe['true_relative_flux_var']

        fractional_error = self.gulls_dataframe['measured_relative_flux_error'] / self.gulls_dataframe['measured_relative_flux']
        self.gulls_dataframe['measured_relative_flux_var_error'] = fractional_error * self.gulls_dataframe['measured_relative_flux_var']

        #self.gulls_dataframe['variability_flux'] = self.renormalized_stellar_flux
        #self.gulls_dataframe['variability_flux_bol'] = self.spot_lc.flux
        return None

    def get_delta_chi2(self):
        chi2_novar = np.sum(((self.gulls_dataframe['true_relative_flux'] - self.gulls_dataframe['measured_relative_flux']) /
                             self.gulls_dataframe['measured_relative_flux_error']) ** 2)
        chi2_var = np.sum(((self.gulls_dataframe['true_relative_flux'] - self.gulls_dataframe['measured_relative_flux_var']) /
                           self.gulls_dataframe['measured_relative_flux_error']) ** 2)

        delta_chi2 = chi2_var - chi2_novar
        return delta_chi2

    def set_vars_tonan(self):
        '''
        Utility function for filling variability related columns with NaNs.
        Used for consistency of output format.
        :return: None
        '''
        self.gulls_dataframe['variability_flux'] = np.nan
        self.gulls_dataframe['variability_flux_bol'] = np.nan
        self.gulls_dataframe['true_relative_flux_var'] = np.nan
        self.gulls_dataframe['true_relative_flux_var_error'] = np.nan
        self.gulls_dataframe['measured_relative_flux_var']=np.nan
        self.gulls_dataframe['measured_relative_flux_var_error']=np.nan
        #print('This happened')

        return None

    def write_lightcurve_file(self, nlines=14):
        output_lightcurve_name = f'{self.input_lightcurve_path.split("/")[-1].split(".")[0]}_var.det.lc'
        output_lightcurve_path = f'{self.output_directory}/{output_lightcurve_name}'
        # if os.path.exists(output_lightcurve_path):
        #    raise FileExistsError('Your lightcurve already exists in this location')


        # If this LC was not selected for variability, just write nans into variability columns.
        if self.use_var_columns:
            #if LC was selected, check chi2. If over threshold, overwrite vars with nans as well.
            delta_chi2 = self.get_delta_chi2()
            if delta_chi2 > self.chi2_threshold:
                self.use_var_columns = False
                self.set_vars_tonan()
        else:
            self.set_vars_tonan()

        #Now get ready to write the output file.
        # copy lightcurve header
        with open(output_lightcurve_path, 'w') as out:
            with open(self.input_lightcurve_path, 'r') as inp:
                for i in range(nlines):
                    inline = inp.readline()
                    out.write(inline)
            varstring = '#butterpy variables:'
            for par in self.parameters_array:
                varstring += f'{par} '
            varstring += f'{self.ndays} {self.days_offset}\n'
            out.write(varstring)
            usevar_string = f'#usevar:{self.use_var_columns}\n'
            out.write(usevar_string)
            #debugstr = f'#{self.spot_lc.flux.dtypes} {self.gulls_dataframe["variability_flux_bol"].dtypes}\n'
            #out.write(debugstr)
        self.gulls_dataframe.to_csv(output_lightcurve_path, sep=' ', mode='a',na_rep = 'nan',index=False)
        return None




def process_lc(input_path,output_path):

    time0 = time.time()
    #print(input_path)
    ndays = 5000
    g2bp = GullsButterpyInterface(input_lightcurve_path=input_path,
                                  output_directory=output_path,ndays=ndays,var_fraction=0.5)
    #parameters_array = [2, 5, 2, 45, 20, 12, 75]
    parameters_array = g2bp.draw_parameters()#draw parameters no matter what for consistency
    #only simulate the lc if it was selected to be variable
    if g2bp.use_var_columns:
        g2bp.create_butterpy_realization()
        g2bp.apply_blackbody_to_spot_lc()
        if g2bp.N_sources == 1:
            g2bp.apply_renormalized_flux_to_event()
        elif g2bp.N_sources == 2:
            g2bp.apply_renormalized_flux_to_event_2source()
        else:
            raise ValueError("Invalid number of sources!")
        #g2bp.create_full_butterpy_lightcurve(parameters_array)
    g2bp.write_lightcurve_file()
    runtime = time.time() - time0
    return runtime


if __name__ == "__main__":
    input_path = '/Users/jmbrashe/Data_Challenge/for_Jon'
    output_path = '../small_test_outputs'
    lightcurves = glob.glob(input_path + '/*.lc')
    output_generator = Parallel(n_jobs=5,verbose=1)(delayed(process_lc)(lightcurves[i],output_path) for i in range(len(lightcurves)))
    times = np.array(output_generator)
    print(np.median(times),np.std(times))
    np.savetxt(f'{output_path}/times.txt',times)
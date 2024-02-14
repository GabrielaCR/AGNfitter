

"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PLOTandWRITE_AGNfitter.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This script contains all functions used in order to visualize the output of the sampling.
Plotting and writing. 
These functions need the chains saved in files samples_mcmc.sav and samples_bur-in.sav.
This script includes:

- main function
- class OUTPUT
- class CHAIN
- class FLUXESARRAYS
- functions SED_plotting_settings, SED_colors

"""
#PYTHON IMPORTS
import matplotlib.pyplot as plt
from matplotlib import rc, ticker
import matplotlib.ticker
import os
import math 
import numpy as np
from . import corner #Author: Dan Foreman-Mackey (danfm@nyu.edu)
import scipy
from astropy import units as u
import pandas as pd 


#AGNfitter IMPORTS
from . import MODEL_AGNfitter as model
from . import PARAMETERSPACE_AGNfitter as parspace
import pickle



def main(data, models, P, out, models_settings, mc_settings):

    """
    Main function of PLOTandWRITE_AGNfitter. Output depends of settings in RUN_AGNfitter.

    ##input:
    - data object
    - dictionary P (parameter space settings)
    - dictionary out (output settings)
    
    """
    if mc_settings['sampling_algorithm'] == 'emcee':
        chain_burnin = CHAIN(data.output_folder+str(data.name)+ '/samples_burn1-2-3.sav', data.output_folder+str(data.name)+ '/samples_burn1-2-3.sav', out, mc_settings)
        chain_mcmc = CHAIN(data.output_folder+str(data.name)+ '/samples_mcmc.sav', data.output_folder+str(data.name)+ '/samples_mcmc.sav', out, mc_settings)
        chain_mcmc.props()
        print( '_________________________________')
        print( 'Properties of the sampling results:')
        print( '- Mean acceptance fraction', chain_mcmc.mean_accept)
        print( '- Mean autocorrelation time', chain_mcmc.mean_autocorr)

    elif mc_settings['sampling_algorithm'] == 'ultranest':
        chain_mcmc = CHAIN(data.output_folder+str(data.name)+ '/ultranest/chains/weighted_post.txt', data.output_folder+str(data.name)+ '/ultranest/chains/run.txt',  out, mc_settings) #equal_weighted_post.txt
    else:
        print('Unknown algorithm, please select one of the available: ultranest or emcee')

    output = OUTPUT(chain_mcmc, data, models, P, models_settings, mc_settings)

    if out['plot_posteriortriangle'] :
        fig = output.plot_PDFtriangle('10pars', P['names'])
        fig.savefig(data.output_folder+str(data.name)+'/PDFtriangle_10pars.' + out['plot_format'])
        plt.close(fig)

    if out['plot_posteriortrianglewithluminosities']: 
        fig = output.plot_PDFtriangle('int_lums', out['intlum_names']) 
        fig.savefig(data.output_folder+str(data.name)+'/PDFtriangle_intlums.' + out['plot_format'])
        plt.close(fig)

    if out['save_posterior_luminosities']: 
        int_lums_array =output.save_realizations()
        np.savetxt(data.output_folder+str(data.name)+'/posterior_intlums_'+str(data.name)+'.txt', int_lums_array.T, fmt= "%1.4f" ,header= str(out['intlum_names']))

    if out['save_posteriors']: 
        posteriors,posteriors_header = output.save_posteriors(P)
        np.savetxt(data.output_folder + str(data.name)+'/posteriors_'+str(data.name)+'.txt' , posteriors, delimiter = " ",fmt= "%1.4f" ,header= posteriors_header)

    if out['writepar_meanwitherrors']:
        outputvalues, outputvalues_header = output.write_parameters_outputvalues(P)
        comments_ouput= ' # Output for source ' +str(data.name) + '\n' +' Rows are: 2.5, 16, 50, 84, 97.5 percentiles # '+'\n'+ '-----------------------------------------------------'+'\n' 
        np.savetxt(data.output_folder + str(data.name)+'/parameter_outvalues_'+str(data.name)+'.txt' , outputvalues, delimiter = " ",fmt= "%1.4f" ,header= outputvalues_header, comments =comments_ouput)

    if ((out['saveSEDrealizations']) or (out['plotSEDrealizations']) or (out['saveSEDresiduals'])):
        fig, save_SEDS, save_residuals = output.plot_manyrealizations_SED(plot_residuals=out['plot_residuals'])
        n_rp=out['realizations2plot']

        if models.settings['RADIO'] == True:
            SEDs_header = '#freq '+' '.join(['SBnuLnu'+str(i) for i in range(n_rp)]) +' ' +' '.join(['BBnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['GAnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['RADnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOTALnuLnu'+str(i) for i in range(n_rp)]) +' '+' '.join(['BBnuLnu_deredd'+str(i) for i in range(n_rp)]) 
        else:
            SEDs_header = '#freq '+' '.join(['SBnuLnu'+str(i) for i in range(n_rp)]) +' ' +' '.join(['BBnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['GAnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOTALnuLnu'+str(i) for i in range(n_rp)]) +' '+' '.join(['BBnuLnu_deredd'+str(i) for i in range(n_rp)]) 

        if out['saveSEDrealizations']:
            np.savetxt(data.output_folder + str(data.name)+'/output_SEDs100_'+str(data.name)+'.txt' , save_SEDS, delimiter = " ",fmt= "%1.4f" ,header= SEDs_header, comments='')
        if out['saveSEDresiduals']:
            res_header = '#freq res'
            np.savetxt(data.output_folder + str(data.name)+'/output_residuals100_'+str(data.name)+'.txt' , save_residuals, delimiter = " ",fmt= "%1.4f" ,header= res_header, comments='')
        if out['plotSEDrealizations']:
            fig.savefig(data.output_folder+str(data.name)+'/SED_manyrealizations_' +str(data.name)+ '.'+out['plot_format'], bbox_inches='tight', dpi=250)

        plt.close(fig)
  
"""=========================================================="""


class OUTPUT:

    """
    Class OUTPUT

    Includes the functions that return all output products.
    You can call all chain and data properties, since it inherits chain and data classes.

    ##input: 
    - object of the CHAIN class, object of DATA class
    """    

    def __init__(self, chain_obj, data_obj, model_obj, P, models_settings, mc_settings):

        self.chain = chain_obj
        if mc_settings['sampling_algorithm'] == 'emcee':
            self.chain.props()
   
        self.out = chain_obj.out
        self.data = data_obj
        self.models = model_obj
        self.models_settings = models_settings
        self.mc_settings = mc_settings
        self.z=self.data.z
        fluxobj_withintlums = FLUXES_ARRAYS(chain_obj, P,  self.out,'int_lums', self.models_settings, self.mc_settings)
        fluxobj_4SEDplots = FLUXES_ARRAYS(chain_obj, P, self.out,'plot', self.models_settings, self.mc_settings)

        fluxobj_withintlums.fluxes( self.data, self.models)
        self.nuLnus = fluxobj_withintlums.nuLnus4plotting
        self.allnus = fluxobj_withintlums.all_nus_rest
        self.int_lums = fluxobj_withintlums.int_lums
        self.int_lums_best = fluxobj_withintlums.int_lums_best

        if ((self.out['plotSEDrealizations']) or (self.out['saveSEDrealizations']) or (self.out['saveSEDresiduals'])):
            fluxobj_4SEDplots.fluxes(self.data, self.models)
            self.nuLnus = fluxobj_4SEDplots.nuLnus4plotting
            self.filtered_modelpoints_nuLnu = fluxobj_4SEDplots.filtered_modelpoints_nuLnu
            self.allnus = fluxobj_4SEDplots.all_nus_rest

    def write_parameters_outputvalues(self, P):        

        Mstar, SFR_opt = model.stellar_info_array(self.chain.flatchain, self.data, self.models, self.out['realizations2int'], self.mc_settings)
        column_names = np.transpose(np.array(["P025","P16","P50","P84","P975"], dtype='|S3'))
        chain_pars = np.column_stack((self.chain.flatchain_sorted, Mstar, SFR_opt))              

        SFR_IR = model.sfr_IR(self.int_lums[0]) #check that ['intlum_names'][0] is always L_IR(8-100)               

        chain_others =np.column_stack((self.int_lums.T, SFR_IR))        
        if self.mc_settings['sampling_algorithm'] == 'ultranest':
            outputvalues = np.column_stack((np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_pars, [2.5,16, 50, 84,97.5], axis=0))))), np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_others, [2.5,16, 50, 84,97.5], axis=0))))),
                                        np.transpose(np.percentile(self.chain.lnprob_flat, [2.5,16, 50, 84,97.5], axis=0)),
                                        np.transpose(np.percentile(self.chain.logz_flat, [2.5,16, 50, 84,97.5], axis=0)) )) 
            outputvalues_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', 'log_like', 'logz'))] ) 

        elif self.mc_settings['sampling_algorithm'] == 'emcee':
            outputvalues = np.column_stack((np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_pars, [2.5,16, 50, 84,97.5], axis=0))))), np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_others, [2.5,16, 50, 84,97.5], axis=0))))),
                                        np.transpose(np.percentile(self.chain.lnprob_flat, [2.5,16, 50, 84,97.5], axis=0)) ))  

            outputvalues_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', '-ln_like'))] )

        return outputvalues, outputvalues_header

    def save_posteriors(self, P):   
        Mstar, SFR_opt = model.stellar_info_array(self.chain.flatchain, self.data, self.models, self.out['realizations2int'], self.mc_settings)
        column_names = np.transpose(np.array(["P025","P16","P50","P84","P975"], dtype='|S3'))
        chain_pars = np.column_stack((self.chain.flatchain, Mstar, SFR_opt))      
        
        if self.out['calc_intlum']:            

            SFR_IR = model.sfr_IR(self.int_lums[0]) #check that ['intlum_names'][0] is always L_IR(8-100)        

            chain_others =np.column_stack((self.int_lums.T, SFR_IR))
            if self.mc_settings['sampling_algorithm'] == 'ultranest':
                outputvalues = np.column_stack((np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_pars, [2.5,16, 50, 84,97.5], axis=0))))), np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_others, [2.5,16, 50, 84,97.5], axis=0))))), np.transpose(np.percentile(self.chain.lnprob_flat, [2.5,16, 50, 84,97.5], axis=0)), 
                                            np.transpose(np.percentile(self.chain.logz_flat, [2.5,16, 50, 84,97.5], axis=0)) ))  
   
                outputvalues_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', 'log_like', 'logz'))] )
            elif self.mc_settings['sampling_algorithm'] == 'emcee':
                outputvalues = np.column_stack((np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_pars, [2.5,16, 50, 84,97.5], axis=0))))), np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_others, [2.5,16, 50, 84,97.5], axis=0))))), np.transpose(np.percentile(self.chain.lnprob_flat, [2.5,16, 50, 84,97.5], axis=0)) ))  
   
                outputvalues_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', '-ln_like'))] )

            if self.out['save_posteriors']:  
                nsample, npar = self.chain.flatchain.shape
                if self.mc_settings['sampling_algorithm'] == 'ultranest':
                    posteriors = np.column_stack((chain_pars[np.random.choice(nsample, (self.out['realizations2int']))], chain_others, self.chain.lnprob_flat.iloc[np.random.choice(nsample, (self.out['realizations2int']))],
self.chain.logz_flat.iloc[np.random.choice(nsample, (self.out['realizations2int']))])) 
                    posteriors_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', 'log_like', 'logz'))] )
                elif self.mc_settings['sampling_algorithm'] == 'emcee':
                    posteriors = np.column_stack((chain_pars[np.random.choice(nsample, (self.out['realizations2int'])),:], chain_others, self.chain.lnprob_flat[np.random.choice(nsample, (self.out['realizations2int']))] )) 
                    posteriors_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', '-ln_like'))] )
                
                return posteriors,posteriors_header
        else:
            print('Error: save_posteriors=True requires calc_intlum=True.')

    def save_realizations(self):        
        return self.int_lums

    def plot_PDFtriangle(self,parameterset, labels):        

        if parameterset=='10pars':
            figure = corner.corner(self.chain.flatchain,levels=[0.68,0.95],  labels= labels, plot_contours=True, plot_datapoints = False, show_titles=True, quantiles=[0.16, 0.50, 0.84])
        elif parameterset == 'int_lums':
            figure = corner.corner(self.int_lums.T, levels=[0.68,0.95], labels= labels,   plot_contours=True, plot_datapoints = False, show_titles=True, quantiles=[0.16, 0.50, 0.84])
        return figure


    def plot_manyrealizations_SED(self,plot_residuals=True):    

        #reading from valid data from object data
        ydata = self.data.fluxes[self.data.fluxes>0.]
        yerror = self.data.fluxerrs[self.data.fluxes>0.]
        yndflags = self.data.ndflag[self.data.fluxes>0.]
        Nrealizations = self.out['realizations2plot']

        #Data frequencies (obs and rest), and model frequencies
        data_nus_obs = 10**self.data.nus[self.data.fluxes>0.]
        data_nus_rest = data_nus_obs * (1+self.z) 
        data_nus = np.log10(data_nus_rest)

        all_nus =self.allnus
        all_nus_rest = 10**all_nus 
        all_nus_obs =  10**all_nus / (1+self.z) #observed

        distance= model.z2Dlum(self.z)
        lumfactor = (4. * math.pi * distance**2.)
        data_nuLnu_rest = ydata* data_nus_obs *lumfactor
        data_errors_rest= yerror * data_nus_obs * lumfactor

        if self.models_settings['RADIO'] == True:
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, RADnuLnu, TOTALnuLnu, BBnuLnu_deredd = self.nuLnus
            save_SEDs= np.column_stack((all_nus,SBnuLnu.T,BBnuLnu.T, GAnuLnu.T, TOnuLnu.T, RADnuLnu.T, TOTALnuLnu.T, BBnuLnu_deredd.T  ))
        else:
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, TOTALnuLnu, BBnuLnu_deredd = self.nuLnus
            save_SEDs= np.column_stack((all_nus,SBnuLnu.T,BBnuLnu.T, GAnuLnu.T, TOnuLnu.T, TOTALnuLnu.T, BBnuLnu_deredd.T  ))

        #plotting settings
        if plot_residuals:
            fig, ax1, ax2, axr = SED_plotting_settings(all_nus_rest, data_nuLnu_rest, self.allnus, self.models, self.out, True)
        else:
            fig, ax1, ax2 = SED_plotting_settings(all_nus_rest, data_nuLnu_rest, self.allnus, self.models, self.out)
        SBcolor, BBcolor, GAcolor, TOcolor, RADcolor, TOTALcolor= SED_colors(combination = 'a')


        lw= 2 #1.5


        alp = 0.25
        mec='None'
        area_sb = np.zeros((len(SBnuLnu[0]),3))
        area_bb = np.zeros((len(BBnuLnu[0]),3))
        area_ga = np.zeros((len(GAnuLnu[0]),3))
        area_to = np.zeros((len(TOnuLnu[0]),3))
        if self.models_settings['RADIO'] == True:
            area_rad = np.zeros((len(RADnuLnu[0]),3))
        area_total = np.zeros((len(TOTALnuLnu[0]),3))
        if Nrealizations == 1:
            alp = 1.0

        if self.out['band_indicators'] == True:
            if self.models_settings['RADIO'] == True:
                ax1.axvspan(9, 11.5, 0.965, 1, color='gray', alpha=0.15, lw=0) #Radio
            ax1.axvspan(11.5, 12.47, 0.965, 1, color='gray', alpha=0.6, lw=0) #FIR
            ax1.axvspan(12.47, 13.48, 0.965, 1, color='gray', alpha=0.15, lw=0) #MIR
            ax1.axvspan(13.48, 14.68, 0.965, 1, color='gray', alpha=0.6, lw=0) #NIR
            ax1.axvspan(14.68, 14.87, 0.965, 1, color='gray', alpha=0.15, lw=0) #OPT
            ax1.axvspan(14.87, 16.48, 0.965, 1, color='gray', alpha=0.6, lw=0) #UV
            if self.models_settings['XRAYS'] != False:
                ax1.axvspan(16.48, 19.48,0.965, 1, color='gray', alpha=0.15, lw=0) #x-ray
            
        for i in range(Nrealizations):
            # last one is the max likelihood fit  
            if i == Nrealizations-1:
                alp = 1
                lw = 2
                mec='k'
                save_residuals= np.column_stack((data_nus, np.array(data_nuLnu_rest-self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.])/data_errors_rest))
                if self.out['plot_best_fit'] == True:
                    p2=ax1.plot(all_nus, SBnuLnu[i], marker="None", linewidth=lw, label="Starburst", color= SBcolor, alpha = alp)
                    p4=ax1.plot(all_nus, GAnuLnu[i],marker="None", linewidth=lw, label="Stellar population",color=GAcolor, alpha = alp)
                    if self.models_settings['turn_on_AGN'] == True:
                        p3=ax1.plot(all_nus[all_nus < 16.671], BBnuLnu[i][all_nus < 16.671], marker="None", linewidth=lw, label="Accretion disk",color= BBcolor, alpha = alp)
                        p3xr=ax1.plot(all_nus[all_nus > 16.685], BBnuLnu[i][all_nus > 16.685], marker="None", linewidth=lw,color= BBcolor, alpha = alp)
                        p5=ax1.plot( all_nus, TOnuLnu[i], marker="None",  linewidth=lw, label="Torus",color= TOcolor ,alpha = alp)                    

                        if self.models_settings['RADIO'] == True:
                            p6=ax1.plot( all_nus, RADnuLnu[i], marker="None",  linewidth=lw, label="Radio from AGN",color= RADcolor ,alpha = alp)
                    p1=ax1.plot( all_nus[all_nus < 16.671], TOTALnuLnu[i][all_nus < 16.671], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp)
                    p1xr=ax1.plot( all_nus[all_nus > 16.685], TOTALnuLnu[i][all_nus > 16.685], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp)
                    ax1.legend(fontsize = 13, shadow = True, loc = 6) #NEW

                if self.out['realizations_format'] == 'shaded':
                    for k in range(len(SBnuLnu[0])):
                        area_sb[k] = SBnuLnu[:, k].min(), SBnuLnu[:, k].mean(), SBnuLnu[:, k].max()
                    p2 = ax1.fill_between(all_nus, area_sb[:, 0], area_sb[:, 2], color= SBcolor, alpha= 0.3)  

                    for k in range(len(GAnuLnu[0])):
                        area_ga[k] = GAnuLnu[:, k].min(), GAnuLnu[:, k].mean(), GAnuLnu[:, k].max()
                    p4 = ax1.fill_between(all_nus, area_ga[:, 0], area_ga[:, 2], color= GAcolor, alpha= 0.3) 

                    if self.models_settings['turn_on_AGN'] == True:
                        for k in range(len(BBnuLnu[0])):
                            area_bb[k] = BBnuLnu[:, k].min(), BBnuLnu[:, k].mean(), BBnuLnu[:, k].max()
                        p3 = ax1.fill_between(all_nus, area_bb[:, 0], area_bb[:, 2], color= BBcolor, alpha= 0.3)  

                        for k in range(len(TOnuLnu[0])):
                            area_to[k] = TOnuLnu[:, k].min(), TOnuLnu[:, k].mean(), TOnuLnu[:, k].max()
                        p5 = ax1.fill_between(all_nus, area_to[:, 0], area_to[:, 2], color= TOcolor, alpha= 0.3)  

                        if self.models_settings['RADIO'] == True:
                            for k in range(len(RADnuLnu[0])):
                                area_rad[k] = RADnuLnu[:, k].min(), RADnuLnu[:, k].mean(), RADnuLnu[:, k].max()
                            p6 = ax1.fill_between(all_nus, area_rad[:, 0], area_rad[:, 2], color= RADcolor, alpha= 0.3)  

                    for k in range(len(TOTALnuLnu[0])):
                        area_total[k] = TOTALnuLnu[:, k].min(), TOTALnuLnu[:, k].mean(), TOTALnuLnu[:, k].max()
                    #p1 = ax1.fill_between(all_nus, area_total[:, 0], area_total[:, 2], color= TOTALcolor, alpha= 0.3)
  
                    if self.out['plot_median_fit'] == True:
                        p2=ax1.plot(all_nus, area_sb[:, 1], marker="None", linewidth=lw, label="Starburst", color= SBcolor, alpha = alp, linestyle = '-')
                        p4=ax1.plot(all_nus, area_ga[:,1], marker="None", linewidth=lw, label="Stellar population",color=GAcolor, alpha = alp, linestyle = '-')
                        if self.models_settings['turn_on_AGN'] == True:
                            p3=ax1.plot(all_nus[all_nus < 16.671], area_bb[:, 1][all_nus < 16.671], marker="None", linewidth=lw, label="Accretion disk",color= BBcolor, alpha = alp, linestyle = '-')
                            p3xr=ax1.plot(all_nus[all_nus > 16.685], area_bb[:, 1][all_nus > 16.685], marker="None", linewidth=lw,color= BBcolor, alpha = alp, linestyle = '-')
                            p5=ax1.plot( all_nus, area_to[:,1], marker="None",  linewidth=lw, label="Torus",color= TOcolor ,alpha = alp, linestyle = '-')
                            if self.models_settings['RADIO'] == True:
                                p6=ax1.plot( all_nus, area_rad[:, 1], marker="None",  linewidth=lw, label="Radio from AGN",color= RADcolor ,alpha = alp, linestyle = '-')
                        p1=ax1.plot( all_nus[all_nus < 16.671], area_total[:, 1][all_nus < 16.671], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp, linestyle = '-')
                        p1xr=ax1.plot( all_nus[all_nus > 16.685], area_total[:, 1][all_nus > 16.685], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp, linestyle = '-')
                        ax1.legend(fontsize = 13, shadow = True, loc = 6) #NEW
 
            if i != Nrealizations -1 and self.out['realizations_format'] == 'curves':
                p2=ax1.plot(all_nus, SBnuLnu[i], marker="None", linewidth=lw,  color= SBcolor, alpha = alp)
                p3=ax1.plot(all_nus, BBnuLnu[i], marker="None", linewidth=lw,color= BBcolor, alpha = alp)
                cut_nu_xr = [all_nus[j] for j in range(len(BBnuLnu[i])-2) if (np.log10(BBnuLnu[i][j+1])-np.log10(BBnuLnu[i][j]) > 2) and (np.log10(BBnuLnu[i][j+2])-np.log10(BBnuLnu[i][j+1]) < 0.5) and (all_nus[j] > 15)]  #16.685
                cut_nu_bb = [all_nus[j] for j in range(len(BBnuLnu[i])-1) if (np.log10(BBnuLnu[i][j])-np.log10(BBnuLnu[i][j+1]) > 0.4) and (all_nus[j] > 15) and (all_nus[j] < 17)] #16.671
                p3=ax1.plot(all_nus[all_nus < cut_nu_bb], BBnuLnu[i][all_nus < cut_nu_bb], marker="None", linewidth=lw,color= BBcolor, alpha = alp) 
                p3xr=ax1.plot(all_nus[all_nus > cut_nu_xr], BBnuLnu[i][all_nus > cut_nu_xr], marker="None", linewidth=lw,color= BBcolor, alpha = alp)
                p4=ax1.plot(all_nus, GAnuLnu[i],marker="None", linewidth=lw, color=GAcolor, alpha = alp) 
                p5=ax1.plot( all_nus, TOnuLnu[i], marker="None",  linewidth=lw, color= TOcolor ,alpha = alp)

                if self.models_settings['RADIO'] == True and i == Nrealizations -1:
                    p6=ax1.plot( all_nus, RADnuLnu[i], marker="None",  linewidth=lw, color= RADcolor ,alpha = alp)

                #p1=ax1.plot( all_nus, TOTALnuLnu[i], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp)
                p1=ax1.plot( all_nus[all_nus < 16.671], TOTALnuLnu[i][all_nus < 16.671], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp)
                p1xr=ax1.plot( all_nus[all_nus > 16.685], TOTALnuLnu[i][all_nus > 16.685], marker="None", linewidth=lw,  color= TOTALcolor, alpha= alp)

            #If some Xray prior was applied, plot the Xray emission according to the UV or midIR flux from models in each realization
            if self.models_settings['XRAYS'] == 'Prior_UV':
                bbb_nuLnu_2500Angs = BBnuLnu_deredd[i][(15.04 < all_nus) & (all_nus < 15.15 )][0]
                bbb_flux_2500Angs = bbb_nuLnu_2500Angs/(10**(all_nus[BBnuLnu_deredd[i] == bbb_nuLnu_2500Angs])) 
                Fnu_2kev = 10**(0.6430*np.log10(bbb_flux_2500Angs) + 6.8734)                   #UV-Xray correlation by Just et al. 2007

                #Proportionality constant a to scale x-ray power-law in 2keV to the value found with alpha_OX-L_2500
                h = 4.135667731*1e-15*1e-3                                                     #eV/Hz --> keV/Hz
                nu_2kev = 4.83598*1e17 
                a = Fnu_2kev/((h*nu_2kev)**(-1.8+1)*np.e**(-nu_2kev/(7.2540*1e19)))
                xray_nu = np.logspace(16.685, 19.7, 1000)                                          #with a hole between BB template and X-Rays
                xray_Fnu = a*(h*xray_nu)**(-1.8+1)*np.e**(-xray_nu/(7.2540*1e19))
                xray_nuLnu = xray_Fnu*xray_nu
                p8=ax1.plot(np.log10(np.logspace(16.685, 19.7, 1000)), xray_nuLnu, marker="None", linewidth=lw, linestyle = '--', color= BBcolor, alpha = alp) 

            if self.models_settings['XRAYS'] == 'Prior_midIR':
                to_nuLnu_6microns = TOnuLnu[i][(13.59897 < all_nus) & (all_nus < 13.79897)][0]
                x = np.log10(to_nuLnu_6microns/1e41)
                L2_10keV_model = 10**(40.981 + 1.024*x - 0.047*x**2)                            #midIR-Xray correlation by Stern 2015
                F2_10keV = L2_10keV_model/(10**17.9061)

                #Proportionality constant a to scale x-ray power-law in 2keV to the value found with alpha_OX-L_2500
                h = 4.135667731*1e-15*1e-3                                                      #eV/Hz --> keV/Hz
                nu_2_10kev = 10**17.9061
                a = F2_10keV/((h*nu_2_10kev)**(-1.8+1)*np.e**(-nu_2_10kev/(7.2540*1e19)))
                xray_nu = np.logspace(16.685, 19.7, 1000)                                           #with a hole between BB template and X-Rays
                xray_Fnu = a*(h*xray_nu)**(-1.8+1)*np.e**(-xray_nu/(7.2540*1e19))
                xray_nuLnu = xray_Fnu*xray_nu
                p8=ax1.plot(np.log10(np.logspace(16.685, 19.7, 1000)), xray_nuLnu, marker="None", linewidth=lw, linestyle = '--', color= TOcolor, alpha = alp)

            det = [yndflags==1]
            upp = [yndflags==0]
            det2 = [yndflags==1 & (data_nus < 15.38)| (data_nus > 16.685)]

            p7 = ax1.plot(data_nus[tuple(det2)], self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.][tuple(det2)],   marker='o', linestyle="None",markersize=5, color="red", alpha =alp)
            
            if plot_residuals: # and i == Nrealizations -1:
                p6r = axr.plot(data_nus[tuple(det2)], (data_nuLnu_rest[tuple(det2)]-self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.][tuple(det2)])/data_errors_rest[tuple(det2)],   marker='o', mec=mec, linestyle="None",markersize=5, color="red", alpha =alp)
                #p6r = axr.plot(data_nus[tuple(det2)], ((data_nuLnu_rest[tuple(det2)]-self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.][tuple(det2)])/data_nuLnu_rest[tuple(det2)])*100,   marker='o', mec=mec, linestyle="None",markersize=5, color="red", alpha =alp)
            
        upplimits = ax1.errorbar(data_nus[tuple(upp)], 2.*data_nuLnu_rest[tuple(upp)], yerr= data_errors_rest[tuple(upp)]/2, uplims = True, linestyle='',  markersize=5, color="black")
        (_, caps, _) = ax1.errorbar(data_nus[tuple(det)], data_nuLnu_rest[tuple(det)], yerr= data_errors_rest[tuple(det)], capsize=4, linestyle="None", linewidth=2,  marker='o',markersize=5, color="black", alpha = 1)


        ax1.text(0.04, 0.8, r'id='+str(self.data.name)+r', z ='+ str(self.z), ha='left', transform=ax1.transAxes, fontsize = 19 )
        if plot_residuals:
            if self.mc_settings['sampling_algorithm'] == 'ultranest':
                ax1.text(0.96, 0.8, 'max log-likelihood = {ml:.1f}'.format(ml=np.max(self.chain.lnprob_flat)), ha='right', transform=ax1.transAxes, fontsize = 19 )
            elif self.mc_settings['sampling_algorithm'] == 'emcee':
                ax1.text(0.96, 0.8, 'max ln-likelihood = {ml:.1f}'.format(ml=np.max(self.chain.lnprob_flat)), ha='right', transform=ax1.transAxes, fontsize = 19 )
        print(' => SEDs of '+ str(Nrealizations)+' different realization were plotted.')

        return fig, save_SEDs, save_residuals


"""=========================================================="""


class CHAIN:

    """
    Class CHAIN

    ##input: 
    - name of file, where chain was saved
    - dictionary of ouput setting: out

    ##bugs: 

    """     

    def __init__(self, outputfilename, runfile, out, mc_settings):
            self.outputfilename = outputfilename
            self.runfile = runfile
            self.out = out
            self.mc_settings = mc_settings

    def props(self):
        if os.path.lexists(self.outputfilename):
            f = open(self.outputfilename, 'rb')

            if self.mc_settings['sampling_algorithm'] == 'ultranest':
                samples = pd.read_csv(f, sep = ' ')
                f.close()
                f2 = open(self.runfile, 'rb')
                runfile = pd.read_csv(f2, sep = ' ')
                f2.close()

                cumsum = np.cumsum(np.array(samples['weight']))
                mask = cumsum > 1e-4
                size_mask = int(np.sum(mask) - (np.sum(mask) % 100))     

                self.lnprob = samples.iloc[-size_mask:]
                self.chain = samples.iloc[-size_mask:, 2:]
                self.logz = runfile.iloc[-size_mask:]

                self.lnprob_flat = self.lnprob['logl'] 
                self.logz_flat = self.logz['logz']

                #sort parameter vector for likelihood
                idx_sorted = np.argsort(np.array(-self.lnprob_flat), axis=0)
                lnprob_sorted = self.lnprob.sort_values(by = ['logl'], ascending=False)['logl']
                self.lnprob_max = lnprob_sorted.iloc[0]
                self.flatchain = self.chain
                chain_length = int(len(self.flatchain))

                self.flatchain_sorted = self.flatchain.iloc[idx_sorted]
                self.best_fit_pars = self.flatchain.iloc[idx_sorted[0]].values

            elif self.mc_settings['sampling_algorithm'] == 'emcee':
                samples = pickle.load(f, encoding='latin1')
                f.close()
                self.chain = samples['chain']
                nwalkers, nsamples, npar = samples['chain'].shape

                Ns, Nt = self.out['Nsample'], self.out['Nthinning'] 
                self.lnprob = samples['lnprob']
                self.lnprob_flat = samples['lnprob'][:,int(nsamples/2):int(nsamples):Nt].ravel() #[:,0:Ns*Nt:Nt]

                isort = (- self.lnprob_flat).argsort() #sort parameter vector for likelihood
                lnprob_sorted = np.reshape(self.lnprob_flat[isort],(-1,1))
                self.lnprob_max = lnprob_sorted[0]

                self.flatchain = samples['chain'][:,int(nsamples/2):int(nsamples):Nt,:].reshape(-1, npar) #[:,0:Ns*Nt:Nt,:]
                chain_length = int(len(self.flatchain))

                self.flatchain_sorted = self.flatchain[isort]
                self.best_fit_pars = self.flatchain[isort[0]]

                self.mean_accept =  samples['accept'].mean()
                self.mean_autocorr = samples['acor'].mean()

        else:

            'Error: The sampling has not been perfomed yet, or the chains were not saved properly.'

    def plot_trace(self, P, nwplot=50):

        """ Plot the sample trace for a subset of walkers for each parameter.
        """
        #-- Latex -------------------------------------------------
        rc('text', usetex=True)
        rc('font', family='serif')
        rc('axes', linewidth=2)
        #-------------------------------------------------------------
        self.props()

        self.nwalkers, nsample, npar = self.chain.shape
        nrows = npar + 1
        ncols =1     

        def fig_axes(nrows, ncols, npar, width=13):
            fig = plt.figure(figsize=(width, nrows*1.8))#*nrows/ncols))      
            fig.subplots_adjust(hspace=0.9)
            axes = [fig.add_subplot(nrows, ncols, i+1) for i in range(npar)]
            return fig, axes

        fig, axes = fig_axes(nrows, ncols, npar+1)

        nwplot = min(nsample, nwplot)
        for i in range(npar):
            ax = axes[i]
            for j in range(0, self.nwalkers, max(1, self.nwalkers // nwplot)):
                ax.plot(self.chain[j,:,i], lw=0.5,  color = 'black', alpha = 0.3)
            ax.set_title(r'\textit{Parameter : }'+P['names'][i], fontsize=12)  
            ax.set_xlabel(r'\textit{Steps}', fontsize=12)
            ax.set_ylabel(r'\textit{Walkers}',fontsize=12)

        ax = axes[-1]
        for j in range(0, self.nwalkers, max(1, self.nwalkers // nwplot)):
            ax.plot(self.lnprob[j,:], lw=0.5, color = 'black', alpha = 0.3)
        ax.set_title(r'\textit{Likelihood}', fontsize=12)   
        ax.set_xlabel(r'\textit{Steps}', fontsize=12)
        ax.set_ylabel(r'\textit{Walkers}',fontsize=12)

        return fig, nwplot

"""=========================================================="""


class FLUXES_ARRAYS:

    """
    This class constructs the luminosities arrays for many realizations from the parameter values
    Output is returned by FLUXES_ARRAYS.fluxes().
    
    ## inputs:
    - object of class CHAIN
    - dictionary of output settings, out
    - str giving output_type: ['plot', 'intlum',  'bestfit']

    ## output:
    - frequencies and nuLnus + ['filteredpoints', 'integrated luminosities', - ]
    """


    def __init__(self, chain_obj, P, out, output_type, models_settings, mc_settings):
        self.chain_obj = chain_obj
        self.output_type = output_type
        self.out = out
        self.models_settings=models_settings
        self.mc = mc_settings
        self.P = P

    def fluxes(self, data, models):    

        """
        This is the main function of the class.
        """
        self.chain_obj.props()

        SBFnu_list = []
        BBFnu_list = []
        GAFnu_list= []
        TOFnu_list = []
        TOTALFnu_list = []
        BBFnu_deredd_list = []
        if self.output_type == 'plot':
            filtered_modelpoints_list = []

        gal_obj,sb_obj,tor_obj, bbb_obj, agnrad_obj = models.dictkey_arrays_4plot

        # Take the  4 dictionaries for plotting. Dicts are described in DICTIONARIES_AGNfitter.py
        MD= models.dict_modelfluxes
        nsample, npar = self.chain_obj.flatchain.shape
        source = data.name

        if self.output_type == 'plot':
            if self.mc['sampling_algorithm'] == 'ultranest':
                par = self.chain_obj.flatchain.iloc[np.random.choice(nsample, (self.out['realizations2plot'])),:].copy()#.T 
                par_best = self.chain_obj.best_fit_pars
                par.loc[-1]=par_best
            elif self.mc['sampling_algorithm'] == 'emcee':
                par = self.chain_obj.flatchain[np.random.choice(nsample, (self.out['realizations2plot'])),:]#.T   
                #par = self.chain_obj.flatchain[np.random.choice(np.arange(nsample-1000, nsample), (20*self.out['realizations2plot'])),:]#.T 200/1000 works    
                par_best = self.chain_obj.best_fit_pars
                par[-1]=par_best
  
            realization_nr=self.out['realizations2plot']
        
        elif self.output_type == 'int_lums':
            if self.mc['sampling_algorithm'] == 'ultranest':
                par = self.chain_obj.flatchain.iloc[np.random.choice(nsample, (self.out['realizations2int'])),:].copy()#.T
                par_best = self.chain_obj.best_fit_pars
                par.loc[-1] = par_best
            elif self.mc['sampling_algorithm'] == 'emcee':
                par = self.chain_obj.flatchain[np.random.choice(nsample, (self.out['realizations2int'])),:]#.T
                par_best = self.chain_obj.best_fit_pars
                par[-1]=par_best
            realization_nr=self.out['realizations2int']
        
        elif self.output_type == 'best_fit':
            par = self.chain_obj.best_fit_pars

        if self.models_settings['BBB'] =='D12_S' or self.models_settings['BBB'] =='D12_K' or self.models_settings['XRAYS'] != False:
            ### extend SED to X-rays if there are models with Xrays or Xray priors were applied
            lognu_max = 19
        else:
            lognu_max = 16.5

        if self.models_settings['RADIO']==True:
            ## extend SED to radio
            RADFnu_list = []
            lognu_min = 9
        else:
            lognu_min = 10.5

        self.all_nus_rest = np.arange(lognu_min, lognu_max, 0.001) 

        for g in range(realization_nr):
            ## Pick dictionary key-values, nearest to the MCMC- parameter values
            if self.mc['sampling_algorithm'] == 'ultranest':
                gal_obj.pick_nD(par.iloc[g][self.P['idxs'][0]:self.P['idxs'][1]]) 
                sb_obj.pick_nD(par.iloc[g][self.P['idxs'][1]:self.P['idxs'][2]])
                tor_obj.pick_nD(par.iloc[g][self.P['idxs'][2]:self.P['idxs'][3]]) 
                bbb_obj.pick_nD(par.iloc[g][self.P['idxs'][3]:self.P['idxs'][4]])

                if models.settings['BBB']=='R06' or models.settings['BBB']=='THB21': 
                    if models.settings['RADIO'] == True:
                        GA, SB, TO, BB, RAD = par.iloc[g][-5:]
                    else:
                        GA, SB, TO, BB = par.iloc[g][-4:]   
                else:
                    if models.settings['RADIO'] == True:
                        GA, SB, TO, RAD = par.iloc[g][-4:]
                    elif models.settings['RADIO'] == False:
                        GA, SB, TO = par.iloc[g][-3:]
                    BB = 0.
  
            elif self.mc['sampling_algorithm'] == 'emcee':
                gal_obj.pick_nD(par[g][self.P['idxs'][0]:self.P['idxs'][1]]) 
                sb_obj.pick_nD(par[g][self.P['idxs'][1]:self.P['idxs'][2]])  
                tor_obj.pick_nD(par[g][self.P['idxs'][2]:self.P['idxs'][3]])
                bbb_obj.pick_nD(par[g][self.P['idxs'][3]:self.P['idxs'][4]])

                if models.settings['BBB']=='R06' or models.settings['BBB']=='THB21': 
                    if models.settings['RADIO'] == True:
                        GA, SB, TO, BB, RAD = par[g][-5:]
                    else:
                        GA, SB, TO, BB = par[g][-4:]   
                else:
                    if models.settings['RADIO'] == True:
                        GA, SB, TO, RAD = par[g][-4:]
                    elif models.settings['RADIO'] == False:
                        GA, SB, TO = par[g][-3:]
                    BB = 0.

            all_gal_nus, gal_Fnus = gal_obj.get_fluxes(gal_obj.matched_parkeys)
            all_sb_nus, sb_Fnus= sb_obj.get_fluxes(sb_obj.matched_parkeys) 
            all_tor_nus, tor_Fnus= tor_obj.get_fluxes(tor_obj.matched_parkeys)
            all_bbb_nus, bbb_Fnus = bbb_obj.get_fluxes(bbb_obj.matched_parkeys)

            #Produce model fluxes at all_nus_rest for plotting, through interpolation

            GAinterp = scipy.interpolate.interp1d(all_gal_nus, np.log10(gal_Fnus.flatten()),  kind = 'linear', bounds_error=False, fill_value=-100.)
            all_gal_Fnus = 10**GAinterp(self.all_nus_rest)

            if models.settings['RADIO'] == True:

                if (agnrad_obj.pars_modelkeys != ['-99.9']).all() :    #If there is a radio model with fitting parameters
                    if self.mc['sampling_algorithm'] == 'ultranest':
                        agnrad_obj.pick_nD(par.iloc[g][self.P['idxs'][4]:self.P['idxs'][5]])
                    elif self.mc['sampling_algorithm'] == 'emcee':
                        agnrad_obj.pick_nD(par[g][self.P['idxs'][4]:self.P['idxs'][5]])
                    all_agnrad_nus, agnrad_Fnus = agnrad_obj.get_fluxes(agnrad_obj.matched_parkeys)
                else:                                                                #If there is a radio model with fix parameters
                    all_agnrad_nus, agnrad_Fnus = agnrad_obj.get_fluxes('-99.9')

                RADinterp = scipy.interpolate.interp1d(all_agnrad_nus, np.log10(agnrad_Fnus), bounds_error=False, fill_value=-100)
                all_agnrad_Fnus = 10**RADinterp(self.all_nus_rest)
                all_agnrad_Fnus[self.all_nus_rest>=17.5]= 0
                RADFnu =   all_agnrad_Fnus * 10**float(RAD)
                RADFnu_list.append(RADFnu)

            SBinterp = scipy.interpolate.interp1d(all_sb_nus, np.log10(sb_Fnus.flatten()), bounds_error=False, fill_value=-100)
            all_sb_Fnus = 10**SBinterp(self.all_nus_rest)

            #If the UV-Xray correlation was applied, interpolate in a separately the emission of accretion disk and Xrays
            if models.settings['BBB'] !='KD18' and models.settings['XRAYS']==True:
                BBinterp1 = scipy.interpolate.interp1d(all_bbb_nus[all_bbb_nus < 16.685], np.log10(bbb_Fnus.flatten()[all_bbb_nus < 16.685]), bounds_error=False, fill_value=-100) ##
                BBinterp2 = scipy.interpolate.interp1d(all_bbb_nus[all_bbb_nus >= 16.685], np.log10(bbb_Fnus.flatten()[all_bbb_nus >= 16.685]),  bounds_error=False, fill_value=-100) 
                all_bbb_Fnus = np.concatenate((10**BBinterp1(self.all_nus_rest[ self.all_nus_rest < 16.685]), 10**BBinterp2(self.all_nus_rest[ self.all_nus_rest >= 16.685])))
            else:
                BBinterp = scipy.interpolate.interp1d(all_bbb_nus, np.log10(bbb_Fnus.flatten()), bounds_error=False, fill_value=-100)
                all_bbb_Fnus = 10**BBinterp(self.all_nus_rest)

            ### Plot dereddened
            if (models.settings['BBB']=='R06' or models.settings['BBB']=='THB21') and models.settings['XRAYS'] != True: 
                all_bbb_nus, bbb_Fnus_deredd = MD.BBBFdict_4plot['0.0']
                BBderedinterp = scipy.interpolate.interp1d(all_bbb_nus, np.log10(bbb_Fnus_deredd.flatten()), bounds_error=False, fill_value=-100)
                all_bbb_Fnus_deredd = 10**BBderedinterp(self.all_nus_rest)

            elif models.settings['BBB']=='SN12' and models.settings['XRAYS'] != True: 
                all_bbb_nus, bbb_Fnus_deredd = MD.BBBFdict[tuple(np.append(bbb_obj.matched_parkeys[:-1], 0.0))] 
                BBderedinterp = scipy.interpolate.interp1d(all_bbb_nus, np.log10(bbb_Fnus_deredd.flatten()), bounds_error=False, fill_value=-100)
                all_bbb_Fnus_deredd = 10**BBderedinterp(self.all_nus_rest)

            elif models.settings['XRAYS'] == True: 
                EBVbbb_pos = bbb_obj.par_names.index('EBVbbb')
                params = bbb_obj.matched_parkeys_grid
                params[EBVbbb_pos] = str(0.0)
                all_bbb_nus, bbb_Fnus_deredd = MD.BBBFdict_4plot[tuple(params)]                    #Intrinsic fluxes without reddening
                BBinterp1 = scipy.interpolate.interp1d(all_bbb_nus[all_bbb_nus < 16.685], np.log10(bbb_Fnus_deredd.flatten()[all_bbb_nus < 16.685]), bounds_error=False, fill_value=-100) ##
                BBinterp2 = scipy.interpolate.interp1d(all_bbb_nus[all_bbb_nus >= 16.685], np.log10(bbb_Fnus_deredd.flatten()[all_bbb_nus >= 16.685]),  bounds_error=False, fill_value=-100) 
                all_bbb_Fnus_deredd = np.concatenate((10**BBinterp1(self.all_nus_rest[ self.all_nus_rest < 16.685]), 10**BBinterp2(self.all_nus_rest[ self.all_nus_rest >= 16.685])))
            else:
                all_bbb_Fnus_deredd = all_bbb_Fnus

            TOinterp = scipy.interpolate.interp1d(all_tor_nus, np.log10(tor_Fnus.flatten()),  kind = 'linear', bounds_error=False, fill_value=-100) #quadratic
            all_tor_Fnus = 10**TOinterp(self.all_nus_rest)   
            all_tor_Fnus[self.all_nus_rest>16]= 0
            all_tor_Fnus[self.all_nus_rest<11.7]= 0

            if self.output_type == 'plot':
                if self.mc['sampling_algorithm'] == 'ultranest':
                    par2= par.iloc[g]
                elif self.mc['sampling_algorithm'] == 'emcee':
                    par2= par[g]
                filtered_modelpoints, _ = parspace.ymodel(data.nus,data.z, data.dlum, models, self.P, *par2)
            #Using the costumized normalization 
            SBFnu =   all_sb_Fnus *10**float(SB) 
            if (models.settings['BBB']=='R06' or models.settings['BBB']=='THB21'): 
                BBFnu = all_bbb_Fnus * 10**float(BB) 
                BBFnu_deredd = all_bbb_Fnus_deredd * 10**float(BB)
            else:
                BBFnu = (all_bbb_Fnus /(4*math.pi*data.dlum**2)) * 10**float(BB) 
                BBFnu_deredd = (all_bbb_Fnus_deredd /(4*math.pi*data.dlum**2)) * 10**float(BB) 


            GAFnu =   all_gal_Fnus * 10**float(GA) 
            TOFnu =   all_tor_Fnus * 10**float(TO)

            TOTALFnu =  SBFnu + BBFnu + GAFnu + TOFnu

            if models.settings['RADIO'] == True:
                TOTALFnu +=  RADFnu
            
            #Append to the list for all realizations
            SBFnu_list.append(SBFnu)
            BBFnu_list.append(BBFnu)
            GAFnu_list.append(GAFnu)
            TOFnu_list.append(TOFnu)
            TOTALFnu_list.append(TOTALFnu)
            BBFnu_deredd_list.append(BBFnu_deredd)
            #Only if SED plotting: do the same with the  modelled flux values at each data point 
            if self.output_type == 'plot':
                filtered_modelpoints_list.append(filtered_modelpoints)
              
        #Convert lists into Numpy arrays
        SBFnu_array = np.array(SBFnu_list)
        BBFnu_array = np.array(BBFnu_list)
        GAFnu_array = np.array(GAFnu_list)
        TOFnu_array = np.array(TOFnu_list)
        TOTALFnu_array = np.array(TOTALFnu_list)
        BBFnu_array_deredd = np.array(BBFnu_deredd_list)   

        #Put them all together to transport
        if models.settings['RADIO'] == True:
            RADFnu_array = np.array(RADFnu_list)
            FLUXES4plotting = (SBFnu_array, BBFnu_array, GAFnu_array, TOFnu_array, RADFnu_array, TOTALFnu_array,BBFnu_array_deredd)
        elif models.settings['RADIO'] == False:
            FLUXES4plotting = (SBFnu_array, BBFnu_array, GAFnu_array, TOFnu_array, TOTALFnu_array,BBFnu_array_deredd)

        #Convert Fluxes to nuLnu
        self.nuLnus4plotting = self.FLUXES2nuLnu_4plotting(self.all_nus_rest, FLUXES4plotting, data.z)

        #Only if SED plotting:
        if self.output_type == 'plot':
            filtered_modelpoints = np.array(filtered_modelpoints_list)
            distance= model.z2Dlum(data.z)
            lumfactor = (4. * math.pi * distance**2.)
            self.filtered_modelpoints_nuLnu = (filtered_modelpoints *lumfactor* 10**(data.nus))
        #Only if calculating integrated luminosities:    
        if self.output_type == 'int_lums':
            #Convert Fluxes to nuLnu
            self.int_lums= np.log10(self.integrated_luminosities(self.out ,self.all_nus_rest, self.nuLnus4plotting))
            #self.int_lums = int_lums[:,:-1]  # all except last one which is best fit
            self.int_lums_best = self.int_lums[:,-1]  # last one, not yet working

    def FLUXES2nuLnu_4plotting(self, all_nus_rest, FLUXES4plotting, z):

        """
        Converts FLUXES4plotting into nuLnu_4plotting.

        ##input: 
        - all_nus_rest (give in 10^lognu, not log.)
        - FLUXES4plotting : fluxes for the four models corresponding
                            to each element of the total chain
        - source redshift z                    
        """

        all_nus_obs = 10**all_nus_rest /(1+z) 
        distance= model.z2Dlum(z)
        lumfactor = (4. * math.pi * distance**2.)
        if len(FLUXES4plotting) == 7:       #If the AGN radio component was included
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, RADnuLnu, TOTALnuLnu, BBnuLnu_deredd = [ f *lumfactor*all_nus_obs for f in FLUXES4plotting]
            return SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, RADnuLnu, TOTALnuLnu, BBnuLnu_deredd
        else:
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, TOTALnuLnu, BBnuLnu_deredd = [ f *lumfactor*all_nus_obs for f in FLUXES4plotting]
            return SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, TOTALnuLnu, BBnuLnu_deredd


    def integrated_luminosities(self,out ,all_nus_rest, nuLnus4plotting):

        """
        Calculates the integrated luminosities for 
        all model templates chosen by the user in out['intlum_models'], 
        within the integration ranges given by out['intlum_freqranges'].

        ##input: 
        - settings dictionary out[]
        - all_nus_rest
        - nuLnus4plotting: nu*luminosities for the four models corresponding
                            to each element of the total chain
        """
        if len(nuLnus4plotting) == 7:      #If the AGN radio component was included
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, RADnuLnu, TOTALnuLnu, BBnuLnu_deredd =nuLnus4plotting
        elif len(nuLnus4plotting) == 6:
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, TOTALnuLnu, BBnuLnu_deredd =nuLnus4plotting

        out['intlum_freqranges'] = (out['intlum_freqranges']*out['intlum_freqranges_unit']).to(u.Hz, equivalencies=u.spectral())
        int_lums = []
        for m in range(len(out['intlum_models'])):

            if out['intlum_models'][m] == 'sb':    
                nuLnu= SBnuLnu
            elif (out['intlum_models'][m] == 'bbb') or (out['intlum_models'][m] == 'Lx0.5-2keV') or (out['intlum_models'][m] == 'Lx2-8keV'):    
                nuLnu= BBnuLnu
            elif out['intlum_models'][m] == 'bbbdered':    
                nuLnu=BBnuLnu_deredd
            elif out['intlum_models'][m] == 'gal':    
                nuLnu=GAnuLnu
            elif out['intlum_models'][m] == 'tor':    
                nuLnu=TOnuLnu
            elif out['intlum_models'][m] == 'agn_rad':    
                nuLnu=RADnuLnu
            elif (out['intlum_models'][m] == 'AGNfrac_IR') or (out['intlum_models'][m] =='tor+bbb') or (out['intlum_models'][m] == 'AGNfracTO') or (out['intlum_models'][m] == 'AGNfracGA') or (out['intlum_models'][m] == 'AGNfracSB') or (out['intlum_models'][m] == 'AGNfrac_opt') or (out['intlum_models'][m] =='gal+bbb') or (out['intlum_models'][m] == 'AGNfrac_rad') or (out['intlum_models'][m].find('AGNfrac') == 0):    
                nuLnuto=TOnuLnu
                nuLnusb=SBnuLnu
                nuLnuga=GAnuLnu
                nuLnubb=BBnuLnu #_deredd  #BBnuLnu

                if (out['intlum_models'][m] == 'AGNfrac_rad') or (out['intlum_models'][m] == 'AGNfrac1-10GHz') or (out['intlum_models'][m] == 'AGNfrac10-30GHz') or (out['intlum_models'][m] == 'AGNfrac50-200GHz'):
                    nuLnurad=RADnuLnu  #In the cases in which the AGNfrac accounts for radio fluxes
                
                # In these cases, the AGNfraction and the AGN luminosity values are calculated and saved at the same time
                if out['intlum_freqranges'][m][0] ==out['intlum_freqranges'][m][1]: ### AGNfraction or AGN luminosity for monochromatic luminosities
                    index  = (np.abs(all_nus_rest - np.log10(out['intlum_freqranges'][m][0].value))).argmin()
                    Lnuto = nuLnuto[:,index].ravel() 
                    Lnusb = nuLnusb[:,index].ravel()
                    Lnuga = nuLnuga[:,index].ravel()
                    Lnubb = nuLnubb[:,index].ravel()

                    if out['intlum_models'][m] == 'AGNfrac':
                        AGNfrac = (Lnuto)/(Lnuto+Lnusb+Lnuga)
                        int_lums.append(10**AGNfrac)
                    elif out['intlum_models'][m] == 'tor+bbb':
                        int_lums.append(Lnuto+Lnubb)
                    elif out['intlum_models'][m] == 'AGNfrac_rad':
                        Lnurad = nuLnurad[:,index].ravel()
                        AGNfrac = Lnurad/(Lnuto+Lnusb+Lnurad)
                        int_lums.append(10**AGNfrac)

                else:                                                       ### AGNfraction or AGN luminosity for luminosities in a band
                    index  = ((all_nus_rest >= np.log10(out['intlum_freqranges'][m][1].value)) & (all_nus_rest<= np.log10(out['intlum_freqranges'][m][0].value)))            
                    all_nus_rest_int = 10**(all_nus_rest[index])
                    Lnuto = nuLnuto[:,index] / all_nus_rest_int
                    Lnusb = nuLnusb[:,index] / all_nus_rest_int
                    Lnuga = nuLnuga[:,index] / all_nus_rest_int
                    if (out['intlum_models'][m] == 'AGNfrac_rad') or (out['intlum_models'][m] == 'AGNfrac1-10GHz') or (out['intlum_models'][m] == 'AGNfrac10-30GHz') or (out['intlum_models'][m] == 'AGNfrac50-200GHz'):
                        Lnurad = nuLnurad[:,index] / all_nus_rest_int
                        Lnurad_int = scipy.integrate.trapz(Lnurad, x=all_nus_rest_int)
                    Lnubb = nuLnubb[:,index] / all_nus_rest_int
                    Lnuto_int = scipy.integrate.trapz(Lnuto, x=all_nus_rest_int)
                    Lnusb_int = scipy.integrate.trapz(Lnusb, x=all_nus_rest_int)
                    Lnuga_int = scipy.integrate.trapz(Lnuga, x=all_nus_rest_int)
                    Lnubb_int = scipy.integrate.trapz(Lnubb, x=all_nus_rest_int)

                    if out['intlum_models'][m] == 'AGNfrac_IR':                
                        AGNfrac = (Lnuto_int)/(Lnuto_int+Lnusb_int) 
                        int_lums.append(10**AGNfrac)
                    elif out['intlum_models'][m] == 'AGNfracTO':                
                        AGNfrac = (Lnuto_int)/(Lnuto_int+Lnusb_int+Lnuga_int)
                        int_lums.append(10**AGNfrac)
                    elif out['intlum_models'][m] == 'AGNfracGA':                
                        AGNfrac = (Lnuga_int)/(Lnuto_int+Lnusb_int+Lnuga_int)
                        int_lums.append(10**AGNfrac)
                    elif out['intlum_models'][m] == 'AGNfracSB':                
                        AGNfrac = (Lnusb_int)/(Lnuto_int+Lnusb_int+Lnuga_int)
                        int_lums.append(10**AGNfrac)
                    elif out['intlum_models'][m] == 'AGNfrac_opt':                

                        #AGNfrac = Lnuga_int/Lnubb_int
                        AGNfrac = Lnubb_int/(Lnuga_int+Lnubb_int)

                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac_rad'):
                        AGNfrac = Lnurad_int/(Lnuto_int+Lnusb_int+Lnurad_int)
                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac'):
                        AGNfrac = Lnurad_int/(Lnuto_int+Lnusb_int+Lnurad_int)
                        int_lums.append(10**AGNfrac)
                    elif out['intlum_models'][m] == 'gal+bbb':
                        int_lums.append(Lnuga_int+Lnubb_int)
                    elif out['intlum_models'][m] == 'tor+bbb':
                        int_lums.append(Lnuto_int+Lnubb_int)
                    elif (out['intlum_models'][m] == 'AGNfrac1-10GHz') or (out['intlum_models'][m] == 'AGNfrac10-30GHz') or (out['intlum_models'][m] == 'AGNfrac50-200GHz') :
                        AGNfrac = Lnurad_int/(Lnusb_int+Lnurad_int)
                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac400-1600GHz'):
                        if self.models_settings['RADIO'] == False:
                            Lnurad_int = 0
                        AGNfrac = (Lnurad_int+Lnuto_int)/(Lnusb_int+Lnurad_int+Lnuto_int)
                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac187-75mu'):
                        if self.models_settings['RADIO'] == False:
                            Lnurad_int = 0
                        AGNfrac = (Lnurad_int+Lnuto_int)/(Lnusb_int+Lnurad_int+Lnuto_int+Lnuga_int)
                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac75-20mu') or (out['intlum_models'][m] == 'AGNfrac20-8mu') or (out['intlum_models'][m] == 'AGNfrac8-4mu'):
                        AGNfrac = (Lnuto_int)/(Lnusb_int+Lnuto_int+Lnuga_int+Lnubb_int)
                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac4-2mu') or (out['intlum_models'][m] == 'AGNfrac2-1mu')  or (out['intlum_models'][m] == 'AGNfrac1-0.5mu'):
                        AGNfrac = (Lnuto_int+Lnubb_int)/(Lnusb_int+Lnuto_int+Lnuga_int+Lnubb_int)
                        int_lums.append(10**AGNfrac)
                    elif (out['intlum_models'][m] == 'AGNfrac0.5-0.3mu') or (out['intlum_models'][m] == 'AGNfrac0.3-0.06mu'):
                        AGNfrac = (Lnubb_int)/(Lnuto_int+Lnuga_int+Lnubb_int)
                        int_lums.append(10**AGNfrac)
			
            #Once the nuLnu is defined, the luminosity is calculated and saved                
            if (out['intlum_models'][m] != 'AGNfrac_IR') and (out['intlum_models'][m] !='tor+bbb') and (out['intlum_models'][m] != 'AGNfracTO') and (out['intlum_models'][m] != 'AGNfracGA') and (out['intlum_models'][m] != 'AGNfracSB') and (out['intlum_models'][m] != 'AGNfrac_opt') and  (out['intlum_models'][m] !='gal+bbb') and (out['intlum_models'][m] != 'AGNfrac_rad') and (out['intlum_models'][m].find('AGNfrac') == -1): 
                if out['intlum_freqranges'][m][0] ==out['intlum_freqranges'][m][1]: ### monochromatic luminosities
                    index  = (np.abs(all_nus_rest - np.log10(out['intlum_freqranges'][m][0].value))).argmin()   
                    Lnu_mono = nuLnu[:,index].ravel()
                    int_lums.append(Lnu_mono)
                else:
                    index  = ((all_nus_rest >= np.log10(out['intlum_freqranges'][m][1].value)) & (all_nus_rest<= np.log10(out['intlum_freqranges'][m][0].value)))            
                    all_nus_rest_int = 10**(all_nus_rest[index])
                    Lnu = nuLnu[:,index] / all_nus_rest_int
                    Lnu_int = scipy.integrate.trapz(Lnu, x=all_nus_rest_int)
                    int_lums.append(Lnu_int)

        return np.array(int_lums)



"""
Some stand-alone functions on the SED plot format
"""

def SED_plotting_settings(x, ydata, modeldata, models, out, plot_residuals= False):

    """
    This function produces the setting for the figures for SED plotting.
    **Input:
    - all nus, and data (to make the plot limits depending on the data)
    """
    #-- Latex -------------------------------------------------
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('axes', linewidth=2)
    #-------------------------------------------------------------

    output_plots = []
    fig = plt.figure(figsize=(16,6)) 


    ax1 = fig.add_axes([0.15,0.3,0.8,0.6])  
    ax2 = ax1.twiny()
    
    ax1.set_xlabel(r'rest-frame $\mathbf{log \  \nu} [\mathtt{Hz}] $', fontsize=18)
    ax2.set_xlabel(r'$\mathbf{\lambda} [\mathtt{\mu m}] $', fontsize=18)
    ax1.set_ylabel(r'$\mathbf{\nu L(\nu) [\mathtt{erg \ } \mathtt{ s}^{-1}]}$',fontsize=18)

    ax1.tick_params(axis='both',reset=False,which='major',length=8,width=1.5, labelsize = 18, direction='in',top='on',right='on')
    ax1.tick_params(axis='both',reset=False,which='minor',length=4,width=1.5, labelsize = 18, direction='in',top='on',right='on')
    locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ),  numticks=8) 
    ax1.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs= np.arange(2,10)*0.1, numticks=15)
    ax2.xaxis.set_minor_locator(locmin)
    ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    ax1.set_autoscalex_on(True) 
    ax1.set_autoscaley_on(True) 
    ax1.set_xscale('linear')
    ax1.set_yscale('log')
    mediandata = np.median(ydata)

    ax1.set_ylim(min(ydata)*4*1e-2, max(ydata)*8) 
    ax1.set_xlim(min(np.log10(x)), max(np.log10(x))) 

    if out['band_indicators'] == True:
        max_ann = 10**((np.log10(max(ydata)*8)-np.log10(min(ydata)*4*1e-1))*0.92 + np.log10(min(ydata)*4*1e-1))

        if models.settings['RADIO'] == True:
            ax1.annotate('RADIO', (10, max_ann), color= 'black', rotation = 0, fontsize = 11)

        ax1.annotate('FIR', (12, max_ann), color= 'black', rotation = 0, fontsize = 11)
        ax1.annotate('MIR', (13, max_ann), color= 'black', rotation = 0, fontsize = 11)
        ax1.annotate('NIR', (14, max_ann), color= 'black', rotation = 0, fontsize = 11)
        ax1.annotate('OPT', (14.7, max_ann), color= 'black', rotation = 0, fontsize = 11)
        ax1.annotate('UV', (15.5, max_ann), color= 'black', rotation = 0, fontsize = 11)
        ax1.annotate('X-RAY', (18, max_ann), color= 'black', rotation = 0, fontsize = 11)

    ax2.set_xscale('log')
    ax2.set_yscale('log')


    #ax2.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax2.tick_params(axis='both',reset=False,which='major',length=8,width=1.5, labelsize = 18, direction='in',top='on',right='on')
    ax2.tick_params(axis='both',reset=False,which='minor',length=4,width=1.5, labelsize = 18, direction='in',top='on',right='on')

    x2 = (2.98e14/ x)[::-1] # Wavelenght axis

    ax2.plot(x2, np.ones(len(x2)), alpha=0)
    ax2.invert_xaxis()
    locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ),  numticks=15) 
    ax2.xaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs= np.arange(2,10)*0.1, numticks=15)
    ax2.xaxis.set_minor_locator(locmin)
    ax2.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    output_plots.append(fig)
    output_plots.append(ax1)
    output_plots.append(ax2)


    if plot_residuals==True:
        axr = fig.add_axes([0.15,0.1,0.8,0.2],sharex=ax1)
        axr.set_xlabel(r'rest-frame ${\log \  \nu}$ $[\mathrm{Hz}] $', fontsize=18)
        axr.set_ylabel(r'residual $[\sigma]$', fontsize=18)
        axr.set_autoscalex_on(True) 
        axr.set_xscale('linear')
        axr.minorticks_on()
        axr.tick_params(axis='x',reset=False,which='major',length=8,width=1.5, labelsize = 18, direction='in',right='on')
        axr.tick_params(axis='x',reset=False,which='minor',length=4,width=1.5, labelsize = 18, direction='in',right='on')
        axr.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
        axr.xaxis.set_major_locator(ticker.MultipleLocator(1))
        axr.tick_params(labelsize = 19)
        xr = np.log10(x[::-1]) # frequency axis
        axr.plot(xr, np.zeros(len(xr)), 'gray', alpha=1)
        axr.set_xlim(min(np.log10(x)), max(np.log10(x))) 
        ax1.xaxis.set_visible(False)
        output_plots.append(axr)

    return output_plots


def SED_colors(combination = 'a'):

    if combination=='a':   
        steelblue = '#4682b4'
        darkcyan ='#009acd'
        deepbluesky = '#008b8b'
        seagreen = '#2E8B57'    
        lila = '#68228B'
        darkblue='#123281'

    return seagreen, darkblue, 'orange', lila, darkcyan, 'red'


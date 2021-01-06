

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
#matplotlib.use('Agg')
# import matplotlib
# matplotlib.use('TkAgg')
import sys, os
import math 
import numpy as np
from . import corner #Author: Dan Foreman-Mackey (danfm@nyu.edu)
import scipy
from astropy import units as u
from astropy import constants as const

#AGNfitter IMPORTS
from . import MODEL_AGNfitter as model
from . import PARAMETERSPACE_AGNfitter as parspace
import pickle



###!!!def main(data, P, out, models_settings):
def main(data, models, P, out, models_settings):


    """
    Main function of PLOTandWRITE_AGNfitter. Output depends of settings in RUN_AGNfitter.

    ##input:
    - data object
    - dictionary P (parameter space settings)
    - dictionary out (output settings)
    
    """

    chain_burnin = CHAIN(data.output_folder+str(data.name)+ '/samples_burn1-2-3.sav', out)
    chain_mcmc = CHAIN(data.output_folder+str(data.name)+ '/samples_mcmc.sav',  out)
    chain_mcmc.props()

    print( '_________________________________')
    print( 'Properties of the sampling results:')
    print( '- Mean acceptance fraction', chain_mcmc.mean_accept)
    print( '- Mean autocorrelation time', chain_mcmc.mean_autocorr)

    output = OUTPUT(chain_mcmc, data, models, P, models_settings)

    if out['plot_tracesburn-in']:
        fig, nplot=chain_burnin.plot_trace(P)
        fig.suptitle('Chain traces for %i of %i walkers.' % (nplot,chain_burnin.nwalkers))
        fig.savefig(data.output_folder+str(data.name)+'/traces_burnin.' + out['plot_format'])
        plt.close(fig)

    if out['plot_tracesmcmc']:
        fig, nplot = chain_mcmc.plot_trace(P)
        fig.suptitle('Chain traces for %i of %i walkers.'% (nplot,chain_mcmc.nwalkers))
        fig.savefig(data.output_folder+str(data.name)+'/traces_mcmc.' + out['plot_format'])
        plt.close(fig)

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
            SEDs_header = '#freq '+' '.join(['SBnuLnu'+str(i) for i in range(n_rp)]) +' ' +' '.join(['BBnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['GAnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOTALnuLnu'+str(i) for i in range(n_rp)]) +' '+' '.join(['BBnuLnu_deredd'+str(i) for i in range(n_rp)]) 
        else:
            SEDs_header = '#freq '+' '.join(['SBnuLnu'+str(i) for i in range(n_rp)]) +' ' +' '.join(['BBnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['GAnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['RADnuLnu'+str(i) for i in range(n_rp)])+' '+' '.join(['TOTALnuLnu'+str(i) for i in range(n_rp)]) +' '+' '.join(['BBnuLnu_deredd'+str(i) for i in range(n_rp)]) 

        if out['saveSEDrealizations']:
            np.savetxt(data.output_folder + str(data.name)+'/output_SEDs100_'+str(data.name)+'.txt' , save_SEDS, delimiter = " ",fmt= "%1.4f" ,header= SEDs_header, comments='')
        if out['saveSEDresiduals']:
            #res_header = '#freq '+' '.join(['res'+str(i) for i in range(n_rp)]) 
            res_header = '#freq res'
            np.savetxt(data.output_folder + str(data.name)+'/output_residuals100_'+str(data.name)+'.txt' , save_residuals, delimiter = " ",fmt= "%1.4f" ,header= res_header, comments='')
        if out['plotSEDrealizations']:
            fig.savefig(data.output_folder+str(data.name)+'/SED_manyrealizations_' +str(data.name)+ '.'+out['plot_format'])
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

    def __init__(self, chain_obj, data_obj, model_obj, P, models_settings):

        self.chain = chain_obj
        self.chain.props()

        self.out = chain_obj.out
        self.data = data_obj
        self.models = model_obj
        self.models_settings = models_settings
        self.z=self.data.z
        fluxobj_withintlums = FLUXES_ARRAYS(chain_obj, P,  self.out,'int_lums', self.models_settings)
        fluxobj_4SEDplots = FLUXES_ARRAYS(chain_obj, P, self.out,'plot', self.models_settings)

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

        #Mstar0, SFR_opt0 = model.stellar_info_array(np.array([self.chain.best_fit_pars]), self.data, 1)
        Mstar, SFR_opt = model.stellar_info_array(self.chain.flatchain, self.data, self.models, self.out['realizations2int'])
        column_names = np.transpose(np.array(["P025","P16","P50","P84","P975"], dtype='|S3'))
        #chain_pars_best = np.hstack((self.chain.best_fit_pars, Mstar0, SFR_opt0))
        chain_pars = np.column_stack((self.chain.flatchain_sorted, Mstar, SFR_opt))              

        SFR_IR = model.sfr_IR(self.int_lums[0]) #check that ['intlum_names'][0] is always L_IR(8-100)        
        #SFR_IR_best = model.sfr_IR(np.array([self.int_lums_best[0]])) #check that ['intlum_names'][0] is always L_IR(8-100)         

        chain_others =np.column_stack((self.int_lums.T, SFR_IR))
        #chain_others_best =np.hstack((self.int_lums_best.T, SFR_IR_best))            
        outputvalues = np.column_stack((np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_pars, [2.5,16, 50, 84,97.5], axis=0))))),
                                        np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_others, [2.5,16, 50, 84,97.5], axis=0))))),
                                        np.transpose(np.percentile(self.chain.lnprob_flat, [2.5,16, 50, 84,97.5], axis=0)) ))  

        #outputvalues_best = np.hstack( (chain_pars_best, chain_others_best, np.max(self.chain.lnprob_flat)) )
        #outputvalues = np.vstack((outputvalues, outputvalues_best))
        outputvalues_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', '-ln_like'))] )

        return outputvalues, outputvalues_header

    def save_posteriors(self, P):     

        Mstar, SFR_opt = model.stellar_info_array(self.chain.flatchain, self.data, self.models, self.out['realizations2int'])
        column_names = np.transpose(np.array(["P025","P16","P50","P84","P975"], dtype='|S3'))
        chain_pars = np.column_stack((self.chain.flatchain, Mstar, SFR_opt))        
        
        if self.out['calc_intlum']:            

            SFR_IR = model.sfr_IR(self.int_lums[0]) #check that ['intlum_names'][0] is always L_IR(8-100)        

            chain_others =np.column_stack((self.int_lums.T, SFR_IR))
            outputvalues = np.column_stack((np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_pars, [2.5,16, 50, 84,97.5], axis=0))))),
                                            np.transpose(list(map(lambda v: (v[0],v[1],v[2],v[3],v[4]), zip(*np.percentile(chain_others, [2.5,16, 50, 84,97.5], axis=0))))),
                                            np.transpose(np.percentile(self.chain.lnprob_flat, [2.5,16, 50, 84,97.5], axis=0)) ))  
   
            outputvalues_header= ' '.join([ i for i in np.hstack((P['names'], 'logMstar', 'SFR_opt', self.out['intlum_names'], 'SFR_IR', '-ln_like'))] )

            if self.out['save_posteriors']:  
                nsample, npar = self.chain.flatchain.shape
                posteriors = np.column_stack((chain_pars[np.random.choice(nsample, (self.out['realizations2int'])),:], chain_others, self.chain.lnprob_flat[np.random.choice(nsample, (self.out['realizations2int']))])) 
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
            fig, ax1, ax2, axr = SED_plotting_settings(all_nus_rest, data_nuLnu_rest, self.allnus, True)
        else:
            fig, ax1, ax2 = SED_plotting_settings(all_nus_rest, data_nuLnu_rest, self.allnus)
        SBcolor, BBcolor, GAcolor, TOcolor, RADcolor, TOTALcolor= SED_colors(combination = 'a')
        lw= 1.5


        alp = 0.25
        mec='None'
        if Nrealizations == 1:
            alp = 1.0
        for i in range(Nrealizations):
            
            # last one is the max likelihood fit
            if i == Nrealizations -1:
                alp = 1
                lw = 2
                mec='k'
                save_residuals= np.column_stack((data_nus, np.array(data_nuLnu_rest-self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.])/data_errors_rest))

            #Settings for model lines
            p2=ax1.plot(all_nus, SBnuLnu[i], marker="None", linewidth=lw, label="1 /sigma", color= SBcolor, alpha = alp)
            p3=ax1.plot(all_nus, BBnuLnu[i], marker="None", linewidth=lw, label="1 /sigma",color= BBcolor, alpha = alp)
            p4=ax1.plot(all_nus, GAnuLnu[i],marker="None", linewidth=lw, label="1 /sigma",color=GAcolor, alpha = alp)
            p5=ax1.plot( all_nus, TOnuLnu[i], marker="None",  linewidth=lw, label="1 /sigma",color= TOcolor ,alpha = alp)
            p1= ax1.plot( all_nus, TOTALnuLnu[i], marker="None", linewidth=lw,  label="1 /sigma", color= TOTALcolor, alpha= alp)

            if self.models_settings['RADIO'] == True:
                p6=ax1.plot( all_nus, RADnuLnu[i], marker="None",  linewidth=lw, label="1 /sigma",color= RADcolor ,alpha = alp)

            det = [yndflags==1]
            upp = [yndflags==0]
            
            p7 = ax1.plot(data_nus, self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.],   marker='o', linestyle="None",markersize=5, color="red", alpha =alp)
            
            if plot_residuals:
                p6r = axr.plot(data_nus[tuple(det)], (data_nuLnu_rest[tuple(det)]-self.filtered_modelpoints_nuLnu[i][self.data.fluxes>0.][tuple(det)])/data_errors_rest[tuple(det)],   marker='o', mec=mec, linestyle="None",markersize=5, color="red", alpha =alp)	
            
            upplimits = ax1.errorbar(data_nus[tuple(upp)], 2.*data_nuLnu_rest[tuple(upp)], yerr= data_errors_rest[tuple(upp)]/2, uplims = True, linestyle='',  markersize=5, color="black")
            (_, caps, _) = ax1.errorbar(data_nus[tuple(det)], data_nuLnu_rest[tuple(det)], yerr= data_errors_rest[tuple(det)], capsize=4, linestyle="None", linewidth=1.5,  marker='o',markersize=5, color="black", alpha = 1)


        ax1.text(0.04, 0.92, r'id='+str(self.data.name)+r', z ='+ str(self.z), ha='left', transform=ax1.transAxes )
        if plot_residuals:
            ax1.text(0.96, 0.92, 'max ln-likelihood = {ml:.1f}'.format(ml=np.max(self.chain.lnprob_flat)), ha='right', transform=ax1.transAxes )
        #ax1.annotate(r'XID='+str(self.data.name)+r', z ='+ str(self.z)+'max log-likelihood = {ml:.1f}'.format(ml=np.max(self.chain.lnprob_flat)), xy=(0, 1),  xycoords='axes points', xytext=(20, 310), textcoords='axes points' )#+ ', log $\mathbf{L}_{\mathbf{IR}}$= ' + str(Lir_agn) +', log $\mathbf{L}_{\mathbf{FIR}}$= ' + str(Lfir) + ',  log $\mathbf{L}_{\mathbf{UV}} $= '+ str(Lbol_agn)
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

    def __init__(self, outputfilename, out):
            self.outputfilename = outputfilename
            self.out = out

    def props(self):
        if os.path.lexists(self.outputfilename):
            f = open(self.outputfilename, 'rb')
            samples = pickle.load(f, encoding='latin1')
            f.close()

            self.chain = samples['chain']
            nwalkers, nsamples, npar = samples['chain'].shape

            Ns, Nt = self.out['Nsample'], self.out['Nthinning']        
            self.lnprob = samples['lnprob']
            self.lnprob_flat = samples['lnprob'][:,0:Ns*Nt:Nt].ravel()

            isort = (- self.lnprob_flat).argsort() #sort parameter vector for likelihood
            lnprob_sorted = np.reshape(self.lnprob_flat[isort],(-1,1))
            self.lnprob_max = lnprob_sorted[0]


            self.flatchain = samples['chain'][:,0:Ns*Nt:Nt,:].reshape(-1, npar)
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
        rc('axes', linewidth=1.5)
        #-------------------------------------------------------------
        self.props()

        self.nwalkers, nsample, npar = self.chain.shape
        nrows = npar + 1
        ncols =1     

        def fig_axes(nrows, ncols, npar, width=13):
            fig = plt.figure(figsize=(width, width*1.6))#*nrows/ncols))    
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


    def __init__(self, chain_obj, P, out, output_type, models_settings):
        self.chain_obj = chain_obj
        self.output_type = output_type
        self.out = out
        self.models_settings=models_settings
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
        RADFnu_list = []
        TOTALFnu_list = []
        BBFnu_deredd_list = []
        if self.output_type == 'plot':
            filtered_modelpoints_list = []

        gal_obj,sb_obj,tor_obj, bbb_obj = models.dictkey_arrays_4plot

        # Take the  4 dictionaries for plotting. Dicts are described in DICTIONARIES_AGNfitter.py
        MD= models.dict_modelfluxes
        nsample, npar = self.chain_obj.flatchain.shape
        source = data.name

        if self.output_type == 'plot':
            par = self.chain_obj.flatchain[np.random.choice(nsample, (self.out['realizations2plot'])),:]#.T        
            par_best = self.chain_obj.best_fit_pars
            par[-1]=par_best

            realization_nr=self.out['realizations2plot']
        
        elif self.output_type == 'int_lums':
            par = self.chain_obj.flatchain[np.random.choice(nsample, (self.out['realizations2int'])),:]#.T
            par_best = self.chain_obj.best_fit_pars
            par[-1]=par_best
            realization_nr=self.out['realizations2int']
        
        elif self.output_type == 'best_fit':
            par = self.chain_obj.best_fit_pars

        if self.models_settings['BBB'] =='D12_S' or self.models_settings['BBB'] =='D12_K' or self.models_settings['XRAYS']==True:
            ### extend SED to X-rays
            lognu_max = 19
        else:
            lognu_max = 16.2

        if self.models_settings['RADIO']==True:
            ## extend SED to radio
            lognu_min = 9
        else:
            lognu_min = 11.5

        self.all_nus_rest = np.arange(lognu_min, lognu_max, 0.001) 
        
        for g in range(realization_nr):
            ## Pick dictionary key-values, nearest to the MCMC- parameter values
            gal_obj.pick_nD(par[g][self.P['idxs'][0]:self.P['idxs'][1]]) 

            tor_obj.pick_nD(par[g][self.P['idxs'][2]:self.P['idxs'][3]])     
            all_tor_nus, tor_Fnus= tor_obj.get_fluxes(tor_obj.matched_parkeys)

            sb_obj.pick_nD(par[g][self.P['idxs'][1]:self.P['idxs'][2]])               
            all_sb_nus, sb_Fnus= sb_obj.get_fluxes(sb_obj.matched_parkeys) 

            if models.settings['BBB']=='R06': 
                if models.settings['RADIO'] == True:
                    GA, SB, TO, BB, RAD = par[g][-5:]
                else:
                    GA, SB, TO, BB = par[g][-4:]   
                bbb_obj.pick_nD(par[g][self.P['idxs'][3]:self.P['idxs'][4]])
                all_bbb_nus, bbb_Fnus = bbb_obj.get_fluxes(bbb_obj.matched_parkeys)
            else:
                if models.settings['RADIO'] == True:
                    GA, SB, TO, RAD = par[g][-4:]
                elif models.settings['RADIO'] == False:
                    GA, SB, TO = par[g][-3:]
                BB = 0.
                bbb_obj.pick_nD(par[g][self.P['idxs'][3]:self.P['idxs'][4]])
                all_bbb_nus, bbb_Fnus = bbb_obj.get_fluxes(bbb_obj.matched_parkeys)


            #Produce model fluxes at all_nus_rest for plotting, through interpolation
            #print(gal_obj.get_fluxes(gal_obj.matched_parkeys))
            all_gal_nus, gal_Fnus = gal_obj.get_fluxes(gal_obj.matched_parkeys)
            GAinterp = scipy.interpolate.interp1d(all_gal_nus, gal_Fnus.flatten(), bounds_error=False, fill_value=0.)
            all_gal_Fnus = GAinterp(self.all_nus_rest)

            if models.settings['RADIO'] == True:
                all_agnrad_nus, agnrad_Fnus = MD.AGN_RADFdict_4plot[[i for i in MD.AGN_RADFdict_4plot.keys()][0]] 
                plt.semilogy(all_agnrad_nus, agnrad_Fnus*10**all_agnrad_nus, marker = '.', linestyle = '')
                RADinterp = scipy.interpolate.interp1d(all_agnrad_nus, agnrad_Fnus, bounds_error=False, fill_value=0.)
                all_agnrad_Fnus = RADinterp(self.all_nus_rest)
                all_agnrad_Fnus[self.all_nus_rest>=16]= 0

            SBinterp = scipy.interpolate.interp1d(all_sb_nus, sb_Fnus.flatten(), bounds_error=False, fill_value=0.)
            all_sb_Fnus = SBinterp(self.all_nus_rest)

            if (models.settings['BBB'] =='SN12' or models.settings['BBB'] =='R06') and models.settings['XRAYS']==True:
                BBinterp1 = scipy.interpolate.interp1d(all_bbb_nus[all_bbb_nus < 17], bbb_Fnus.flatten()[all_bbb_nus < 17], bounds_error=False, fill_value=0.) ##
                BBinterp2 = scipy.interpolate.interp1d(all_bbb_nus[all_bbb_nus >= 17], bbb_Fnus.flatten()[all_bbb_nus >= 17],  bounds_error=False, fill_value=0.) 
                all_bbb_Fnus = np.concatenate((BBinterp1(self.all_nus_rest[ self.all_nus_rest < 17]), BBinterp2(self.all_nus_rest[ self.all_nus_rest >= 17])))
            else:
                BBinterp = scipy.interpolate.interp1d(all_bbb_nus, bbb_Fnus.flatten(), bounds_error=False, fill_value=0.)
                all_bbb_Fnus = BBinterp(self.all_nus_rest)

            ### Plot dereddened
            if models.settings['BBB']=='R06' and bbb_obj.par_names == 'EBVbb': # len(bbb_obj.par_names)==1:     #if BB =! 0:
                all_bbb_nus, bbb_Fnus_deredd = MD.BBBFdict_4plot['0.0']
                BBderedinterp = scipy.interpolate.interp1d(all_bbb_nus, bbb_Fnus_deredd.flatten(), bounds_error=False, fill_value=0.)
                all_bbb_Fnus_deredd = BBderedinterp(self.all_nus_rest)
            else:
                all_bbb_Fnus_deredd = all_bbb_Fnus

            TOinterp = scipy.interpolate.interp1d(all_tor_nus, np.log10(tor_Fnus).flatten(), bounds_error=False, fill_value=0.)
            all_tor_Fnus = 10**(TOinterp(self.all_nus_rest))        
            all_tor_Fnus[self.all_nus_rest>16]= 0
            all_tor_Fnus[self.all_nus_rest<11.7]= 0

            if self.output_type == 'plot':
                par2= par[g]
                filtered_modelpoints, _ = parspace.ymodel(data.nus,data.z, data.dlum, models, self.P, *par2)
            #Using the costumized normalization 
            SBFnu =   all_sb_Fnus *10**float(SB) 
            if models.settings['BBB']=='R06': 
                BBFnu = all_bbb_Fnus * 10**float(BB) 
            else:
                BBFnu = (all_bbb_Fnus /(4*math.pi*data.dlum**2)) * 10**float(BB) 

            GAFnu =   all_gal_Fnus * 10**float(GA) 
            TOFnu =   all_tor_Fnus * 10**float(TO)
            BBFnu_deredd = all_bbb_Fnus_deredd * 10**float(BB)

            if models.settings['RADIO'] == True:
                RADFnu =   all_agnrad_Fnus * 10**float(RAD)
                RADFnu_list.append(RADFnu)
                TOTALFnu =  SBFnu + BBFnu + GAFnu + TOFnu + RADFnu

            TOTALFnu =  SBFnu + BBFnu + GAFnu + TOFnu
            
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
        if models.settings['RADIO'] == True:
            RADFnu_array = np.array(RADFnu_list)
        TOTALFnu_array = np.array(TOTALFnu_list)
        BBFnu_array_deredd = np.array(BBFnu_deredd_list)    

        #Put them all together to transport
        if models.settings['RADIO'] == True:
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

        # elif self.output_type == 'best_fit':
        #     self.filtered_modelpoints_nuLnu = self.FLUXES2nuLnu_4plotting(all_nus_rest,  filtered_modelpoints, self.chain_obj.data.z)


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
        if len(FLUXES4plotting) == 7:
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
        if len(nuLnus4plotting) == 7:
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, RADnuLnu, TOTALnuLnu, BBnuLnu_deredd =nuLnus4plotting
        elif len(nuLnus4plotting) == 6:
            SBnuLnu, BBnuLnu, GAnuLnu, TOnuLnu, TOTALnuLnu, BBnuLnu_deredd =nuLnus4plotting

        out['intlum_freqranges'] = (out['intlum_freqranges']*out['intlum_freqranges_unit']).to(u.Hz, equivalencies=u.spectral())
        int_lums = []
        for m in range(len(out['intlum_models'])):

            if out['intlum_models'][m] == 'sb':    
                nuLnu= SBnuLnu
            elif out['intlum_models'][m] == 'bbb':    
                nuLnu= BBnuLnu
            elif out['intlum_models'][m] == 'bbbdered':    
                nuLnu=BBnuLnu_deredd
            elif out['intlum_models'][m] == 'gal':    
                 nuLnu=GAnuLnu
            elif out['intlum_models'][m] == 'tor':    
                 nuLnu=TOnuLnu
            elif out['intlum_models'][m] == 'agn_rad':    
                    nuLnu=RADnuLnu
            elif out['intlum_models'][m] == 'AGNfrac':    
                nuLnuto=TOnuLnu
                nuLnusb=SBnuLnu
                nuLnuga=GAnuLnu

                if out['intlum_freqranges'][m][0] ==out['intlum_freqranges'][m][1]: ### monochromatic luminosities
                    index  = (np.abs(all_nus_rest - np.log10(out['intlum_freqranges'][m][0].value))).argmin()
                    Lnuto = nuLnuto[:,index].ravel() 
                    Lnusb = nuLnusb[:,index].ravel()
                    Lnuga = nuLnuga[:,index].ravel()
                    AGNfrac = Lnuto/(Lnuto+Lnusb+Lnuga)
                else:
                    index  = ((all_nus_rest >= np.log10(out['intlum_freqranges'][m][1].value)) & (all_nus_rest<= np.log10(out['intlum_freqranges'][m][0].value)))            
                    all_nus_rest_int = 10**(all_nus_rest[index])
                    Lnuto = nuLnuto[:,index] / all_nus_rest_int
                    Lnusb = nuLnusb[:,index] / all_nus_rest_int
                    Lnuga = nuLnuga[:,index] / all_nus_rest_int
                    Lnuto_int = scipy.integrate.trapz(Lnuto, x=all_nus_rest_int)
                    Lnusb_int = scipy.integrate.trapz(Lnusb, x=all_nus_rest_int)
                    Lnuga_int = scipy.integrate.trapz(Lnuga, x=all_nus_rest_int)
                    AGNfrac = Lnuto_int/(Lnuto_int+Lnusb_int+Lnuga_int)
                             
            if out['intlum_models'][m] != 'AGNfrac': 

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
            else:
                int_lums.append(AGNfrac)

        return np.array(int_lums)



"""
Some stand-alone functions on the SED plot format
"""

def SED_plotting_settings(x, ydata, modeldata, plot_residuals= False):

    """
    This function produces the setting for the figures for SED plotting.
    **Input:
    - all nus, and data (to make the plot limits depending on the data)
    """
    fig = plt.figure(figsize=(9,5))
    ax1 = fig.add_axes([0.15,0.3,0.8,0.6])
    #ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()

    #-- Latex -------------------------------------------------
    rc('text', usetex=True)
    rc('font', family='serif')
    rc('axes', linewidth=1.5)
    #-------------------------------------------------------------

    #    ax1.set_title(r"\textbf{SED of Type 2}" + r"\textbf{ AGN }"+ "Source Nr. "+ source + "\n . \n . \n ." , fontsize=17, color='k')    
    ax1.set_xlabel(r'rest-frame $\mathbf{log \  \nu} [\mathtt{Hz}] $', fontsize=13)
    ax2.set_xlabel(r'$\mathbf{\lambda} [\mathtt{\mu m}] $', fontsize=13)
    ax1.set_ylabel(r'$\mathbf{\nu L(\nu) [\mathtt{erg \ } \mathtt{ s}^{-1}]}$',fontsize=13)

    ax1.tick_params(axis='both',reset=False,which='major',length=8,width=1.5)
    ax1.tick_params(axis='both',reset=False,which='minor',length=4,width=1.5)

    ax1.set_autoscalex_on(True) 
    ax1.set_autoscaley_on(True) 
    ax1.set_xscale('linear')
    ax1.set_yscale('log')


    mediandata = np.median(ydata)
    ax1.set_ylim(min(ydata)*4*1e-1, max(ydata)*2) #ax1.set_ylim(6*1e37, 1e47) # mediandata /90., mediandata * 50.
    ax1.set_xlim(min(np.log10(x)), max(np.log10(x))) #min(np.log10(x)), max(np.log10(x)) 9,19
    #ax1.set_xlim(min(modeldata), max(modeldata))

    ax2.set_xscale('log')
    ax2.set_yscale('log')


    #ax2.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax2.tick_params(axis='both',reset=False,which='major',length=8,width=1.5)
    ax2.tick_params(axis='both',reset=False,which='minor',length=4,width=1.5)

    x2 = (2.98e14/ x)[::-1] # Wavelenght axis

    ax2.plot(x2, np.ones(len(x2)), alpha=0)
    ax2.invert_xaxis()
    #ax2.set_xticks([1e3, 100., 10.,1., 0.1, 0.01]) 

    if plot_residuals==True:
        axr = fig.add_axes([0.15,0.1,0.8,0.2],sharex=ax1)
        axr.set_xlabel(r'rest-frame ${\log \  \nu}$ $[\mathrm{Hz}] $')
        axr.set_ylabel(r'residual $[\sigma]$')
        axr.set_ylabel(r'residual $[\sigma]$')
        axr.set_autoscalex_on(True) 
        axr.set_xscale('linear')
        axr.minorticks_on()
        xr = np.log10(x[::-1]) # frequency axis
        axr.plot(xr, np.zeros(len(xr)), 'gray', alpha=1)
        axr.set_xlim(min(modeldata), max(modeldata))
        ax1.xaxis.set_visible(False)

        return fig, ax1, ax2, axr
    else:
        return fig, ax1, ax2


def SED_colors(combination = 'a'):

    if combination=='a':   
        steelblue = '#4682b4'
        darkcyan ='#009acd'
        deepbluesky = '#008b8b'
        seagreen = '#2E8B57'    
        lila = '#68228B'
        darkblue='#123281'

    return seagreen, darkblue, 'orange', lila, darkcyan, 'red'


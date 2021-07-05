def dao_normal(wave_len, flux,  final_plt=False, bin_wid=100, fit_order=3):
    
    """ The input for this function needs to be numpy arrays for wavelength (wave_len) and flux.
    
    :param  wave_len: numpy arrays for wavelength (wave_len).    
    :param  flux: numpy arrays for flux.
    :param  bin_wid: The width of the bins is set by bin_wid, which need to be an integer. 
    :param  fit_order: Integer: The order for the polynomial fitting
    """  
    
    
    ####these are empty array in which we will store the median fluxes and wavelengths of each bin 
    
    
    array_size = flux.size  
    
    ##n_bins (integer) is the total number of bins that fit in the array 
    n_bins = array_size // bin_wid 
    
    flux_med = np.zeros(n_bins)
    wave_med = np.zeros(n_bins)
    ###binning and finding the median flux and wavelength each bin. 
    ### iterations and bin size a set by the indices in the inputed arrays. 
    
    
    for i in range(0, n_bins): 
        ###the index upper boundary for a each successive bin. 
        ### The upper and lower boundary for each bin account for the counting initalizing at 0
        n = i * bin_wid  + (bin_wid - 1 ) 
        
        ####the lower bound.
        n_low = i * bin_wid 
    
        ###breaking down each array into segments of the set bin size 
        data_seg = flux[n_low: n + 1] 
        wave_seg = wave_len[n_low: n + 1]
        
        ###extracts the median and stores it defined arrays outside the loop
        med =   np.median(data_seg) 
        flux_med[i] = med
    
        mid_wave = np.median(wave_seg) 
    
        wave_med[i] = mid_wave  
     
    #### polynomial fitting of the median fluxes with respect to the median wavelenths 
    ### the integer at the end of the function idicates the order of the polynomial fit . 
    ### Refer to numpy documentation on these functions as needed. 
    con_fit = np.poly1d(np.polyfit(wave_med, flux_med , fit_order))   
    
    ##the final normalization of the spectrum to the median-derived continuum fit.
    normed_spec = flux / con_fit(wave_len)
    
    if final_plt==True:
        wave_fit = np.linspace(min(wave_len), max(wave_len),2000)
    
        fig, ((ax1), (ax2)) = plt.subplots(2,1, figsize=(10,12))  
        ax1.plot(wave_len, flux, zorder=0)
        ax1.scatter(wave_med, flux_med, color='r', linewidths='5.0', zorder=2)
        ax1.plot(wave_fit, con_fit(wave_fit), color='k', linewidth='2.5')
        ax1.set_ylabel('flux ')
    
        ax2.plot(wave_len, normed_spec)
        ax2.set_xlabel('wavelength (Angstroms)')
        ax2.set_ylabel('normalized flux ')
    
    
     
    ### this function will return the normalized spectrum and the polynomial fits used.     
    if final_plt==True:
        return normed_spec, con_fit, fig, ((ax1), (ax2)) 
    else:
        return normed_spec, con_fit




def load_dao(dao_filename):
    """
    This file is will load in the dao file ( or files like it). It will extrac the
    
    :param  filename: this needs to be inserted as string 
    
    """
    hdu_dao = fits.open(dao_filename) 
    flux_dao = hdu_dao[0].data
    
    header_dao = hdu_dao[0].header
    
    pixel0_dao = header_dao['CRPIX1']
    delta_w_dao = header_dao['CDELT1']
    w0_dao = header_dao['CRVAL1']

    wave_dao = w0_dao + np.arange(0,flux_dao.size, 1)*delta_w_dao 
    
    return wave_dao, flux_dao, header_dao



def dao_wrapper(file_name, name_mods="_norm_test", out_type=None, out_plot_type=None): 
    """file_name_manip: this function assumes that it following the typical format where the only period seprates the filename and type. eg(Path/file.type) 
    it assumes it will be written where it was read"""
    filename_dot_split = file_name.split(".") 
    file_name_slash_split = filename_dot_split[0].split("/")   
    modify_name = file_name_slash_split[-1] + name_mods
    
    recomb = file_name_slash_split[0:-1] 
    recomb.append(modify_name)
    file_out_name = '/'.join(recomb) + out_type
    plot_out_name = '/'.join(recomb) + out_plot_type 
    
    wave, orig_flux, file_header  = load_dao(file_name) 
    
    flux_norm, spec_fit, fig, ((ax1), (ax2))=  dao_normal(wave, orig_flux,  final_plt=True, bin_wid=100, fit_order=3)
    
    
    out_table = Table([wave, flux_norm], names=['lambda', 'flux'])
    ascii.write(out_table, file_out_name ,overwrite=True) 
    fig.savefig(plot_out_name)
    
    
    
    

    
    
    ###calling the file
    
    return None
Tools for adding source variability to GULLS microlensing lightcurves. 

Needs exact version of butterpy given in requirements.txt.

It works on single and double sources, though not extensively tested on double. It adds six columns related to the variability.

```
variability_flux: relative variable flux of source star as observed in different filters

variability_flux_bol: bolometric relative variable flux of source star

true_relative_flux_var: True relative flux of microlensing event with source variability. 

true_relative_flux_var_error: Error for above. Note that for binary sources the way I defined this ceases to make a lot of sense, but this column isn't really needed

measured_relative_flux_var: measured relative flux of microlensing event with source variability. 

measured_relative_flux_var_error: error for above.
```
The measured_relative_flux_var columns are the important ones. For cases where we don't add variability (see below), I'm having it fill these six columns with nan.  

There's also two rows added to the header. One has a list of butterpy parameters used. I'll put a description of those parameters on the repository later.
```
#butterpy variables:0.5951699988308047 1.5404863697028155 0.22774927304115827 27.386934541822455 84.26987126112638 15.54076903240235 63.1829670236636 1.2800184926131124 -0.8057103423299585 3619.95 3131.1607182270486 5000 1409
```
 The other is a flag indicates if the lightcurve has source variability or not. Could be used to skip reading in useless nan columns. 
```
#usevar:True
```

I added a couple of  parameters that can be adjusted as we need.

`var_fraction` - number between 0 and 1 to indicate the fraction of lightcurves that should have any variability. Then for each lightcurve it will generate a random number to decide whether to add variability. Since in most cases the impact of the variability is small anyway I think something like 0.5 will be fine. 

`chi2_threshold` - maximum delta_chi2 that can be introduced by the source variability. Checked at the end of variability computation and if exceeded, variability columns are filled with nan. 

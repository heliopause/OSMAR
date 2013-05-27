readme.txt

 1. run_first on any image set (flipud and copy to working disk)
 2. calibrations (geometric and radiometric)
 3. generate source images
 4. viewing angle mapping
 5. incident angle mapping
 6. process samples
		use all available exposure sets to increase SNR (and reduce saturation)
			this probably needs to be done first
				but possibly after gain correction (think about this)
			see Debevec and Malik (1997)
 7. perform normal estimation to correct for non-planar surface
		basically, find actual incident angle for given measurement
			currently, assumption is that sample is flat, since incident angles
			  were estimated using mirror in sample position
 8. calculate BRDF and display results
		need to properly deal with tilted surfaces (see 7.)
		include ERROR BARS (see Foo thesis, OSE book)
			not sure how to incorporate misfocus (human error)
 9. enable combination of all three wavelengths
 		have to scale the image maps properly (spatially and intensity-wise)
 		should be able to specify a particular visual response
10.	create mathematical models

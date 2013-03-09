readme.txt

 1. run_first on any image set (flipud and copy to working disk)
 2. calibrations (geometric and radiometric)
 3. generate source images
 4. viewing angle mapping
 5. incident angle mapping
		obtain full mapping for projector pixels to incident angle
 6. process samples
		use all available exposure sets to increase SNR (and reduce saturation)
			this probably needs to be done first
				but possibly after gain correction (think about this)
 7. calculate BRDF
 8. display results
		with ERROR BARS (see Foo thesis, OSE book)
			not sure how to incorporate misfocus (human error)
 9. perform normal estimation to correct for non-planar surface
		basically, find actual incident angle for given measurement
			currently, assumption is that sample is flat, since incident angles
			  were estimated using mirror in sample position
10. use all three wavelengths
		could later combine according to a particular visual response (predator!)
# my_thesis
This repository contains the code we used in our work for my Master's thesis and our paper.

This project used some of the codes that I mentioned below:

1. [juliet](https://github.com/nespinoza/juliet): This tool is for modelling the transit signal.

2. [limb-darkening code](https://github.com/nespinoza/limb-darkening): This code is to generate limb-darkening coefficients for each target.

Our primary goal in this study was to calculate the limb-darkening coefficients using the transit light curve and to compare it with the limb-darkening coefficients calculated using theoretical model stellar atmospheres. To calculate the limb-darkening effect from transit light curves, we used the recently released TESS data. For more details see my thesis by clicking [here](https://jayshilpatel.files.wordpress.com/2019/05/dissertation_final.pdf).

As described in thesis, this code will download data products from MAST portal and generate limb darkening coefficients from this data using juliet. Additionally this code will also calculate limb darkening coefficients from ATLAS and PHOENIX model stellar atmospheres using limb-darkening code. This code will then compare them and save result as a errorbar plot in Results/comp_us_and_evidance folder. This would be the main result of this work. This code will also calculate the offset present in LDCs calculated using both ways. Temperature variation in these offsets would be saved in Results/variation_with_temp folder.

To confirm that our fit is indeed good, we also compare several planetary parameters calculated using juliet from TESS data with their existing values from the literature. This comparison will be saved in Results/comp_a_r_p folder.

This code will do all these things autometically. But to do this you have to provide a data file of the targets. This file should contains various planetary as well as stellar parameters. collecting_data.py file can help in doing this.

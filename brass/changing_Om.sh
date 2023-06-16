#! /bin/bash

for Ok in {-20..20..5}; do
	Om=$((31 - $Ok))
	ok=$(bc <<< "scale=2; $Ok/100")
	#phase2
	cp phase3_params_Om031_OL069.c changing_Om_phase2.c
	file=changing_Om_phase2.c
	echo "Ok = $Ok, phase1"
	sed -i -e "/.*Path to data1BAO: \/hom.*/c\\#Path to data1BAO: /home/ssanz/TFG/outputs_santi/rustico/Power_Spectrum_eBOSS_NGC_Om0$Om\_OL069.txt" $file
	sed -i -e "/.*Path to data2BAO: \/hom.*/c\\#Path to data1BAO: /home/ssanz/TFG/outputs_santi/rustico/Power_Spectrum_eBOSS_SGC_Om0$Om\_OL069.txt" $file
	sed -i -e "/.*Maximum number of mc.*/c\\#Maximum number of accepted steps: 1000000" $file
	sed -i -e "/.*sampling step.*/c\\#mcmc sampling step (recommended 1.9): 0.5" $file
	sed -i -e "/.*Path to proposal covariance.*/c\\#Path to proposal covariance: /home/ssanz/TFG/lrg_eboss/output/proposal/mcmc_fid_Om031_OL069.txt" $file
	sed -i -e "/.*Path to smooth Plin:.*/c\\#Path to smooth Plin: /home/ssanz/TFG/outputs_santi/linspace_class/psmlin__linspace_Om031_OL069.txt" $file
	sed -i -e "/.*Path to Olin:.*/c\\#Path to Olin: /home/ssanz/TFG/outputs_santi/linspace_class/Olin__linspace_Om031_OL069.txt" $file
	sed -i -e "/.*Path to output.*/c\\#Path to output: /home/ssanz/TFG/outputs_santi/changing_Om/phase2/" $file
	sed -i -e "/.*Identifier of output:.*/c\\#Identifier of output: phase2_run1_Om0$Om\_OL069" $file
	./file_gcc.out changing_Om_phase2.c
	sed -i -e "/.*Path to proposal covariance.*/c\\#Path to proposal covariance: /home/ssanz/TFG/outputs_santi/changing_Om/phase2/mcmc_phase2_run1_Om0$Om\_OL069.txt" $file
	sed -i -e "/.*Identifier of output:.*/c\\#Identifier of output: phase2_run2_Om0$Om\_OL069" $file
	echo "Ok = $Ok, phase3"
	./file_gcc.out changing_Om_phase2.c
	sed -i -e "/.*Maximum number of mc.*/c\\#Maximum number of accepted steps: 1000000" $file
	sed -i -e "/.*sampling step.*/c\\#mcmc sampling step (recommended 1.9): 1.9" $file
	sed -i -e "/.*Path to proposal covariance.*/c\\#Path to proposal covariance: /home/ssanz/TFG/outputs_santi/changing_Om/phase2/mcmc_phase2_run2_Om0$Om\_OL069.txt" $file
	sed -i -e "/.*Identifier of output:.*/c\\#Identifier of output: phase2_run3_Om0$Om\_OL069" $file
	echo "Ok = $Ok, phase3"
	./file_gcc.out changing_Om_phase2.c
done

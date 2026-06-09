#!/broad/sankaranlab/ajlee/conda_envs/sc_env/bin/python


import os

DIR=os.getcwd()


os.system('mkdir '+DIR+'/input_data')


def make_bed():
	#
	input1=open('/broad/sankaranlab/ajlee/projects/hbf_gwas/finemap_data/hbf-META_finemap_20240712_p5e-8.txt','r') # hg38
	output1=open(DIR+'/input_data/hbf_meta.finemap.hg38.bed','w') 
	all_input1=input1.readlines()
	for line in all_input1[1:]:
		each=line.strip().split('\t')

		chrom='chr'+each[1]
		pt1=each[2]
		pt2=str(int(pt1)+1)

		a1=each[3]
		a2=each[4]

		var_id=chrom+':'+pt1+':'+a1+':'+a2
		pip=each[5]

		new_line=[chrom, pt1, pt2, var_id, pip]
		output1.write('\t'.join(new_line)+'\n')
	
	input1.close()
	output1.close()


def lift_over():
	#
	chain_file='/broad/sankaranlab/ajlee/tools/ucsc/chain_files/hg38ToHg19.over.chain.gz'

	command_line='/broad/sankaranlab/ajlee/tools/ucsc/liftOver '+DIR+'/input_data/hbf_meta.finemap.hg38.bed '+chain_file+' '+DIR+'/input_data/hbf_meta.finemap.hg19.bed '+DIR+'/input_data/hbf_meta.finemap.hg19.unmapped.txt'
	os.system(command_line)


make_bed()
lift_over()




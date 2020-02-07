#!/usr/bin/env python
#module load  Python/2.7.16-gimkl-2018b

import os,vcf,random,operator
#os.system("mkdir sexdet")
#os.system("populations -P output_M2/ -M popmap_allNONEGnowak16-48.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.1  --write-single-snp  -O sexdet") # then filter it without

#########params
!!males = ["WK16-02","WK16-16","WK16-97","WP-23","WP-48","WK16-03","WK16-17","WK16-34","WK16-68","WK16-98","WP-08","WP-41","WK16-04","WK16-54","WK16-86","WP-42","WP-50","WK16-38","WK16-87","WK16-27","WK16-40","WK16-56","WK16-73","WP-03","WP-11","WP-12","WP-20","WP-45","WK16-13","WK16-46","WK16-77","WP-05","WP-38","WP-54","WK16-15","WK16-31","WK16-47","WK16-63","WK16-96","WP-14"]
females = ["WK16-33","WK16-66","WK16-80","WP-07","WP-15","WP-40","WK16-49","WK16-85","WP-16","WP-24","WP-49","WK16-37","WK16-70","WP-01","WP-09","WP-17","WP-34","WK16-06","WK16-23","WK16-55","WK16-72","WP-02","WP-10","WP-18","WP-35","WP-43","WP-51","WK16-11","WK16-90","WP-19","WP-36","WP-44","WP-52","WK16.12","WK16-28","WK16-44","WK16-58","WK16-75","WK16-92","WP-04","WP-37","WP-53","WK16-30","WK16-62","WK16-94","WP-13","WP-21","WP-46","WK16-78","WP-06","WP-22","WP-39","WP-47","WK16-20"]
nsim=0
minpersextobespecific = 5   #number sequenced


#INITIALISE
input_vcf=vcf.Reader(fsock=None, filename="sexdet/populations.snps.vcf", compressed=False, prepend_chr=False, strict_whitespace=False)#open the vcf parser
sample_males= random.sample(males,20)
sample_females= random.sample(females,20)
male_markers =[]
female_markers = []
male_specific = 0
female_specific = 0
def write_output(output_folder = "sexdetsimulations" , nsim=nsim  ,nmale_specific = male_specific,nfemale_specific =female_specific ,males_ID = sample_males,females_ID = sample_females,male_markers = male_markers ,female_markers = female_markers,minpersextobespecific= minpersextobespecific):
	if not os.path.exists(output_folder):os.mkdir(output_folder)
	filename =output_folder+"/sim"+str(nsim)+"_min" +str(minpersextobespecific)+ "tobesexspec_"+str(len(males_ID))+"persex_results.txt"
	output=open(filename,"w")
	output.write("sim:"+str(nsim)+"\n")
	output.write("male specific\t"+str(nmale_specific)+"\nfemale specific\t"+str(nfemale_specific)+"\nmales:\t"+"\t".join(males_ID)+"\nfemales:\t"+"\t".join(females_ID)+"\n")
	output.write("male_markers:\n"+"\n".join(male_markers)+"\n")
	output.write("female_markers:\n"+"\n".join(female_markers)+"\n")
	output.close()
	print("wrote: "+filename)
	os.system("cat "+filename)

while nsim<200:
	input_vcf=vcf.Reader(fsock=None, filename="sexdet/populations.snps.vcf", compressed=False, prepend_chr=False, strict_whitespace=False)#open the vcf parser
	sample_males= random.sample(males,20)
	sample_females= random.sample(females,20)
	nsim+=1
	print("simmulation nr:",nsim,"\n")
	i=0
	male_markers =[]
	female_markers = []
	male_specific = 0
	female_specific = 0
	for record in input_vcf:
		i+=1
		#if i%1000==0: print(i)
		male_set = set()
		female_set  = set()
		number_males = 0
		number_females = 0
		#print(male_set)
		for sample in record.samples:
			if sample.sample in sample_males:
				male_set = male_set |  set(sample["GT"].split("/"))
				male_set = set([x for x in male_set if x !="."])
				if  set(sample["GT"].split("/"))!={"."}: number_males +=1
			elif sample.sample in sample_females:
				#print("here") 
				female_set = female_set | set(sample["GT"].split("/"))
				female_set = set([x for x in female_set if x !="."])
				if  set(sample["GT"].split("/"))!={"."}: number_females +=1
			else:
				pass #not in the sample

		if len(male_set)>0 and len(female_set) ==0 and number_males > minpersextobespecific: 
			male_specific+=1
			male_markers.append(str(record.CHROM)+"_"+str(record.POS))
			print("male specific, number_males",number_males,"number_females",number_females)
		if len(female_set)>0  and len(male_set)==0 and number_females > minpersextobespecific:   
			female_specific+=1
			female_markers.append(str(record.CHROM)+"_"+str(record.POS))
			print("female specific, number_males",number_males,"number_females",number_females)
	#print("sim:"+str(nsim))
	#print("male specific\t"+str(male_specific)+"\nfemale specific\t"+str(female_specific)+"\nmales:\t"+"\t".join(sample_males)+"\nfemales:\t"+"\t".join(sample_females))
	#print("male_markers:\n"+"\n".join(male_markers))
	#print("female_markers:\n"+"\n".join(female_markers))
	write_output(output_folder = "sexdetsimulations" , nsim=nsim  ,nmale_specific = male_specific,nfemale_specific =female_specific ,males_ID = sample_males,females_ID = sample_females,male_markers = male_markers ,female_markers = female_markers,minpersextobespecific= minpersextobespecific)




### compute average of male and female for each sim an ind is included into
filename = "sexdetsimulations/sim200_min5tobesexspec_20persex_results.txt"
#def parse(filename):
class Boot():
	pass

def parse(filename):
	sample = Boot()
	lines = [x.strip() for x in open(filename)]
	sample.sim = lines[0].split(":")[1]
	sample.nmale_specific = lines[1].split("\t")[1]
	sample.nfemale_specific = lines[2].split("\t")[1]
	sample.males = lines[3].split(":")[1].split()
	sample.females = lines[4].split(":")[1].split()
	sample.male_markers = lines[(lines.index("male_markers:")+1):lines.index("female_markers:")]
	sample.female_markers = lines[lines.index("female_markers:")+1:]
	return sample


#find number of specific markers average for each ind 
for male in males:
	n_male_markers_per_sim = []
	n_female_markers_per_sim = []
	for filename in os.listdir("sexdetsimulations"):
		a = parse("sexdetsimulations/"+filename)
		if male in a.males:
			n_male_markers_per_sim.append(int(a.nmale_specific))
			n_female_markers_per_sim.append(int(a.nfemale_specific))
	print ("male",male,"in",len(n_male_markers_per_sim),"sims, average male markers:",sum(n_male_markers_per_sim)/len(n_male_markers_per_sim),"average female markers:",sum(n_female_markers_per_sim)/len(n_female_markers_per_sim))


for female in females:
	n_male_markers_per_sim = []
	n_female_markers_per_sim = []
	for filename in os.listdir("sexdetsimulations"):
		a = parse("sexdetsimulations/"+filename)
		if female in a.females:
			n_male_markers_per_sim.append(int(a.nmale_specific))
			n_female_markers_per_sim.append(int(a.nfemale_specific))
	print ("female",female,"in",len(n_male_markers_per_sim),"sims, average male markers:",sum(n_male_markers_per_sim)/len(n_male_markers_per_sim),"average female markers:",sum(n_female_markers_per_sim)/len(n_female_markers_per_sim))


marker_counts_males = {}
marker_counts_females = {}

for filename in os.listdir("sexdetsimulations"):
		a = parse("sexdetsimulations/"+filename)
		print(a)
		for marker in a.male_markers:
			if marker in marker_counts_males.keys():
				marker_counts_males[marker]+=1
			else:
				marker_counts_males[marker]=1
		for marker in a.female_markers:
			if marker in marker_counts_females.keys():
				marker_counts_females[marker]+=1
			else:
				marker_counts_females[marker]=1

more_common_male = sorted(marker_counts_males.items(), key=operator.itemgetter(1))
more_common_female = sorted(marker_counts_females.items(), key=operator.itemgetter(1))

#more_common_male
('36961_32', 24),
 ('55184_13', 26),
 ('9904_14', 27),
 (' 53604_22', 28),
 ('56150_40', 32),
 ('56042_19', 39),
 ('22460_18', 73)]


#more_common_female

 ('72722_26', 15),
 ('47675_9', 17),
 ('94561_31', 19),
 ('91106_48', 23),
 ('74091_20', 25),
 ('52441_14', 47)]
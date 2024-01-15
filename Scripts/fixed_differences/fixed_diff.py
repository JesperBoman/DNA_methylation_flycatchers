##Pseudocode

#1. unzip the vcfs
#2. Identify monomorphic sites in both species individually
#3. Combine the two datasets of fixed sites within species to obtain fixed differences in a file



#1. unzip the vcfs
import os
os.mkdir("temp")

os.system( "gunzip -dc /crex/proj/sllstore2017033/nobackup/work/alexn/vcf/gt_OC.vcf.gz > temp/gt_OC.vcf")
os.system( "gunzip -dc /crex/proj/sllstore2017033/nobackup/work/alexn/vcf/gt_OP.vcf.gz > temp/gt_OP.vcf")


#2. Identify monomorphic sites in both species Individually


##collared


output= open("temp/fixed_positions_gt_OC.bed","w")#output file
output.write("scaffold\tposition\tspecies_alleles\n")
total_nb_snps = 0
with open("temp/gt_OC.vcf") as f: # open for collared vcf
    for line in f:
        if not line.startswith("#"):
            total_nb_snps+=1
            if total_nb_snps %100000 == 0: print (total_nb_snps)
            if not "./." in line.split(): # "./." is missing data, so if all individuals have genotype
                dp =line.split()[8].split(":").index("DP") # where is the depth information?
                if all([int(x.split(":")[dp])>=5 for x in line.split()[9:]]) and  len(set([x.split(":")[0] for x in line.split()[9:]]) )== 1:# if all individuals have 5 or more reads and the site is monorphic
                    if line.split()[9].split(":")[0]  == "0/0":# if reference allele
                        output.write (line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[3]+"\n") #write in tht output file
                    elif line.split()[9].split(":")[0] == "1/1": # if alt allele
                        output.write (line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[4]+"\n") #write in tht output file
output.close()


##pied

output= open("temp/fixed_positions_gt_OP.bed","w")#output file
output.write("scaffold\tposition\tspecies_alleles\n")
total_nb_snps = 0
with open("temp/gt_OP.vcf") as f: # open for collared vcf
    for line in f:
        if not line.startswith("#"):
            total_nb_snps+=1
            if total_nb_snps %100000 == 0: print (total_nb_snps)
            if not "./." in line.split(): # "./." is missing data, so if all individuals have genotype
                dp =line.split()[8].split(":").index("DP") # where is the depth information?
                if all([int(x.split(":")[dp])>=5 for x in line.split()[9:]]) and  len(set([x.split(":")[0] for x in line.split()[9:]]) )== 1:# if all individuals have 5 or more reads and the site is monorphic
                    if line.split()[9].split(":")[0]  == "0/0":# if reference allele
                        output.write (line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[3]+"\n") #write in tht output file
                    elif line.split()[9].split(":")[0] == "1/1": # if alt allele
                        output.write (line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[4]+"\n") #write in tht output file
output.close()



#combines the two files
###populate a hash with all collared monomorphic sites
total_lines=0
fixed_diff_coll={}
with open("temp/fixed_positions_gt_OC.bed") as f:
    for line in f:
        total_lines+=1
        if total_lines%100000 ==0: print total_lines
        info = line.split()
        fixed_diff_coll["_".join(info[:2])] = info[2]


### go through the pied differences
fixed_differences = 0
output=open("temp/fixed_differences_Oland.txt","w")
output.write("scaffold\tpos\tcollared_allele\tpied_allele\n")
with open("temp/fixed_positions_gt_OP.bed") as f:
    for line in f:
        info = line.split()
        if fixed_diff_coll.has_key("_".join(info[:2])  ) : # if this site is also fixed in the collared
            if fixed_diff_coll[ "_".join(info[:2]) ] != info[2]:  # if it is fixed for a different allele
                output.write(info[0]+"\t"+info[1]+"\t"+fixed_diff_coll[ "_".join(info[:2]) ]+"\t"+info[2]+"\n")#wite to output : scaffold,pos,collared_allele,,pied_allele
                fixed_differences+=1
                print fixed_differences,info[0],info[1]
output.close()

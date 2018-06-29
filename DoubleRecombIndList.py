
# coding: utf-8

# In[17]:


import numpy as np
import re
import pprint


# In[18]:


bedfile = open('/Users/elizabeth/Dropbox/Projects/MGH-Broad/PGC/PTSD/LAI/FixSwitchErrors/Grady/simualtedTruth/bedfiles/Grady_indivs.txt', 'r')
beds = np.genfromtxt(bedfile, dtype=None)


# In[19]:


beds[:2]


# In[21]:


#when tallying up all the simulated individuals, can just maybe export the total double-crossover counts for each indiv
#rather than exporting all the details for each one now that I know the code works
for bed in beds:
    Atractsfile = open(bed + '.A.bed', 'r')
    Btractsfile = open(bed + '.B.bed', 'r')
    out = open(bed + '.DouRecombs_truth.txt', 'w')
    Atracts = np.genfromtxt(Atractsfile, dtype=None)
    Btracts = np.genfromtxt(Btractsfile, dtype=None)

    for line in range(0, len(Atracts)):
        chromA = Atracts[line][0]
        startA = Atracts[line][5]
        endA = Atracts[line][7]
        ancA = Atracts[line][3]

        for item in range(0, len(Btracts)):
            chromB = Btracts[item][0]
            startB = Btracts[item][5]
            endB = Btracts[item][7]
            ancB = Btracts[item][3]

            count=0

            if chromA == chromB:
                if ancA != ancB:
                    #match to ensure that the ancestries are opposing. Right now, only modeling 2 way admixed scenario

                    #if the start and stop positions are within 10 cM of each other, consider this to be the same location
                    if (startA-10) < startB < (startA+10):
                        if (endA-10) < endB < (endA+10):
                            count = count +1
                        if count > 0:
                            output = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(chromA), str(startA), str(endA), str(ancA), str(startB), str(endB),str(ancB))
                            #output the chromosome number, start and end pos of A and B and their respective ancs for sanity check
                            #these are the double recombination events - the ones where there are matching tracts in opposing ancestries
                            out.write(output)
out.close()


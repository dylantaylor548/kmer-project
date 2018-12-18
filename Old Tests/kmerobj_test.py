
import os
from fractions import Fraction
import copy

# This defines the Kmerobj class. Each kmer will be stored as its own Kmerobj.
# The class stores as its attributes the sequence of the kmer, a list of the locations
# in a FASTA file where that kmer occurs (the index of the first base in the kmer),
# the list of sequences (by sequence name) that contain that kmer, and the worth of the kmer
# The defaults of seqlist and worth are [] and 1 respectively, but are entered as None so
# that those attributes are still mutable.
class Kmerobj(object):

    def __init__(self,kmer,locationlist,seqlist=None,worth=None):
        self.kmer = kmer
        if seqlist is None:
            seqlist = []
        self.seqlist = seqlist
        self.locationlist = locationlist
        if worth is None:
            worth = 1
        self.worth = worth

    def kmermin(self,kmerobj):
        if self.kmer < kmerobj.kmer:
            return self
        else:
            return kmerobj

    def seqlistappend(self,kmerobj):
        seqset = set(self.seqlist)
        seqset = seqset | set(kmerobj.seqlist)
        self.seqlist = list(seqset)

        


# This function essentially converts a sequence of bases into a list of Kmerobj objects
# Such that each kmer in the sequence is represented in the list a single time. If the
# kmer occurs more than once in the sequence, the locations of the multiple instances are
# stored in the locationlist attribute. It also takes the name of the sequence and the
# sequence number (order in the FASTA file) as inputs. This is done so that it can be used
# to generate a list of Kmerobj objects from an entire FASTA file.
def kmerizeseq(seq_name,seqnum,sequence,kmer_size):
    kmer_list = []
    sequence = sequence.upper()
    if (kmer_size <= len(sequence) and kmer_size >= 1):
        for start in range(0,len(sequence)-kmer_size+1):
            kmerseq = sequence[start:start+kmer_size]
            kmerinst = Kmerobj(kmer=kmerseq,locationlist=[start+seqnum],seqlist=[seq_name])
            kmer_list.append(kmerinst)
    return kmer_list
            

# This one is a bit of a doozy. The first part of the code essentially acts as the
# original kmerize_directory function did: reading in a FASTA file and storing any line
# starting with '>' as seq_name and the other lines as sequence. We will go over specifics
# as they appear in the code.
def kmerizefasta(file_path,kmer_size):
    tot_kmer_list = []
    seq_cov_dict = {}
    seq_num = 0
    if file_path.endswith('.fasta'):
        with open(file_path) as opened_file:
            for line in opened_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    seq_name = line[1:]
                    continue
                sequence = line

                # At this point the program has stored the sequence name as 'seq_name'
                # and the sequence itself as 'sequence'. It proceeds to kmerize 'sequence'
                # using the above function and store the resulting list. For each Kmerobj
                # in the list, it goes through the growing total kmer list and checks if any
                # of those Kmerobj objects have the same kmer attribute as the Kmerobj just
                # created. If it finds a match, it will append the seqlist attribute of the
                # recently created Kmerobj to the seqlist attribute of the matched object.
                # It will also append the location list of the recent kmerobj to the matched
                # object. Finally, it will increase the worth of the matched object by 1.
                # If no match is found, it will store the kmerobj as a new object in the total
                # kmer list. It will return the sequence coverage dictionary and the total
                # kmer list as a tuple.
                seq_cov_dict.setdefault(seq_name,0)
                seq_kmer_list = kmerizeseq(seq_name,seq_num*kmer_size,sequence,kmer_size)
                for kmerobj in seq_kmer_list:
                    match = False
                    position = 0
                    while((not match) and position < len(tot_kmer_list)):
                        if kmerobj.kmer == tot_kmer_list[position].kmer:
                            match = True
                        else:
                            position += 1
                    if position == len(tot_kmer_list):
                        tot_kmer_list.append(kmerobj)
                        
                    else:
                        tot_kmer_list[position].seqlistappend(kmerobj)
                        tot_kmer_list[position].locationlist += kmerobj.locationlist
                        tot_kmer_list[position].worth = len(tot_kmer_list[position].seqlist)
                seq_num += 1
                print("We just finished sequence #" + str(seq_num))
        return seq_cov_dict, tot_kmer_list
    else:
        print("That's not a fasta file. Try again.")


# This is analogous to the findkmax function in the original code. Essentially,
# it goes through each kmerobj in a generated list of kmerobj and compares it to
# the currently selected kmax (i.e. the kmer with the highest worth). If a kmerobj
# has a higher worth, it will me make the kmax that kmerobj. Ties are awarded to
# the kmerobj representing more sequences, then lowest lexographically.
def findkmax(tot_kmer_list):
    kmax = None
    for kmerobj in tot_kmer_list:
        if kmax == None:
            kmax = kmerobj
        elif kmerobj.worth > kmax.worth:
            kmax = kmerobj
        elif kmerobj.worth == kmax.worth:
            if len(kmerobj.seqlist) > len(kmax.seqlist):
                kmax = kmerobj
            elif len(kmerobj.seqlist) == len(kmax.seqlist):
                kmax = kmerobj.kmermin(kmax)
    return kmax


# This is analogous to the checkcoverage function in the orginal code. For a given kmerobj,
# it will check each sequence in that kmerobj's seqlist attribute and determine if any of
# them have coverage = cutoff. If all of the sequences have coverage greater than the cutoff
# it will return a value of True. If not, it will return a value of False.
def checkcoverage(kmerobj,seq_cov_dict,cutoffc):
    for seq in seq_cov_dict:
        if seq in kmerobj.seqlist:
            if seq_cov_dict[seq]-1 < cutoffc:
                return False
    return True


# This is analogous to the rep_kmers_indict function in the original code. Given a list of
# kmerobj, a sequence coverage dictionary, kmer size, and cuttoff value (c), it will return
# the "shortest" list of kmerobj that represents each sequence in the dictionary at least c
# times. Additionally, it will eliminate any kmers that overlap eachother in the sequence.
def rep_list_noverlap(tot_kmer_list,seq_cov_dict,kmer_size,cutoffc):
    replist = []
    tot_locs = set()

    while not all(count >= cutoffc for count in seq_cov_dict.values()):
        if tot_kmer_list == []:
            # If the program has gone through all the kmers and hasn't reached coverage for
            # all sequences, it will let the user know the issue.
            print("Your value of c is too high; not enough unique kmers were found.")
            return
        else:
            kmer_max = findkmax(tot_kmer_list)
            # ^^^ this stores the kmerobj with the highest worth in the complete list
            kmax_locs = set(kmer_max.locationlist)
            kmax_range = set()
            for i in range(0,kmer_size):
                kmax_range = kmax_range | set([loc+i for loc in kmer_max.locationlist])
                kmax_range = kmax_range | set([loc-i for loc in kmer_max.locationlist])
            overlaplocs = tot_locs & kmax_range
            # ^^^ this creates a set 
            if overlaplocs == set():
                tot_locs = tot_locs | kmax_locs
                replist.append(kmer_max)
                print("Highest coverage kmer is: '" + kmer_max.kmer + "', with worth: " + str(kmer_max.worth) + " and coverage: " + str(len(kmer_max.seqlist)))
                for seq in kmer_max.seqlist:
                    for kmerobj in tot_kmer_list:
                        if (seq in kmerobj.seqlist) and (seq_cov_dict[seq] < cutoffc) and (kmerobj != kmer_max):
                            kmerobj.worth -= Fraction(1,cutoffc)
                    seq_cov_dict[seq] += 1
                tot_kmer_list.remove(kmer_max)
            else:
                print("kmer overlaps one already chosen")
                tot_kmer_list.remove(kmer_max)

    print("\nChecking output for redundant kmers...")
    replist_copy = copy.deepcopy(replist)
    for kmerobj in replist_copy:
        if checkcoverage(kmerobj,seq_cov_dict,cutoffc):
            print(kmerobj.kmer + " is redundant.")
            replist.remove(kmerobj)

    return replist

#########################################################################################

fastafile = input("Please input the path to the .fasta file of interest ")
c = int(input("Please input the minimum kmer coverage for each of your sequences: "))
kmer_len = int(input("Please select your desired kmer length: "))

coverage_dict, kmerizedfile = kmerizefasta(fastafile,kmer_len)

kmerobj_rep_list = rep_list_noverlap(kmerizedfile,coverage_dict,kmer_len,c)

print("There are " + str(len(kmerobj_rep_list)) + " kmers in the representative list.")

#########################################################################################

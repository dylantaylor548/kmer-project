
import os

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

    def seqlistlen(self):
        return len(self.seqlist)

    def kmermin(self,kmerobj):
        if self.kmer < kmerobj.kmer:
            return self
        else:
            return kmerobj


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
                        tot_kmer_list[position].seqlist.append(kmerobj.seqlist)
                        tot_kmer_list[position].locationlist.append(kmerobj.locationlist)
                        tot_kmer_list[position].worth += 1
                seq_num += 1
                print("We just finished sequence #" + str(seq_num))
        return seq_cov_dict, tot_kmer_list
    else:
        print("That's not a fasta file. Try again.")


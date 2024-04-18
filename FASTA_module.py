#!usr/bin/env python
import sys
#Import everything from the Sequence module
from Sequence_module import *
#Create class called Fasta which inherits the methods from Sequence module
#And is used for reading and getting all the information from fasta files
class Fasta(Sequence):
    '''
    Fasta class
    initialize variables
    write methods for reading fasta file, calculating average GC content of the file
    calculating minimum and maximum length of sequences in file
    '''
    def __init__(self,file_name):
        #initialize file name
        self.file_name = file_name
        # initialize an empty list to contain the sequence objects
        self.fasta_reads = list()


    def read_fasta(self, file_name):
        '''
        Define a method to read a fasta file and store it in an object
        '''
        with open(file_name,'r') as fasta_file:
            ##with the file open in a read mode
            ## initially assign the seq_id value as none
            seq_id = None
            ## loop through every line of the fasta file
            for line in fasta_file:
                #check if the line starts with '>'
                if line.startswith('>'):
                    #if the above condition is true, check if there is a seq_id
                    if seq_id is not None:
                        #if there is a seq_id, then populate existing seq_id to seq_object
                        seq_obj = Sequence(seq_id, seq)
                         #Append the seq_object to a list
                        self.fasta_reads.append(seq_obj)
                        #Assign/reassign the new seq_id to the variable and strip off the new line character
                    seq_id = line.rstrip()
                    ##Assign/reassign empty string to the sequence
                    seq = ''
                ##if the line does not startwith '>', then it is a sequence and it executes the else block
                else:
                    #append the line being read to the seq to include multiple lines of sequences
                    seq = seq+line.rstrip()
            #populate the instances of the Sequence object with seq_id and seq
            seq_obj = Sequence(seq_id, seq)
            #append the seq object to the list
            self.fasta_reads.append(seq_obj)


    def average_gc(self):
        '''
        Define a method to calculate average GC content of the file
        '''
        #initialize the avg_gc value to 0
        avg_gc = 0
        #initialize total gc value to 0
        total_gc = 0
        #loop over the reads in the list of fasta_reads
        for read in self.fasta_reads:
            #get the gc content of each read using the method defined in sequence module
            # add it to avg_gc value and assign it to total gc value
            total_gc = avg_gc + read.get_gc_content()
            #calculate the avg_gc value by dividing the total gc by the length of the fasta_reads
        avg_gc = total_gc / (len(self.fasta_reads))
        #return the avg_gc content value of the file
        return avg_gc

    def avg_length(self):
        '''
        Define a method to find the average length of sequences of a file
        '''
        #initialize the average length value to 0
        avg_len = 0
        #loop over the list of reads
        for read in self.fasta_reads:
            # calculate the average length of sequences by adding all the sequence lengths
            # and dividing it by the total length of the reads
            avg_len = avg_len + read.get_length()
            #return the total counts of the reads in the file
        return avg_len/self.get_read_count()


    def max_length(self):
        '''
        Define a method to identify the sequence of maximum length
        '''
        ## add a break to loop by stating that if the list of fasta_reads does not contain anything
        ## then return none
        if len(self.fasta_reads) == 0:
            return None
        ##initialize the maximum length value to the first sequence of the first object in the list of reads
        max_len = len(self.fasta_reads[1].get_sequence())
        ##loop over the reads in the list
        for read in self.fasta_reads:
            #for every read get the length of the sequence using the getter method defined in main module
            #assign it to a variable
            seq_length = len(read.get_sequence())
            #check if the sequence length is greater than the maximum length
            if seq_length > max_len:
                ##if the condition is true, re-assign the maximum length value to the sequence length
                max_len = seq_length
                #return the value of maximum length of all the sequences
        return max_len


    def min_length(self):
        '''
        Define a method to identify the sequence of minimum length
        '''
        ## add a break to loop by stating that if the list of fasta_reads does not contain anything
        ## then return none
        if len(self.fasta_reads) == 0:
            return None
        ##initialize the minimum length value to the first sequence of the first object in the list of reads
        min_len = len(self.fasta_reads[1].get_sequence())
        #loop over the reads in the list
        for read in self.fasta_reads:
            # for every read get the length of the sequence using the getter method defined in main module
            # assign it to a variable
            seq_len = len(read.get_sequence())
            #check if the sequence length is less than or equal to the minimum length
            if seq_len <= min_len:
                #if the condition is true, re-assign the minimum length value to the sequence length
                min_len = seq_len
                #if the above condition is not true continue searching the reads using the else block
                #if there is no else condition, the execution is stopped after the first condition is
                #found to be false
            else:
                continue
                #return the minimum length of all the sequences
        return min_len


    def get_read_count(self):
        '''
        Define a method to get the read count of all the reads in the list of fasta reads
        '''
        return len(self.fasta_reads)

    def print_metrics(self):
        '''
        Define a method to print the basic metrics of the input file
        '''
        #print the total read count, maximum, minimum and average lengths of sequences
        # of the input file
        print("The total read count is: ", self.get_read_count(), '\n',
              "The max length of seq is: ", self.max_length(), '\n',
              "The min length of seq is: ", self.min_length(), '\n',
              "The avg length of seq is: ", self.avg_length())
        return self.print_metrics


    def write_output(self,op_file_name, len_thresh=None):
        '''
        Define a method to write the output to a file
        If a minimum length threshold is provided
        write all the sequences at or above the threshold to the file
        '''
        # open the output file in the write mode and assign it to a variable
        op_file = open(op_file_name, 'w')
        #loop over the reads in the list
        for read in self.fasta_reads:
            #check if the read length is greater than or equal to the minimum length threshold
            if read.get_length() >= len_thresh:
                #if the condition is true then write the corressponding seq_id and the sequence
                # to the file
                print(read.get_seqid(), file=op_file)
                print(read.get_sequence(), file=op_file)

















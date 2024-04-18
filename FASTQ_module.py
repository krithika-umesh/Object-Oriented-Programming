#!usr/bin/env python
#Import everything from the Sequence module
from Sequence_module import*
#Create class called Fastq which inherits the methods from Sequence module
#And is used for reading and getting all the information from fastq files
class Fastq(Sequence):
    '''
    Fastq class
    initialize variables
    write methods for reading fastq file, calculating average GC content of the file
    calculating minimum and maximum length of sequences in file
    '''
    def __init__(self,file_name):
        #initialize file name
        self.file_name = file_name
        #initialize an empty list to contain the sequence objects
        self.fastq_reads = list()


    def read_fastq(self,file_name):
        '''
        Define a method to read fastq files and store it in an object
        '''
        ## open the file in read mode and assign it to a variable
        fastq_file = open(file_name, 'r')
        #in a while loop go over all sequences until there are none left
        while True:
            #assign the first line of the file to seq_id.
            #replace the '@' with '>' in the file to convert it to fasta format
            #strip off \n character
            seq_id = fastq_file.readline().replace('@','>').rstrip()
            #put a check to break the loop if not seq_id is found
            #which indicates the end of the file
            if not seq_id: break
            #assign the next line to seq and strip off \n
            seq = fastq_file.readline().rstrip()
            #assign the third line to info and strip off \n
            info = fastq_file.readline().rstrip()
            #assign the last line to qual and strip pff the \n
            qual = fastq_file.readline().rstrip()

            #populate the instances of sequence object with
            #seq_id, seq, info, qual
            seq_obj = Sequence(seq_id, seq, info, qual)
            #Append the sequence object to fastq_reads list
            self.fastq_reads.append(seq_obj)


    def average_gc(self):
        '''
        Define a method to calculate average GC content of the file
        '''
        # initialize the avg_gc value to 0
        avg_gc = 0
        # loop over the reads in the list of fastq_reads
        for read in self.fastq_reads:
            # get the gc content of each read using the method defined in sequence module
            # add it to avg_gc value and assign it to total gc value
            total_gc = avg_gc + read.get_gc_content()
            # calculate the avg_gc value by dividing the total gc by the length of the fastq_reads
            avg_gc = total_gc / (len(self.fastq_reads))
            # return the avg_gc content value of the file
        return avg_gc


    def get_read_count(self):
        '''
         Define a method to get the read count of all the reads in the list of fastq reads
        '''
        return len(self.fastq_reads)


    def avg_length(self):
        '''
         Define a method to find the average length of sequences of a file
         '''
        # initialize the average length value to 0
        avg_len = 0
        # loop over the list of reads
        for read in self.fastq_reads:
            # calculate the average length of sequences by adding all the sequence lengths
            # and dividing it by the total length of the reads
            avg_len = avg_len + read.get_length()
            # return the total counts of the reads in the file
        return avg_len / self.get_read_count()


    def max_length(self):
        '''
         Define a method to identify the sequence of maximum length
        '''
        ## add a break to loop by stating that if the list of fasta_reads does not contain anything
        ## then return none
        if len(self.fastq_reads) == 0:
            return None
        ##initialize the maximum length value to the first sequence of the first object in the list of reads
        max_len = len(self.fastq_reads[1].get_sequence())
        ##loop over the reads in the list
        for read in self.fastq_reads:
            # for every read get the length of the sequence using the getter method defined in main module
            # assign it to a variable
            seq_length = len(read.get_sequence())
            # check if the sequence length is greater than the maximum length
            if seq_length > max_len:
                ##if the condition is true, re-assign the maximum length value to the sequence length
                max_len = seq_length
                # return the value of maximum length of all the sequences
        return max_len

    def min_length(self):
        '''
         Define a method to identify the sequence of minimum length
        '''
        ## add a break to loop by stating that if the list of fasta_reads does not contain anything
        ## then return none
        if len(self.fastq_reads) == 0:
            return None
        ##initialize the minimum length value to the first sequence of the first object in the list of reads
        min_len = len(self.fastq_reads[1].get_sequence())
        # loop over the reads in the list
        for read in self.fastq_reads:
            # for every read get the length of the sequence using the getter method defined in main module
            # assign it to a variable
            seq_len = len(read.get_sequence())
            # check if the sequence length is less than or equal to the minimum length
            if seq_len <= min_len:
                # if the condition is true, re-assign the minimum length value to the sequence length
                min_len = seq_len
                # if the above condition is not true continue searching the reads using the else block
                # if there is no else condition, the execution is stopped after the first condition is
                # found to be false
            else:
                continue
                # return the minimum length of all the sequences
        return min_len

    def print_metrics(self):
        '''
          Define a method to print the basic metrics of the input file
        '''
        # print the total read count, maximum, minimum, average lengths and GC content of sequences
        # of the input file
        print("The total read count is: ", self.get_read_count(), '\n',
              "The max length of seq is: ", self.max_length(), '\n',
              "The min length of seq is: ", self.min_length(), '\n',
              "The avg length of seq is: ", self.avg_length(), '\n',
              "The average GC content is: ", self.average_gc())
        return self.print_metrics


    def write_output(self, op_file_name, len_thresh=0,qual_thresh=0):
        '''
         Define a method to write the output to a file
         If a minimum length threshold is provided
         write all the sequences at or above the threshold to the file
        '''
        #assign the default length and quality threshold value to 0
        # open the output file in the write mode and assign it to a variable
        op_file = open(op_file_name, 'w')
        # loop over the reads in the list
        for read in self.fastq_reads:
            # check if the read length is greater than or equal to the minimum length threshold
            #and check if the quality of the read is greater than or equal to the threshold
            if (read.get_qual_score() >= qual_thresh) and (read.get_length() >= len_thresh):
                # if both the conditions are true then write the corresponding seq_id and the sequence
                # to the file
                op_file.write(read.get_seqid() + '\n')
                op_file.write(read.get_sequence() + '\n')








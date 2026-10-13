infile = open("Mef2d_KO/Mef2d_KO_summits.bed") #opens input
outfile = open("Mef2d_KO_summits_great.bed", "w") #creates the output file

for x in infile: #goes through the file on eline at a time

    fields = x.strip().split("\t") #removes spaces and splits colums
    chrom = fields[0] #colum 1 is the chromosome/ sets the chrom
    #keeps chromes 1-19 and chrome x
    if chrom == "chr1" or chrom == "chr2" or chrom == "chr3" or chrom == "chr4" or chrom == "chr5" or chrom == "chr6" or chrom == "chr7" or chrom == "chr8" or chrom == "chr9" or chrom == "chr10" or chrom == "chr11" or chrom == "chr12" or chrom == "chr13" or chrom == "chr14" or chrom == "chr15" or chrom == "chr16" or chrom == "chr17" or chrom == "chr18" or chrom == "chr19" or chrom == "chrX":

        start = int(fields[1]) - 500 #gets start position and subtracts 500
        end = int(fields[2]) + 500 # gets end position and adds 500

        if start < 0: #makes sure it isnt negative
            start = 0

        name = fields[3] #assign name in column 4
        score = float(fields[4]) # assign score in colum 5
        score = round(score) # ROUND and go to integer
        outfile.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\t" + str(score) + "\n") #writes new line in output file

infile.close() #closes input file
outfile.close() # closes outpit file

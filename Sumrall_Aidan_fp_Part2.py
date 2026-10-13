infile = open("greatExportAll.tsv.txt") #opens input great file
outfile = open("Mef2d_targets.csv", "w") #opens output file
genes = {} #creates empty dictionary

for gene in infile: # reads file line by line

    line = gene.strip() #removes extra spaces

    if line == "": #skip a line if it is empty
        continue

    if line.startswith("Region"): # skip headers
        continue

    fields = line.split("\t") # split columns
    peak = fields[0] #gets the peak name
    gene_info = fields[1] # gets gene info
    gene_list = gene_info.split(",") #it splits the gene by comma

    for gene in gene_list: # loops through each gene

        gene_item = gene.strip() # removes extra spaces and characters
        if "(" not in gene_item:
            continue
        parts = gene_item.split("(") #splits at the parathesis
        gene = parts[0].strip() #gets the genes name
        distance = parts[1] # gets distance
        distance = distance.replace(")", "") # removes )
        distance = distance.replace("+", "") # removes the +
        distance = distance.strip() # remove spaces
        distance = int(distance) # convert to an integer

        if abs(distance) <= 2000: #keeps only genes within 2000 base pairs as requested

            if gene not in genes: #checks if gene is already in dictionary

                genes[gene] = peak + "," + str(distance) # if it isnt it saves peak and disance

            else: # if it is in list

                old_info = genes[gene].split(",") #splits old info
                old_distance = int(old_info[1]) #gets the old distance info

                if abs(distance) < abs(old_distance): # compares the distances

                    genes[gene] = peak + "," + str(distance) # replaces it with whichever is closer to the peak

infile.close() # Closes the input file :P
gene_names = list(genes.keys()) #gets gene names
gene_names.sort() #sorts the gene names alphabetically

for x in gene_names: #write the output

    outfile.write(x) #writes gene name
    outfile.write(",") # writes a comma
    outfile.write(genes[x]) # writes the saved peak and distance for that gene
    outfile.write("\n") # new line to move onto next gene

outfile.close() # close output

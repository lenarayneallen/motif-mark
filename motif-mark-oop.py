#!/usr/bin/env python3

import cairo
import argparse

#ARGPARSE
def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", help="FASTA file name", required=True, type=str)
    parser.add_argument("-m", help="Motifs file name", required=True, type=str)
    return parser.parse_args() 

args = get_args()
f = args.f
m = args.m

#DEFINING CLASSES ------------------------------------------------

class InputFile:
    def __init__(self, file_handle):
        '''Creates object of class InputFile'''
        self.filehandle = file_handle
        self.filetype = None
    ## Methods ##
    def get_filetype(self):
        '''Determines if input file is a FASTA file or a motif file'''
        if self.filehandle.endswith(".fasta") == True:
            self.filetype = "fasta"
        elif self.filehandle.endswith(".txt") == True:
            self.filetype = "motifs"
    def generate_motif_list(self):
        '''Generates a list of motifs in a motif file'''
        motiflist = []
        if self.filetype == "motifs":
            #GENERATE MOTIF OBJECTS
            with open(self.filehandle, "r") as mf:
                for line in mf:
                    line = line.strip()
                    motiflist.append(line)
        return(motiflist)                   
    def generate_fasta_dictionary(self):
        '''Generates a dictionary of FASTA sequences from FASTA file
        where keys = gene names and values = [sequence, (chrom, start_pos, end_pos, strandedness)]'''
        if self.filetype == "fasta":
            #GENERATE SEQUENCE OBJECTS
            with open(self.filehandle, "r") as ff:
                fasta_dict = {}
                currentheader = None
                currentsequences = []
                for line in ff:
                    line = line.strip()
                    if line.startswith(">") == True and currentheader == None:
                        #handle > lines
                        currentheader = line.split()
                    elif line.startswith(">") != True:
                        currentsequences.append(line)
                    elif line.startswith(">") == True and currentheader != None:
                        #get gene, chrom, pos, strand from header
                        genename = currentheader[0][1:]
                        chrom_and_loc = currentheader[1].split(":")
                        chrom = chrom_and_loc[0]
                        loc = chrom_and_loc[1].split("-")
                        startpos = loc[0]
                        endpos = loc[1]
                        if len(currentheader) == 4:
                            strandedness = "reverse complement"
                        else:
                            strandedness = ""

                        #concat sequences in current sequences list
                        currentsequences = "".join(currentsequences)

                        #initialize fasta dictionary
                        fasta_dict[genename] = [currentsequences, (chrom, startpos, endpos, strandedness)]
                        
                        
                        #clear currentsequences
                        currentsequences = []
                        currentheader = currentheader = line.split()
                    
                    
                #get gene, chrom, pos, strand from header
                genename = currentheader[0][1:]
                chrom_and_loc = currentheader[1].split(":")
                chrom = chrom_and_loc[0]
                loc = chrom_and_loc[1].split("-")
                startpos = loc[0]
                endpos = loc[1]
                if len(currentheader) == 4:
                    strandedness = "reverse complement"
                else:
                    strandedness = ""

                #concatenate sequences in current sequences list
                currentsequences = "".join(currentsequences)

                #add to fasta dictionary
                fasta_dict[genename] = [currentsequences, (chrom, startpos, endpos, strandedness)]
                
            return(fasta_dict)

class Sequence:
    def __init__(self, sequence, gene_name, chromosome, strandedness, start_pos, end_pos):
        '''Creates an object of class Sequence'''
        self.sequence = sequence
        self.genename = gene_name
        self.chromosome = chromosome
        self.strandedness = strandedness
        self.startpos = start_pos
        self.endpos = end_pos
        self.seqtype = None
        self.introns = None
        self.exon = None
        self.motifs = None


    def find_exon_pos(self):  
        '''Gets 1-based start position of exon. (1-based as this information is only
        used for the pycario drawing, which needs to be to scale) '''
        exonstart = None
        exonend = None
        #find where capital bases start and end
        for i, base in enumerate(self.sequence):
            if base.isupper() == True and exonstart == None:
                exonstart = i+1
            elif exonstart != None and exonend == None and base.islower() == True:
                exonend = i+1
        self.exon_startpos = exonstart
        self.exon_endpos = exonend
    
    def get_exon_objects(self):
        '''generate object of class Exon to hold exon sequence
        saved as attribute of Sequence object (.exon)'''
        self.exon = Exon(self.sequence[self.exon_startpos:self.exon_endpos])
    
    def get_introns(self):
        '''generate intron objects'''
        intronlist = []
        intronlist.append(self.sequence[0:self.exon_startpos])
        intronlist.append(self.sequence[self.exon_endpos:])
        self.introns = intronlist
        
    def find_motifs(self, motifobjs):
        '''Find location of motifs (provided as list of Motif objects) in the sequence'''
        
        #initialize dictionary to hold motifs as keys and (sequence substrings matching motif, (start,end)) as values
        motifdict = {}
        #for motif
        for motifobj in motifobjs:
            motifdict[motifobj.rawmotif] = []
            sequence = self.sequence
            #for base in whole sequence
            for i in range(len(sequence)):
                #get slice the length of the rawmotif
                current_slice = sequence[i:(i+len(motifobj.rawmotif))]
                count = 0
                #iterate through current slice
                for j, base in enumerate(current_slice):
                    #if base in motif is the same as the current base in the current window of the sequence (case insensitive), add 1 to the counter
                    if motifobj.rawmotif[j] == base or motifobj.rawmotif[j].lower() == base or motifobj.rawmotif[j].upper() == base :
                        count += 1
                    #else if concoordant base in the motif is an "ambiguous" base, add 1 to the counter if matching base in the sequence is appropriate
                    elif motifobj.rawmotif[j] in motifobj.ambiguitydict and base in motifobj.ambiguitydict[motifobj.rawmotif[j]]:
                        count += 1
                #if the count is equal to the length of the motif, all bases satisfy the correct conditions, and thus the location of the motif in the sequence can be added to the motif dictionary
                if count == len(motifobj.rawmotif):
                    motifdict[motifobj.rawmotif].append((current_slice, (i,i+len(motifobj.rawmotif))))
        self.motifs = motifdict
    

class Motif():
    def __init__(self,motifseq):
        '''Creates an object of class Motif'''
        self.rawmotif = motifseq
        #True or False
        self.ambiguous_nucleotides = None
    def ambig_nts(self):
        '''Generates a dictionary where keys = all possible ambiguous bases, and values = bases represented by that ambiguous base'''
        ambiguity_dict = {"Y":("C","T","c", "t"), "y":("c","t","C","T"), "W":("A","T","a","t"), "w":("A","T","a","t"),"S":("C","G","c","g"),"s":("C","G","c","g"),"M":("A","C","a","c"),"m":("A","C","a","c"),"K":("G","T","g","t"),"k":("G","T","g","t"),"R":("A","G","a","g"),"r":("A","G","a","g"),"B":("C","G","T","c","g","t"),"b":("C","G","T","c","g","t"),"D":("A","G","T","a","g","t"),"d":("A","G","T","a","g","t"),"H":("A","C","T","a","c","t"),
                          "h":("A","C","T","a","c","t"), "V":("A","C","G","a","c","g"), "v":("A","C","G","a","c","g"), "N":("A","C","G","T","a","c","g","t"),"n":("A","C","G","T","a","c","g","t"), "U":("U","T","u","t"), "u":("U","T","u","t")}
        self.ambiguitydict = ambiguity_dict


class Exon():
    def __init__(self, sequence):
        '''Generates an object of class Exon'''
        self.sequence = sequence
        
class DrawingMaterials():
    def __init__(self,sequencelist):
        '''Generates an object of class DrawingMaterials, used to create final pycario output image'''
        self.sequences_to_draw = sequencelist
    
    def draw_pycario_image(self):
        '''Draws pycario figure'''
        numseqs = len(self.sequences_to_draw)
        imageheight = (numseqs+1) * 200
        seqlengths = [len(seq.sequence) for seq in self.sequences_to_draw]
        longestseq = max(seqlengths)
        imagewidth = longestseq + 500
        WIDTH, HEIGHT = imagewidth,imageheight
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
        context = cairo.Context(surface)
        context.set_source_rgb(1,1,1)
        context.paint()

        #draw intron lines 
        for i,seq in enumerate(self.sequences_to_draw):
            seqlen = len(seq.sequence)
            context.set_line_width(3)
            context.set_source_rgb(0,0,0)
            context.move_to(50,(200*(i+1)))
            context.line_to(seqlen+50, (200*(i+1)))
            context.stroke()
            
        
        #draw exon boxes
        for i,seq in enumerate(self.sequences_to_draw):
            exonstart = seq.exon_startpos
            exonend = seq.exon_endpos
            exonlength = exonend - exonstart

            #exon box
            context.set_source_rgb(0.988, 0.843, 0.612)
            context.rectangle(float(exonstart+50), 175*(i+1)+(25*i), float(exonlength), 50)        
            context.fill()

        #draw color-coded motifs
        for i, seq in enumerate(self.sequences_to_draw):
            #color-coding
            motifs_colors = [(0.561, 0.098, 0.408), (0.937, 0.392, 0.38),(0.247, 0.62, 0.439),(0.341, 0.384, 0.835) ]
            motifs = seq.motifs
            #iterate through motifs
            for j, motif in enumerate(motifs):
                motifcolor = motifs_colors[j]
                context.set_source_rgb(motifcolor[0], motifcolor[1], motifcolor[2])
                for location in motifs[motif]:
                    start = location[1][0]
                    end = location[1][1]
                    context.rectangle(float(start+50), 175*(i+1)+(25*i), float(end-start), 50)        
                    context.fill()
                
        #draw labels
        for i, seq in enumerate(self.sequences_to_draw):
            gene_name = seq.genename
            chrom = seq.chromosome
            chrom = chrom.split("r")
            chrom = chrom[1]
            startpos = seq.startpos
            endpos = seq.endpos
            strand = seq.strandedness


            context.set_source_rgb(0,0,0)
            context.set_font_size(22)
            context.select_font_face("Helvetica",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
            context.move_to(50,150*(i+1)+(50*i))
            if strand != "reverse complement":
                context.show_text(f"{gene_name} chr{chrom}:{startpos} - {endpos}")
            else:
                context.show_text(f"{gene_name} chr{chrom}:{startpos} - {endpos} ({strand})")
            context.stroke()


        #MAKE KEY ------------------------------------------------
        context.set_source_rgb(1, 0.918, 0.788)
        context.rectangle(float(longestseq+100), 135, 300, 260)        
        context.fill()

        #draw key box/title
        context.set_source_rgb(0,0,0)
        context.set_font_size(30)
        context.select_font_face("Helvetica",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        context.move_to(float(longestseq+180),185)
        context.show_text("Motif key:")

        #draw color coded motif boxes
        for i, seq in enumerate(self.sequences_to_draw):
            motifs_colors = [(0.561, 0.098, 0.408), (0.937, 0.392, 0.38),(0.247, 0.62, 0.439),(0.341, 0.384, 0.835) ]
            motifs = seq.motifs
            for j, motif in enumerate(motifs):

                #draw box for key
                context.set_source_rgb(motifs_colors[j][0], motifs_colors[j][1], motifs_colors[j][2])
                context.rectangle(float(longestseq+110), 165+((j+1)*45), float(30), 30)        
                context.fill()

                #draw text for label
                context.set_source_rgb(0,0,0)
                context.set_font_size(24)
                context.select_font_face("Helvetica",cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
                context.move_to(float(longestseq+155), 190+((j+1)*45))
                context.show_text(f"{motif.upper()}")

                
            break

        #write to png output
        outfilename = f[:-6]
        surface.write_to_png(f"{outfilename}.png")


    

#create InputFile object for FASTA file
input_file_fasta = InputFile(f)
input_file_fasta.get_filetype()
fasta_dictionary = input_file_fasta.generate_fasta_dictionary()

#create InputFile object for motif file
input_file_motif = InputFile(m)
input_file_motif.get_filetype()
motif_list = input_file_motif.generate_motif_list()

#create list of Sequence objects 
sequenceobjects = [Sequence(fasta_dictionary[key][0], key, fasta_dictionary[key][1][0],fasta_dictionary[key][1][3],  fasta_dictionary[key][1][1], fasta_dictionary[key][1][2]) for key in fasta_dictionary]
#create list of Motif objects
motifobjects = [Motif(x) for x in motif_list]

#get ambiguous base dictionaries
for m_obj in motifobjects:
    m_obj.ambig_nts()

#process sequence objects
for s_obj in sequenceobjects:
    #find exon positions
    s_obj.find_exon_pos()
    #get exon objects
    s_obj.get_exon_objects()
    #get introns
    s_obj.get_introns()
    #find motifs
    s_obj.find_motifs(motifobjects)

#use list of processed sequence objects as input to create an object of class DrawingMaterials
drawingstuff = DrawingMaterials(sequenceobjects)
#draw the pycario figure
drawingstuff.draw_pycario_image()


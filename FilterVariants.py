'''
Created by Ryan O. Schenck
Date: 14 August 2018

Settings: rescomp1 python3.6 interpreter
Deployment: Mapping to /well/leedham/users/rschenck/AbbyMice/VariantCalls
Deploy on cluster: Tools -> Deployment -> Configuration (Set mapping in here)
'''

import sys
import os
import argparse
import logging
import glob
import subprocess
import gzip
import pysam
try:
    from pyfiglet import Figlet
except:
    pass



def VariantFilterBanner():
    try:
        print('\n')
        f = Figlet(font='colossal')
        print(f.renderText('Variant Filter'))
    except:
        pass

def Parser():
    # get user variables
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="inDir", default='./', type=str, help="Input directory of vcf files with .passed.vcf extension, in the case of Strelka2 it expects the directory containing all information. Default: ./")
    parser.add_argument("-c", dest="caller", default=1, type=int, help="Callers to be used, 1) MuTect2 2) Strelka2 3) Octopus 4) Shearwater ML 5) DRAGEN; Default: 1) MuTect2")
    parser.add_argument("-ievs", dest="ignoreLowEVS", default=False, action="store_true", help="Whether to ignore the LowEVS filter on Strelka2 output. Default=False. (Only applies for Strelka2 caller)")
    parser.add_argument("-nr", dest="minnr", default=10, type=int, help="Minimum reads in normal. Default: 10")
    parser.add_argument("-vr", dest="minvr", default=2, type=int, help="Minimum variant reads. Default: 2")
    parser.add_argument("-vd", dest="varsamdepth", default=10, type=int, help="Total reads at variant site sum(ref,alt). Default=10")
    parser.add_argument("-af", dest="minvaf", default=0.1, type=float, help="Minimum variant allele frequency. Default=0.1")
    parser.add_argument("-e", dest="maxevents", default=1, type=int, help="Maximum number of events allowed at a loci. This ignores germline for example. Default=1") # TODO change code to handle germline variants and subsequent mutation. Works for MuTect currently.
    parser.add_argument("-l", dest="clusterFlank", default=10, type=int, help="Flanking region to omit clustered events. Default: 10")
    parser.add_argument("-t", dest="table", default=False, action="store_true", help="Flag to build an informative table. Useful for quick R processing without additional Bioconductor packages.")
    parser.add_argument("-u", "--usesamplemap", dest="usesamplemap", default=False, action="store_true",help="Boolean indicating whether to use same sample filter options. Default: True, must specify --samplemap <file>.")
    parser.add_argument("--samplemap", default='./nonenonenonenone?!', dest='samplemap', type=str, help="File with sample mapping in comma separated text file. Default: None, required if -u ")
    parser.add_argument("-uvr", dest="uminvr", default=2, help="Minimum variant reads if -u setting is present. Default: 2")
    parser.add_argument("-uaf", dest="uminvaf", default=0.01, type=float, help="Minimum variant allele frequency if -u setting is present. Default: 0.005")
    parser.add_argument("-v", dest="verbose", default=True, action="store_false", help="Verbosity setting. Use -v to suppress some stdout.")
    parser.add_argument("-a", dest="ann", default=False, action="store_true", help="Annotate using snpEff. This requires -ref and -j. Default=F")
    parser.add_argument("-ref", dest="ref", default=None, type=str, help="Reference genome to be used for snpEff annotations (See http://snpeff.sourceforge.net/index.html). Default=None")
    parser.add_argument("-j", dest="jarfile", default=None, type=str, help="snpEff.jar executable (See http://snpeff.sourceforge.net/index.html) Default=None")
    parser.add_argument("-reffasta", dest="refgenome", default='/well/leedham/users/rschenck/References/Ensembl/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa', type=str, help="Reference genome that the variants are called with. Default=/well/leedham/users/rschenck/References/Ensembl/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa")
    Options = parser.parse_args()
    # Run options tests
    if Options.inDir[-1] is not '/':
        Options.inDir += '/'
    if os.path.isdir(Options.inDir)==False:
        logging.error("Unable to locate input directory. Please check -i option.")
        sys.exit()
    else:
        logging.info("Input directory found: %s"%(Options.inDir))
    if Options.usesamplemap and os.path.exists(Options.samplemap) == False:
        logging.error("If --usesamplemap/-u is set then a sample map must be provided.")
        sys.exit()
    elif Options.usesamplemap==False:
        logging.info("No sample map provided. Will not apply same specimen filters...")
    else:
        logging.info("Performing joint filtering.")
    if Options.minvaf<0.0 or Options.minvaf>1. or Options.uminvaf<0.0 or Options.uminvaf>1.0:
        logging.error("Variant allele frequency options (-af or -uaf) are outside acceptable ranges.")
        sys.exit()
    elif Options.verbose:
        logging.info("Minimum variant reads: %s"%(Options.minvr))
        logging.info("Minimum normal reads: %s"%(Options.minnr))
        logging.info("Minimum variant allele frequency: %s"%(Options.minvaf))
        logging.info("Minimum clustered events: %sbp"%(Options.clusterFlank))
        if Options.usesamplemap:
            logging.info("Minimum variant reads in at least one sample: %s"%(Options.uminvr))
            logging.info("Minimum variant allele freq in at least one sample: %s"%(Options.uminvaf))
        if Options.ann:
            logging.info("Will attempt to perform annotations using snpEff using reference %s"%(Options.ref))
    if Options.ann:
        if os.path.exists(Options.jarfile)==False or Options.ref is None:
            logging.error("If -a is specified a jarfile (-j) and reference (-ref) for snpEff must be used.")
            sys.exit()
    logging.info("Beginning...")
    return(Options)

def SetupEnvironment():
    '''
    Setup logging function and print banner.
    :return: None
    '''
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s-%(levelname)s: %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)
    logging.info("Thank you for using Variant Filter. Created by Ryan Schenck...")
    VariantFilterBanner()

class MuTect2Sample:

    def __init__(self, vcf, Options, theFlags):
        logging.info("Setting filter flags for %s..."%(os.path.basename(vcf)))
        self.flagentries = theFlags
        self.sampleidentifier = os.path.basename(vcf)
        self.samplePath = os.path.dirname(vcf) + '/'
        with open(vcf, 'r') as inFile:
            self.passedlines = [line.replace('\n','') for line in inFile]
        # Step 1: Read in muts and get tumor and normal column indices
        self.normal_col = None
        self.tumor_col = None
        self.header,self.rawmuts = self._ParseFileAndBuildHeader(Options)
        # Step 1a: Get positions for cluster events
        self.chromsPos = PrepareClusteredEventCheck(self.rawmuts)

        # Step 2: Construct flags for specific attributes to filter on.
        self.flags = ['']*len(self.rawmuts)

        # Step 3: Set flags for normal filtering (@,?,*,^)
        self._normalcheck(Options)

        # Step 4: Set flags for variant filtering (&,!,`)
        self._variantcheck(Options)

        if Options.usesamplemap:
            self.jointsampleref = self._constructjointrefdict()

    def _ParseFileAndBuildHeader(self, Options):
        argparsevars = {'minnr': '-nr', 'minvr': '-vr', 'minvaf': '-af', 'clusterFlank': '-l', 'table': '-t',
                        'inDir': '-i', 'caller': '-c', 'usesamplemap': '--usesamplemap', 'samplemap': '--samplemap',
                        'uminvr': '-uvr', 'uminvaf': '-uaf', 'verbose': '-v', 'maxevents':'-e', 'varsamdepth':'-vd',
                        'ann':'-a', 'ref':'-ref','jarfile':'-jar', 'refgenome':'-reffasta', 'ignoreLowEVS':'-ievs'}
        commandLine = 'FilterVariants.py ' + ' '.join([' '.join([argparsevars[arg], str(getattr(Options, arg))]) for arg in vars(Options)])
        outline = '##FilterVariants=<ID=FilterVariants,CommandLine="%s">'%(commandLine)
        newHeader = ''
        rawMuts = []
        for line in self.passedlines:
            if line.startswith('#'):
                if line.startswith('##FORMAT=<ID=SA_POST_PROB'):
                    newHeader+=outline+'\n'
                elif line.startswith('##normal_sample'):
                    normcolumn = line.split('=')[1]
                    newHeader+=line+'\n'
                elif line.startswith('##tumor_sample'):
                    tumorcolumn = line.split('=')[1]
                    newHeader+=line+'\n'
                elif line.startswith('#CHROM'):
                    chromline = line.split('\t')
                    newHeader+=line+'\n'
                else:
                    newHeader+=line+'\n'
            else:
                rawMuts.append(line)

        # Determine and set the column number for normal and blood...since mutect randomly switches this shit! :(
        self.normal_col = chromline.index(normcolumn)
        self.tumor_col = chromline.index(tumorcolumn)

        return(newHeader, rawMuts)

    def _normalcheck(self, Options):
        normalchecks = {'@':0,'?':0,'*':0,'^':0}
        for i,mut in enumerate(self.rawmuts):
            m = mut.split('\t')
            mnorm = dict(zip(m[8].split(':'), m[self.normal_col].split(':')))
            if int(mnorm['AD'].split(',')[1]) > 0:
                self.flags[i] += '@'
                normalchecks['@']+=1
            if int(mnorm['AD'].split(',')[0]) < Options.minnr:
                self.flags[i] += '*'
                normalchecks['*']+=1
            if min([float(af) for af in mnorm['AF'].split(',')]) > 0.05:
                self.flags[i] += '?'
                normalchecks['?']+=1
            if len(m[4].split(','))>1:
                self.flags[i] += '^'
                normalchecks['^']+= 1

        if Options.verbose:
            logging.info('Normal sample check metrics...')
            for item in normalchecks:
                if normalchecks[item] > 0:
                    logging.info('%s: %s (%.1f%%)'%(self.flagentries[item], normalchecks[item], normalchecks[item]/float(len(self.rawmuts))*100))

    def _variantcheck(self, Options):
        variantchecks = {'&':0,'!':0,'`':0,';':0}
        for i,mut in enumerate(self.rawmuts):
            m = mut.split('\t')
            mvar = dict(zip(m[8].split(':'), m[self.tumor_col].split(':')))
            # check overall depth Options.varsamdepth
            if sum([int(d) for d in mvar['AD'].split(',')]) < Options.varsamdepth:
                self.flags[i] += '`'
                variantchecks['`']+=1
            # check variant read depth Options.minvr
            if int(mvar['AD'].split(',')[1]) < Options.minvr and len(mvar['AD'].split(','))==2:
                self.flags[i] += '&'
                variantchecks['&'] += 1
            elif len(mvar['AD'].split(','))>2:
                checkmultipleEvents = [int(val) for val in mvar['AD'].split(',')[1:]]
                if min(checkmultipleEvents) < Options.minvr:
                    self.flags[i] += '&' # TODO Need to check for one of these events (>1 event) being true with proper support later on...
                    variantchecks['&'] += 1
            elif int(mvar['AD'].split(',')[1]) >= Options.minvr and len(mvar['AD'].split(','))==2:
                pass

            # check variant allele frequency Options.minvaf
            if len(mvar['AD'].split(','))==2 and int(mvar['AD'].split(',')[1])/float(int(mvar['AD'].split(',')[0])+int(mvar['AD'].split(',')[1])) < Options.minvaf:
                self.flags[i] += '!'
                variantchecks['!'] += 1
            elif len(mvar['AD'].split(','))>2:
                checkmultipleEvents = [int(val) for val in mvar['AD'].split(',')[1:]]
                totalReads = sum(checkmultipleEvents) + int(mvar['AD'].split(',')[0])
                multipleVAFs = [val/float(totalReads) for val in checkmultipleEvents]
                if min(multipleVAFs) < Options.minvaf:
                    self.flags[i] += '!'
                    variantchecks['!'] += 1

            # check clustered event Options.clusterFlank
            if ClusterMutCount(int(m[1]),self.chromsPos[m[0]],Options) > 0:
                self.flags[i] += ';'
                variantchecks[';'] += 1

        if Options.verbose:
            logging.info('Variant sample check metrics...')
            for item in variantchecks:
                if variantchecks[item] > 0:
                    logging.info('%s: %s (%.1f%%)'%(self.flagentries[item], variantchecks[item], variantchecks[item]/float(len(self.rawmuts))*100))

    def _constructjointrefdict(self):
        '''
        Constructs a dictionary with the following information {mutpos: [index,flags]}
        :return: a nested dictionary
        '''
        theRef = {}
        for i,m in enumerate(self.rawmuts):
            key = m.split('\t')[0] + ':' + m.split('\t')[1] + '_' + m.split('\t')[3] + '/' + m.split('\t')[4]
            theRef.update({key:[i,self.flags[i]]})
        return(theRef)

    def _applyCNVcorrection(self):
        #TODO once CNV calls have been made and implemented
        pass

class MuTect2Filter:

    def __init__(self, Options, inFiles, theFlags):
        self.allMutectSamples = []
        # Step 1: Perform individual sample based filtering
        self._ApplyMutect2Filter(Options, inFiles, theFlags)

        # Step 2: Perform joint based filtering
        if Options.usesamplemap:
            self._Mutect2JointFiltering(Options, inFiles[0], theFlags)
        else:
            pass

        logging.info("Writing filtered vcf files...")
        for sam in self.allMutectSamples:
            WriteFilteredVCFFiles(sam.samplePath, sam.sampleidentifier, sam.header, sam.rawmuts, sam.flags)
        if Options.ann:
            logging.info("Annotating vcf files using snpEff...")
            for sam in self.allMutectSamples:
                AnnotateVCFFiles(sam.samplePath, sam.sampleidentifier, Options)

        if Options.table:
            self._BuildMutect2Table(Options)

    def _ApplyMutect2Filter(self, Options, inFiles, theFlags):
        if Options.usesamplemap:  # with sample mapping
            for theFile in inFiles[1]:
                self.allMutectSamples.append(MuTect2Sample(theFile, Options, theFlags))
        else:  # no sample mapping
            for theFile in inFiles:
                self.allMutectSamples.append(MuTect2Sample(theFile, Options, theFlags))

    def _Mutect2JointFiltering(self, Options, samplemap, theFlags):
        logging.info('...')
        logging.info('...')
        logging.info("Performing joint filtering by updating filter flags...")
        # Responsible for modifying the flags for each of the variant filtered sets for the final output.
        for sam in samplemap:
            if len(sam)>1:
                sampleIDs = []
                for indsam in sam:
                    for i, sample in enumerate(self.allMutectSamples):
                        if indsam in sample.sampleidentifier:
                            sampleIDs.append(i)

                # Perform joint filter flag adjustments

                #~~~~~ Get Candidate Sites ~~~~#
                # logging.info("Gathering shared mutation loci within specimen...")
                # Step 1: loop over each sample set of positions w/ mutations and find the matching ones
                matchingPositions = {}
                for theid in sampleIDs:
                    for pos in self.allMutectSamples[theid].jointsampleref:
                        try:
                            matchingPositions[pos] += 1
                        except KeyError:
                            matchingPositions.update({pos:1})
                # Step 2: Trim to shared positions w/ muations to examine
                candidateSites = []
                for pos in matchingPositions:
                    if matchingPositions[pos] > 1:
                        candidateSites.append(pos)

                for pos in candidateSites:
                    flagsets = []
                    for theid in sampleIDs:
                        try:
                            # Must do it this way to avoid errors later on. This means that >2 samples in a specimen and mutation shared between < total samples
                            a = theid
                            b = self.allMutectSamples[theid].sampleidentifier
                            c = self.allMutectSamples[theid].jointsampleref[pos]
                            d = self.allMutectSamples[theid].rawmuts[self.allMutectSamples[theid].jointsampleref[pos][0]]
                            flagsets.append(a) # id of the sample inside of self.allMutectSamples
                            flagsets.append(b) # Get Identifier
                            flagsets.append(c) # Get Index and Flagset
                            flagsets.append(d) # Get full mutation information
                        except KeyError:
                            pass

                    # Check if a flagset needs to be jointly filtered
                    checkcount = 0
                    flagspresent = []
                    for z in range(0,len(flagsets),4):
                        if flagsets[z+2][1]!='':
                            checkcount+=1
                            flagspresent.append(flagsets[z+2][1])

                    if checkcount > 0:
                        # Check if all flags are the same within the first (cond1 and cond2) or truly different...then adjust flags
                        if (len(flagspresent)==len(flagsets)/4 and len(list(set(flagspresent)))!=1) or len(flagspresent)==1:
                            if any(i in '.'.join(flagspresent) for i in ['!','&']):
                                # ~~~~ Adjust flags ~~~~#
                                allVafs = []
                                allnvr = []
                                for z in range(0, len(flagsets), 4):
                                    m = flagsets[z+3].split('\t')
                                    mvar = dict(zip(m[8].split(':'), m[self.allMutectSamples[flagsets[z]].tumor_col].split(':')))
                                    allnvr.append(int(mvar['AD'].split(',')[1]))
                                    totaldepth = float(sum([int(val) for val in mvar['AD'].split(',')]))
                                    allVafs.append(int(mvar['AD'].split(',')[1])/totaldepth)
                                # Step 1: Adjust for variant supporting reads flags (&)
                                if all(vr >= Options.uminvr for vr in allnvr) and any(vr>= Options.minvr for vr in allnvr):
                                    for z in range(0, len(flagsets), 4):
                                        self.allMutectSamples[flagsets[z]].flags[flagsets[z+2][0]] = self.allMutectSamples[flagsets[z]].flags[flagsets[z+2][0]].replace('&','')
                                # Step 2: Adjust for variant allele frequencies flags (!)
                                if all(vaf >= Options.uminvaf for vaf in allVafs) and any(vaf >= Options.minvaf for vaf in allVafs):
                                    for z in range(0, len(flagsets), 4):
                                        self.allMutectSamples[flagsets[z]].flags[flagsets[z+2][0]] = self.allMutectSamples[flagsets[z]].flags[flagsets[z+2][0]].replace('!','')
                            else:
                                pass # Don't do any flag adjustments as none of the joint filtering options apply

        logging.info("Filter flags have been updated...")
        logging.info("Joint filtering completed...")
        logging.info('...')
        logging.info('...')

    def _BuildMutect2Table(self, Options):
        if Options.ann==False:
            fullFail = ['@','*','&','!','^','`',';']  # All except '?'

            allMuts = []
            for sam in self.allMutectSamples:
                with open(sam.samplePath + sam.sampleidentifier.replace("passed.vcf","passed.variantfilter.vcf"), 'r') as inFile:
                    lines = inFile.readlines()
                for line in lines:
                    if line.startswith('#')==False:
                        allMuts.append(line)
            trinucs = GetTrinucs(allMuts, Options)

            outTableName = self.allMutectSamples[0].samplePath + 'Mutect2.passed.variantfilter.digest.txt'
            header = "sample\tchr\tpos\tref_nt\tref_nt_pred\tmut_nt\treads\treads_fw\treads_rv\tdepth\tvaf\ttrinuc\tgene_name\timpact\tprotein_change\tnt_change\tprotein_id\n"
            outLines = []
            for sam in self.allMutectSamples:
                sample = sam.sampleidentifier.replace('.passed.vcf','')
                for i, flagIDs in enumerate(sam.flags):
                    if any(f in flagIDs for f in fullFail) == False:
                        m = sam.rawmuts[i].split('\t') # Get the mutation
                        # Check for multiple events
                        if len(m[4].split(',')) > 1:
                            for z in range(0,len(m[4].split(','))):
                                chr = m[0]
                                pos = m[1]
                                ref_nt = m[3]
                                ref_nt_pred = m[3]
                                mut_nt = m[4].split(',')[z]
                                mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
                                reads = int(mvar['AD'].split(',')[z+1])
                                reads_fw = int(mvar['F1R2'].split(',')[z+1])
                                reads_rv = int(mvar['F2R1'].split(',')[z+1])
                                depth = sum([int(reads) for reads in mvar['AD'].split(',')])
                                vaf = reads / float(depth)
                                if len(mut_nt) == 1 and len(ref_nt) == 1:
                                    try:
                                        trinuc = trinucs[str(chr) + ':' + str(pos)]
                                    except KeyError:
                                        trinuc = '.'
                                        logging.warning("Unable to locate trinucleotide context...")
                                else:
                                    trinuc = '.'
                                gene_name = '.'
                                impact = '.'
                                protein_change = '.'
                                nt_change = '.'
                                protein_id = '.'
                                l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
                                     vaf, trinuc, gene_name, impact, protein_change, nt_change, protein_id]
                                strLine = [str(item) for item in l]
                                theLine = '\t'.join(strLine)
                                outLines.append(theLine)
                        else:
                            chr = m[0]
                            pos = m[1]
                            ref_nt = m[3]
                            ref_nt_pred = m[3]
                            mut_nt = m[4]
                            mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
                            reads = int(mvar['AD'].split(',')[1])
                            reads_fw = int(mvar['F1R2'].split(',')[1])
                            reads_rv = int(mvar['F2R1'].split(',')[1])
                            depth = sum([int(reads) for reads in mvar['AD'].split(',')])
                            vaf = reads/float(depth)
                            if len(mut_nt)==1 and len(ref_nt)==1:
                                try:
                                    trinuc = trinucs[str(chr)+':'+str(pos)]
                                except KeyError:
                                    trinuc = '.'
                                    logging.warning("Unable to locate trinucleotide context...")
                            else:
                                trinuc = '.'
                            gene_name = '.'
                            impact = '.'
                            protein_change = '.'
                            nt_change = '.'
                            protein_id = '.'
                            l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
                                 vaf, trinuc, gene_name, impact, protein_change, nt_change, protein_id]
                            strLine = [str(item) for item in l]
                            theLine = '\t'.join(strLine)
                            outLines.append(theLine)
            with open(outTableName, 'w') as outFile:
                outFile.write(header)
                outFile.write('\n'.join(outLines) + '\n')
        else: # If annotations are done...
            # Get trinucleotide context information
            allMuts = []
            for sam in self.allMutectSamples:
                with open(sam.samplePath + sam.sampleidentifier.replace("passed.vcf", "passed.variantfilter.ann.vcf"),
                          'r') as inFile:
                    lines = inFile.readlines()
                for line in lines:
                    if line.startswith('#') == False:
                        allMuts.append(line)
            trinucs = GetTrinucs(allMuts, Options)

            annMuts = {}
            for sam in self.allMutectSamples:
                with open(sam.samplePath + sam.sampleidentifier.replace("passed.vcf","passed.variantfilter.ann.vcf"), 'r') as inFile:
                    lines = [line.replace('\n','') for line in inFile.readlines()]
                allLines = []
                for line in lines:
                    if line.startswith('#') == False:
                        allLines.append(line)
                annMuts.update({sam.sampleidentifier:allLines})

            outTableName = self.allMutectSamples[0].samplePath + 'Mutect2.passed.variantfilter.ann.digest.txt'
            header = "sample\tchr\tpos\tref_nt\tref_nt_pred\tmut_nt\treads\treads_fw\treads_rv\tdepth\tvaf\ttrinuc\tgene_name\timpact\tprotein_change\tnt_change\ttranscript_id\tgene_id\n"
            outLines = []
            for sam in self.allMutectSamples:
                sample = sam.sampleidentifier.replace('.passed.vcf','')
                muts = annMuts[sam.sampleidentifier]
                for mut in muts:
                    m = mut.split('\t')  # Get the mutation
                    info = {}
                    for item in m[7].split(';'):
                        infopair = item.split('=')
                        if len(infopair) != 2:
                            info.update({infopair[0]:'.'})
                        else:
                            info.update({infopair[0]:infopair[1]})
                    info = info['ANN'].split(',')[0]
                    # Check for multiple events
                    if len(m[4].split(',')) > 1:
                        for z in range(0, len(m[4].split(','))):
                            chr = m[0]
                            pos = m[1]
                            ref_nt = m[3]
                            ref_nt_pred = m[3]
                            mut_nt = m[4].split(',')[z]
                            mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
                            reads = int(mvar['AD'].split(',')[z + 1])
                            reads_fw = int(mvar['F1R2'].split(',')[z + 1])
                            reads_rv = int(mvar['F2R1'].split(',')[z + 1])
                            depth = sum([int(reads) for reads in mvar['AD'].split(',')])
                            vaf = reads / float(depth)
                            if len(mut_nt) == 1 and len(ref_nt) == 1:
                                try:
                                    trinuc = trinucs[str(chr) + ':' + str(pos)]
                                except KeyError:
                                    trinuc = '.'
                                    logging.warning("Unable to locate trinucleotide context...")
                            else:
                                trinuc = '.'
                            gene_name = info.split('|')[3]
                            impact = info.split('|')[1].split('&')[0]
                            protein_change = info.split('|')[10]
                            nt_change = info.split('|')[9]
                            transcript_id = info.split('|')[6]
                            gene_id = info.split('|')[4]
                            l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
                                 vaf, trinuc, gene_name, impact, protein_change, nt_change, transcript_id,gene_id]
                            strLine = [str(item) for item in l]
                            theLine = '\t'.join(strLine)
                            outLines.append(theLine)
                    else:
                        chr = m[0]
                        pos = m[1]
                        ref_nt = m[3]
                        ref_nt_pred = m[3]
                        mut_nt = m[4]
                        mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
                        reads = int(mvar['AD'].split(',')[1])
                        reads_fw = int(mvar['F1R2'].split(',')[1])
                        reads_rv = int(mvar['F2R1'].split(',')[1])
                        depth = sum([int(reads) for reads in mvar['AD'].split(',')])
                        vaf = reads / float(depth)
                        if len(mut_nt) == 1 and len(ref_nt) == 1:
                            try:
                                trinuc = trinucs[str(chr) + ':' + str(pos)]
                            except KeyError:
                                trinuc = '.'
                                logging.warning("Unable to locate trinucleotide context...")
                        else:
                            trinuc = '.'
                        gene_name = info.split('|')[3]
                        impact = info.split('|')[1].split('&')[0]
                        protein_change = info.split('|')[10]
                        nt_change = info.split('|')[9]
                        transcript_id = info.split('|')[6]
                        gene_id = info.split('|')[4]
                        l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
                             vaf, trinuc, gene_name, impact, protein_change, nt_change, transcript_id, gene_id]
                        strLine = [str(item) for item in l]
                        theLine = '\t'.join(strLine)
                        outLines.append(theLine)
                with open(outTableName, 'w') as outFile:
                    outFile.write(header)
                    outFile.write('\n'.join(outLines) + '\n')

class Strelka2Sample:

    def __init__(self, Options, theFlags, SamNam, germDir, somDir):
        self.SamNam = SamNam
        self.GermlineDir = germDir
        self.somDir = somDir
        self.bamFile = None
        self.rawGermMuts = []
        self.filteredGermMuts = []
        self.rawGermMutHead = []
        self.rawGermLineFlags = []
        self.rawSNVS = []
        self.rawSNVSHead = []
        self.SNVSHeadout = ''
        self.rawSNVSflags = []
        self.rawIndels = []
        self.rawIndelsHead = []
        self.rawIndelsflags = []
        self.rawPassIndels = []

        # Step 1: Read in all corresponding files for GermMuts and Somatic SNVs
        self._readFiles()
        self._newHeader(Options)

        # Step 1a: Get positions for cluster events
        self.chromsPos = PrepareClusteredEventCheck(self.rawSNVS)
        self.chromsPosGerm = PrepareClusteredEventCheck(self.rawGermMuts)

        # Step 2: Perform Filtering
        logging.info("Filtering calls for %s..." % (self.SamNam))
        logging.info("GERMLINE filtering...")
        self._filterRawGermlineVariant(Options, theFlags)
        logging.info("SOMATIC filtering...")
        self._filterRawSomaticSNVs(Options, theFlags)
        self._filterRawIndelSNVs(Options, theFlags)

        if Options.usesamplemap:
            self.germjointsampleref = self._constructjointrefdictGerm()
            self.snvjointsampleref = self._constructjointrefdictSomatic()

    def _readFiles(self):
        germFile = self.GermlineDir + '/results/variants/variants.vcf.gz'
        somSNVS = self.somDir + '/results/variants/somatic.snvs.vcf.gz'
        somINDELS = self.somDir + '/results/variants/somatic.indels.vcf.gz'

        logging.info("Reading somatic SNVs...")
        with gzip.open(somSNVS, 'rb') as inSom:
            SomLines = [line.decode('UTF-8').replace('\n', '') for line in inSom.readlines()]
        for line in SomLines:
            if line.startswith('#'):
                self.rawSNVSHead.append(line)
            else:
                self.rawSNVS.append(line)

        logging.info("Reading somatic Indels...")
        with gzip.open(somINDELS, 'rb') as inSomIndels:
            SomLinesIndels = [line.decode('UTF-8').replace('\n', '') for line in inSomIndels.readlines()]
        for line in SomLinesIndels:
            if line.startswith('#'):
                self.rawIndelsHead.append(line)
            else:
                self.rawIndels.append(line)

        logging.info("Reading germline variant calls...")
        with gzip.open(germFile, 'rb') as inGerm:
            GermLine = [line.decode('UTF-8').replace('\n', '') for line in inGerm.readlines()]
        for line in GermLine:
            if line.startswith('#'):
                self.rawGermMutHead.append(line)
            else:
                if line.split('\t')[6] == "PASS":
                    self.rawGermMuts.append(line)

    def _newHeader(self, Options):
        argparsevars = {'minnr': '-nr', 'minvr': '-vr', 'minvaf': '-af', 'clusterFlank': '-l', 'table': '-t',
                        'inDir': '-i', 'caller': '-c', 'usesamplemap': '--usesamplemap', 'samplemap': '--samplemap',
                        'uminvr': '-uvr', 'uminvaf': '-uaf', 'verbose': '-v', 'maxevents': '-e', 'varsamdepth': '-vd',
                        'ann': '-a', 'ref': '-ref', 'jarfile': '-jar', 'refgenome': '-reffasta',
                        'ignoreLowEVS': '-ievs'}
        commandLine = 'FilterVariants.py ' + ' '.join(
            [' '.join([argparsevars[arg], str(getattr(Options, arg))]) for arg in vars(Options)])
        outline = '##FilterVariants=<ID=FilterVariants,CommandLine="%s">' % (commandLine)
        newHeader = ''
        for line in self.rawSNVSHead:
            if line.startswith('#'):
                if line.startswith('##reference='):
                    newHeader+=outline+'\n'
                    newHeader += line + '\n'
                elif line.startswith("##cmdline="):
                    self.bamFile = line.split("--tumorBam=")[1].split(" ")[0]
                    newHeader += line + '\n'
                else:
                    newHeader+=line+'\n'
        self.SNVSHeadout = newHeader

    def _filterRawGermlineVariant(self, Options, theFlags):
        vafs = []
        candidates = 0
        checks = {'^':0, '`':0, '&':0, '!':0, ';':0}
        self.rawGermLineFlags = ['' for i in range(0,len(self.rawGermMuts))]
        for i, line in enumerate(self.rawGermMuts):
            # if 'PASS' == line.split('\t')[6]:
            line = line.split('\t')

            # Filter multiple events at same loci
            if len(line[4].split(','))>1:
                self.rawGermLineFlags[i] += '^'
                checks['^'] += 1
            field = line[9].split(':')
            # Filter depth less than
            if int(field[3])<Options.minnr:
                self.rawGermLineFlags[i] += '`'
                checks['`'] += 1
            # Filter Variant Depth
            if int(field[5].split(',')[1]) < 2:
                self.rawGermLineFlags[i] += '&'
                checks['&'] += 1
            elif int(field[6].split(',')[1]) < 2 or int(field[6].split(',')[1]) < 2:
                self.rawGermLineFlags[i] += '&'
                checks['&'] += 1
            else:
                pass
            if int(field[3]) > 0:
                if int(field[5].split(',')[1])/float(field[3]) < 0.1:
                    self.rawGermLineFlags[i] += '!'
                    checks['!'] += 1

            # check clustered event Options.clusterFlank
            if ClusterMutCount(int(line[1]), self.chromsPosGerm[line[0]], Options) > 0:
                self.rawGermLineFlags[i] += ';'
                checks[';'] += 1

            # Total candidates
            if self.rawGermLineFlags[i]=='':
                candidates += 1
                vaf = int(field[5].split(',')[1])/float(field[3])
                vafs.append(vaf)
                self.filteredGermMuts.append('\t'.join(line))

        if Options.verbose:
            for item in checks:
                if checks[item] > 0:
                    logging.info('%s: %s (%.1f%%)'%(theFlags[item], checks[item], checks[item]/float(len(self.rawGermMuts))*100))

            logging.info('Mean VAF: %.2f ;  Min VAF: %.2f ; Max VAF: %.2f'%( sum(vafs)/float(len(vafs)) , min(vafs), max(vafs)))

            logging.info("Total candidate germline mutations: %s"%(candidates))

    def _filterRawSomaticSNVs(self, Options, theFlags):
        vafs = []
        candidates = 0
        checks = {'^': 0, '`': 0, '&': 0, '!': 0, '*':0, '@':0, ';':0}
        base = {'A':'AU','C':'CU','G':'GU','T':'TU'}
        self.rawSNVSflags = ['' for i in range(0, len(self.rawSNVS))]
        keys = self.rawSNVSHead[len(self.rawSNVSHead)-1].split('\t')
        for i, line in enumerate(self.rawSNVS):
            line = line.split('\t')

            lineDat = dict(zip(keys, line))
            alt = lineDat['ALT']
            ref = lineDat['REF']

            # Filter multiple events at same loci
            if len(lineDat['ALT'].split(',')) > Options.maxevents:
                self.rawSNVSflags[i] += '^'
                checks['^'] += 1
            formatKey = lineDat['FORMAT'].split(':')
            normVals = dict(zip(formatKey, lineDat['NORMAL'].split(':')))
            tumorVals = dict(zip(formatKey, lineDat['TUMOR'].split(':')))
            # Filter depth less than opt in normal
            if int(normVals['DP']) < Options.minnr:
                self.rawSNVSflags[i] += '*'
                checks['*'] += 1
            # Filter variant reads in normal check
            if sum([int(t) for t in normVals[base[alt]].split(',')]) > 0:
                self.rawSNVSflags[i] += '@'
                checks['@'] += 1

            # Filter Tumor Total Depth
            if int(tumorVals['DP']) < Options.varsamdepth:
                self.rawSNVSflags[i] += '`'
                checks['`'] += 1

            # Filter Variant Reads Depth
            if int(tumorVals[base[alt]].split(',')[1]) < Options.minvr:
                self.rawSNVSflags[i] += '&'
                checks['&'] += 1

            # check clustered event Options.clusterFlank
            if ClusterMutCount(int(lineDat['POS']), self.chromsPos[lineDat['#CHROM']], Options) > 0:
                self.rawSNVSflags[i] += ';'
                checks[';'] += 1

            # Filter for VAF check
            if float(tumorVals['DP']) > 0:
                if int(tumorVals[base[alt]].split(',')[0])/float(int(tumorVals['DP'])-int(tumorVals['FDP'])) < Options.minvaf:
                    self.rawSNVSflags[i] += '!'
                    checks['!'] += 1

                # Total candidates
                if self.rawSNVSflags[i] == '':
                    candidates += 1
                    vaf = int(tumorVals[base[alt]].split(',')[0])/float(int(tumorVals['DP'])-int(tumorVals['FDP']))

                    if vaf>1.0:
                        print(vaf)
                        print(tumorVals)
                        print(line)
                        print(self.rawSNVSflags[i])
                        logging.error("Appears to be a problem with FilterVariant.py...see source code or use a different method.")
                        sys.exit()

                    vafs.append(vaf)

        if Options.verbose:
            for item in checks:
                if checks[item] > 0:
                    logging.info(
                        '%s: %s (%.1f%%)' % (theFlags[item], checks[item], checks[item] / float(len(self.rawSNVS)) * 100))

            logging.info('Mean VAF: %.2f ;  Min VAF: %.2f ; Max VAF: %.2f' % (
            sum(vafs) / float(len(vafs)), min(vafs), max(vafs)))

            logging.info("Total candidate somatic SNVs: %s" % (candidates))

    def _filterRawIndelSNVs(self, Options, theFlags):
        vafs = []
        candidates = 0
        checks = {'^': 0, '`': 0, '&': 0, '!': 0, '*':0, '@':0, ';':0}
        base = {'A':'AU','C':'CU','G':'GU','T':'TU'}
        self.rawIndelsflags = ['' for i in range(0, len(self.rawIndels))]
        keys = self.rawIndelsHead[len(self.rawIndelsHead)-1].split('\t')
        for i, line in enumerate(self.rawIndels):
            if line.split('\t')[6] == "PASS":
                self.rawPassIndels.append(line)
                #TODO implement way of filtering Strelka2 Indels...

    def _constructjointrefdictGerm(self):
        '''
        Constructs a dictionary with the following information {mutpos: [index,flags]}
        :return: a nested dictionary
        '''
        theRef = {}
        for i,m in enumerate(self.filteredGermMuts):
            key = m.split('\t')[0] + ':' + m.split('\t')[1] + '_' + m.split('\t')[3] + '/' + m.split('\t')[4]
            theRef.update({key:[i,'']})
        return(theRef)

    def _constructjointrefdictSomatic(self):
        '''
        Constructs a dictionary with the following information {mutpos: [index,flags]}
        :return: a nested dictionary
        '''
        theRef = {}
        for i,m in enumerate(self.rawSNVS):
            key = m.split('\t')[0] + ':' + m.split('\t')[1] + '_' + m.split('\t')[3] + '/' + m.split('\t')[4]
            theRef.update({key:[i,self.rawSNVSflags[i]]})
        return(theRef)

class Strelka2Main:

    def __init__(self, Options, theFlags):
        self.samplemap = []
        self.germlineDirs = {}
        self.somaticDirs = {}
        self.theSamples = []
        self.allSampleClasses = []
        self.germlineFilters = {}

        if Options.usesamplemap == True:
            self._sampleMapping(Options)

        # Step 1: Find all of the samples and set the above three items
        self._getSamples(Options)

        # Step 2: Read in all of the information for the files and setup classes
        self._setupSample(Options, theFlags)

        # Step 3: Perform Further Germline filtering...
        self._germlineFilterPart2(Options, theFlags)

        # Step 4: Append flags with Germline Evidence Flag...
        self._AddGermLineFlag(Options, theFlags)

        # Step 5: Update flags based on joint filtering...
        if Options.usesamplemap:
            self._PerformJointFiltering(Options, theFlags)
        else:
            pass

        # Step 6: Build files with headers....
        self._Strelka2WriteVCFFile(Options)

        if Options.ann:
            logging.info("Annotating vcf files using snpEff...")
            for sam in self.allSampleClasses:
                self._strelka2_AnnotateVCFFiles(Options.inDir, sam.SamNam, Options)

        if Options.table:
            self._BuildStrelka2Table(Options)

    def _sampleMapping(self, Options):
        with open(Options.samplemap, 'r') as samMap:
            self.samplemap = [sam.replace('\n', '').split(',') for sam in samMap.readlines()]

    def _getSamples(self, Options):
        '''
        Sets up the class variables for which directories house which files
        :param Options: Execution options.
        :return: None.
        '''
        logging.info("Searching for files...")
        indirs = glob.glob(Options.inDir + '*')
        sams = []
        for item in indirs:
            if 'germline' in item or 'somatic' in item:
                f = os.path.basename(item).split('.')
                sams.append(f[0])

                if 'germline' in item:
                    self.germlineDirs.update({f[0]: item})
                elif 'somatic' in item:
                    self.somaticDirs.update({f[0]: item})
                else:
                    logging.error("Unable to determine what to do with file %s"%(item))
        self.theSamples = list(set(sams))

    def _setupSample(self, Options, theFlags):
        for sam in self.theSamples:
            self.allSampleClasses.append(Strelka2Sample(Options, theFlags, sam, self.germlineDirs[sam], self.somaticDirs[sam]))

    def _germlineFilterPart2(self, Options, theFlags):
        if Options.usesamplemap:
            logging.info('...')
            logging.info('...')
            logging.info("Performing joint filtering on germline variants by updating filter flags...")
            for sam in self.samplemap:
                if len(sam) > 1:
                    sampleIDs = []
                    for indsam in sam:
                        for i, sample in enumerate(self.allSampleClasses):
                            if indsam in sample.SamNam:
                                sampleIDs.append(i)

                    # The mutation must be present in all samples for a germline filter to be applied...
                    # ~~~~~ Get Candidate Sites ~~~~#
                    logging.info("Gathering shared germline variant loci within specimen...")
                    # Step 1: loop over each sample set of positions w/ mutations and find the matching ones
                    matchingPositions = {}
                    for theid in sampleIDs:
                        for pos in self.allSampleClasses[theid].germjointsampleref:
                            try:
                                matchingPositions[pos] += 1
                            except KeyError:
                                matchingPositions.update({pos: 1})
                    # Step 2: Trim to shared positions w/ mutations to examine
                    candidateSites = []
                    for pos in matchingPositions:
                        if matchingPositions[pos] == len(sam):
                            candidateSites.append(pos)

                    logging.info("A total of %s candidate positions for joint filtering on %s..."%(len(candidateSites),','.join(sam)))

                    finalGermlinePos = []
                    for pos in candidateSites:
                        flagsets = []
                        for theid in sampleIDs:
                            try:
                                # Must do it this way to avoid errors later on. This means that >2 samples in a specimen and mutation shared between < total samples
                                a = theid
                                b = self.allSampleClasses[theid].SamNam
                                c = self.allSampleClasses[theid].germjointsampleref[pos]
                                d = self.allSampleClasses[theid].filteredGermMuts[self.allSampleClasses[theid].germjointsampleref[pos][0]]
                                flagsets.append(a)  # id of the sample inside of self.allMutectSamples
                                flagsets.append(b)  # Get Identifier
                                flagsets.append(c)  # Get Index and Flagset
                                flagsets.append(d)  # Get full mutation information
                            except KeyError:
                                pass

                        gpos = []
                        for z in range(0,len(flagsets),4):
                            indpos = '\t'.join(flagsets[z+3].split('\t')[0:5])

                            gpos.append(indpos)

                        gpos = list(set(gpos))
                        if len(gpos)==1:
                            finalGermlinePos.append(gpos[0])
                        else:
                            logging.error("GERMLINE MUTATIONS NOT UNIQUE. DO NOT TRUST RESULTS...")
                            sys.exit()
                    for theSam in sam:
                        self.germlineFilters.update({theSam:finalGermlinePos}) # Holds the final germline filters to be applied...
                else: # No paired for joint filtering
                    # TODO need to look at this section in greater detail when mixed matched and unmatched samples
                    finalGermlinePos = []
                    for i, sample in enumerate(self.allSampleClasses):
                        if sam[0] in sample.SamNam:
                            for m in self.allSampleClasses[i].filteredGermMuts:
                                finalGermlinePos.append('\t'.join(m.split('\t')[0:5]))
                            self.germlineFilters.update({sam[0]: finalGermlinePos})

        else:
            # Use all germline filters without flags, no need to further filter these.
            logging.info('...')
            logging.info('...')
            logging.info("Preparing germline variant extraction for filtering...")
            for i in range(0,len(self.allSampleClasses)):
                finalGermlinePos = []
                for m in self.allSampleClasses[i].filteredGermMuts:
                    finalGermlinePos.append('\t'.join(m.split('\t')[0:5]))
                self.germlineFilters.update({self.allSampleClasses[i].filteredGermMuts.SamNam:finalGermlinePos})

    def _AddGermLineFlag(self, Options, theFlags):
        # Updates the flags with variant called as germline flag
        logging.info("Appending variant flags with germline flag...")
        for i in range(0,len(self.allSampleClasses)):
            count = 0
            for z,m in enumerate(self.allSampleClasses[i].rawSNVS):
                cond = '\t'.join(m.split('\t')[0:5])
                if cond in self.germlineFilters[self.allSampleClasses[i].SamNam]:
                    count += 1
                    self.allSampleClasses[i].rawSNVSflags[z]+="G"
            logging.info("%s variants with germline evidence in %s..."%(count, self.allSampleClasses[i].SamNam))

    def _PerformJointFiltering(self, Options, theFlags):
        base = {'A': 'AU', 'C': 'CU', 'G': 'GU', 'T': 'TU'}
        logging.info('...')
        logging.info('...')
        logging.info("Performing joint filtering by updating filter flags...")
        count = 0
        # Responsible for modifying the flags for each of the variant filtered sets for the final output.
        for sam in self.samplemap:
            if len(sam) > 1:
                sampleIDs = []
                for indsam in sam:
                    for i, sample in enumerate(self.allSampleClasses):
                        if indsam in sample.SamNam:
                            sampleIDs.append(i)

                # Perform joint filter flag adjustments
                # ~~~~~ Get Candidate Sites ~~~~#
                # logging.info("Gathering shared mutation loci within specimen...")
                # Step 1: loop over each sample set of positions w/ mutations and find the matching ones
                matchingPositions = {}
                for theid in sampleIDs:
                    for pos in self.allSampleClasses[theid].snvjointsampleref:
                        try:
                            matchingPositions[pos] += 1
                        except KeyError:
                            matchingPositions.update({pos: 1})
                # Step 2: Trim to shared positions w/ muations to examine
                candidateSites = []
                for pos in matchingPositions:
                    if matchingPositions[pos] > 1:
                        candidateSites.append(pos)

                for pos in candidateSites:
                    flagsets = []
                    for theid in sampleIDs:
                        try:
                            # Must do it this way to avoid errors later on. This means that >2 samples in a specimen and mutation shared between < total samples
                            a = theid
                            b = self.allSampleClasses[theid].SamNam
                            c = self.allSampleClasses[theid].snvjointsampleref[pos]
                            d = self.allSampleClasses[theid].rawSNVS[self.allSampleClasses[theid].snvjointsampleref[pos][0]]
                            e = self.allSampleClasses[theid].rawSNVSHead[len(self.allSampleClasses[theid].rawSNVSHead)-1].split('\t')
                            flagsets.append(a)  # id of the sample inside of self.allMutectSamples
                            flagsets.append(b)  # Get Identifier
                            flagsets.append(c)  # Get Index and Flagset
                            flagsets.append(d)  # Get full mutation information
                            flagsets.append(e)  # Header information
                        except KeyError:
                            pass

                    # Check if a flagset needs to be jointly filtered
                    checkcount = 0
                    flagspresent = []
                    for z in range(0, len(flagsets), 5):
                        if flagsets[z + 2][1] != '':
                            checkcount += 1
                            flagspresent.append(flagsets[z + 2][1])

                    if checkcount > 0:
                        # Check if all flags are the same within the first (cond1 and cond2) or truly different...then adjust flags
                        if (len(flagspresent) == len(flagsets) / 4 and len(list(set(flagspresent))) != 1) or len(flagspresent) == 1:
                            if any(i in '.'.join(flagspresent) for i in ['!', '&']):
                                # ~~~~ Adjust flags ~~~~#
                                allVafs = []
                                allnvr = []

                                for z in range(0, len(flagsets), 5):
                                    m = flagsets[z + 3].split('\t')
                                    m = dict(zip( flagsets[z+4] , m ))

                                    alt = m['ALT']
                                    ref = m['REF']

                                    mvar = dict(zip(m['FORMAT'].split(':'), m['TUMOR'].split(':')))

                                    allnvr.append( int(mvar[base[alt]].split(',')[1]) )
                                    allVafs.append(int(mvar[base[alt]].split(',')[0])/float(int(mvar['DP'])-int(mvar['FDP'])))

                                # Step 1: Adjust for variant supporting reads flags (&)
                                if all(vr >= Options.uminvr for vr in allnvr) and any(vr >= Options.minvr for vr in allnvr):
                                    for z in range(0, len(flagsets), 5):
                                        self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]] = self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]].replace('&', '')
                                        count +=1
                                # Step 2: Adjust for variant allele frequencies flags (!)
                                if all(vaf >= Options.uminvaf for vaf in allVafs) and any(vaf >= Options.minvaf for vaf in allVafs):
                                    for z in range(0, len(flagsets), 5):
                                        self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]] = self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]].replace('!', '')
                                        count +=1
                            else:
                                pass  # Don't do any flag adjustments as none of the joint filtering options apply

        logging.info("A total of %s filter flags have been updated..."%(count))
        logging.info("Joint filtering completed...")
        logging.info('...')
        logging.info('...')

    def _Strelka2WriteVCFFile(self, Options):
        fullFail = ['@', '*', '&', '!', '^', '`', ';']  # All except '?'
        for i in range(0,len(self.allSampleClasses)):
            filtered = 0
            accepted = 0
            germflag = 0
            lowevsAccepted = 0
            finalaccepted = 0
            finalMutsToWrite = []
            for k,m in enumerate(self.allSampleClasses[i].rawSNVS):
                flagIDs = self.allSampleClasses[i].rawSNVSflags[k]
                if any(f in flagIDs for f in fullFail) == False:
                    if 'G' in flagIDs:
                        germflag +=1
                    else:
                        accepted += 1
                        if 'LowEVS' in m.split('\t')[6]:
                            lowevsAccepted += 1
                            if Options.ignoreLowEVS:
                                finalMutsToWrite.append(m)
                        else:
                            finalMutsToWrite.append(m)
                else:
                    filtered +=1
            logging.info('...')
            logging.info("?:?%s total sites rejected by variant filter: %s" % (self.allSampleClasses[i].SamNam, filtered))
            logging.info("?:?%s total sites accepted by variant filter: %s" % (self.allSampleClasses[i].SamNam, accepted))
            logging.info("?:?%s total sites with germline evidence in accepted: %s" % (self.allSampleClasses[i].SamNam, germflag))
            logging.info("?:?%s total sites with LowEVS flag in accepted: %s"  % (self.allSampleClasses[i].SamNam, lowevsAccepted))
            logging.info("?:?%s total SNVs/Indels: %s" % (self.allSampleClasses[i].SamNam, len(finalMutsToWrite)))
            logging.info('...')

            outputName = Options.inDir + self.allSampleClasses[i].SamNam + '.variantfilter.vcf'
            with open(outputName, 'w') as outFile:
                outFile.write(self.allSampleClasses[i].SNVSHeadout)
                for m in finalMutsToWrite:
                    outFile.write(m + '\n')

    def _strelka2_AnnotateVCFFiles(self, inputpath, samname, Options):
        inputName = Options.inDir + samname + '.variantfilter.vcf'
        outputName = Options.inDir + samname + '.variantfilter.ann.vcf'
        cmd = 'java -Xmx8G -jar %s ann -noStats -canon %s %s > %s' % (
        Options.jarfile, Options.ref, inputName, outputName)
        os.system(cmd)
        if os.path.exists(outputName) and os.path.exists(inputName):
            os.remove(inputName)

    def _BuildStrelka2Table(self, Options):
        logging.info("Constructing digested table for samples...")
        logging.info("Extracting forward and reverse reads from bam files...")
        logging.info("Extracting reference base...")
        base = {'A':'AU','C':'CU','G':'GU','T':'TU'}
        if Options.ann==False:
            pass
        #     fullFail = ['@','*','&','!','^','`',';']  # All except '?'
        #
        #     allMuts = []
        #     for sam in self.allMutectSamples:
        #         with open(sam.samplePath + sam.sampleidentifier.replace("passed.vcf","passed.variantfilter.vcf"), 'r') as inFile:
        #             lines = inFile.readlines()
        #         for line in lines:
        #             if line.startswith('#')==False:
        #                 allMuts.append(line)
        #     trinucs = GetTrinucs(allMuts, Options)
        #
        #     outTableName = self.allMutectSamples[0].samplePath + 'Mutect2.passed.variantfilter.digest.txt'
        #     header = "sample\tchr\tpos\tref_nt\tref_nt_pred\tmut_nt\treads\treads_fw\treads_rv\tdepth\tvaf\ttrinuc\tgene_name\timpact\tprotein_change\tnt_change\tprotein_id\n"
        #     outLines = []
        #     for sam in self.allMutectSamples:
        #         sample = sam.sampleidentifier.replace('.passed.vcf','')
        #         for i, flagIDs in enumerate(sam.flags):
        #             if any(f in flagIDs for f in fullFail) == False:
        #                 m = sam.rawmuts[i].split('\t') # Get the mutation
        #                 # Check for multiple events
        #                 if len(m[4].split(',')) > 1:
        #                     for z in range(0,len(m[4].split(','))):
        #                         chr = m[0]
        #                         pos = m[1]
        #                         ref_nt = m[3]
        #                         ref_nt_pred = m[3]
        #                         mut_nt = m[4].split(',')[z]
        #                         mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
        #                         reads = int(mvar['AD'].split(',')[z+1])
        #                         reads_fw = int(mvar['F1R2'].split(',')[z+1])
        #                         reads_rv = int(mvar['F2R1'].split(',')[z+1])
        #                         depth = sum([int(reads) for reads in mvar['AD'].split(',')])
        #                         vaf = reads / float(depth)
        #                         if len(mut_nt) == 1 and len(ref_nt) == 1:
        #                             try:
        #                                 trinuc = trinucs[str(chr) + ':' + str(pos)]
        #                             except KeyError:
        #                                 trinuc = '.'
        #                                 logging.warning("Unable to locate trinucleotide context...")
        #                         else:
        #                             trinuc = '.'
        #                         gene_name = '.'
        #                         impact = '.'
        #                         protein_change = '.'
        #                         nt_change = '.'
        #                         protein_id = '.'
        #                         l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
        #                              vaf, trinuc, gene_name, impact, protein_change, nt_change, protein_id]
        #                         strLine = [str(item) for item in l]
        #                         theLine = '\t'.join(strLine)
        #                         outLines.append(theLine)
        #                 else:
        #                     chr = m[0]
        #                     pos = m[1]
        #                     ref_nt = m[3]
        #                     ref_nt_pred = m[3]
        #                     mut_nt = m[4]
        #                     mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
        #                     reads = int(mvar['AD'].split(',')[1])
        #                     reads_fw = int(mvar['F1R2'].split(',')[1])
        #                     reads_rv = int(mvar['F2R1'].split(',')[1])
        #                     depth = sum([int(reads) for reads in mvar['AD'].split(',')])
        #                     vaf = reads/float(depth)
        #                     if len(mut_nt)==1 and len(ref_nt)==1:
        #                         try:
        #                             trinuc = trinucs[str(chr)+':'+str(pos)]
        #                         except KeyError:
        #                             trinuc = '.'
        #                             logging.warning("Unable to locate trinucleotide context...")
        #                     else:
        #                         trinuc = '.'
        #                     gene_name = '.'
        #                     impact = '.'
        #                     protein_change = '.'
        #                     nt_change = '.'
        #                     protein_id = '.'
        #                     l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
        #                          vaf, trinuc, gene_name, impact, protein_change, nt_change, protein_id]
        #                     strLine = [str(item) for item in l]
        #                     theLine = '\t'.join(strLine)
        #                     outLines.append(theLine)
        #     with open(outTableName, 'w') as outFile:
        #         outFile.write(header)
        #         outFile.write('\n'.join(outLines) + '\n')
        else: # If annotations are done...
            # Get trinucleotide context information
            allMuts = []
            # mutsBySam = dict.fromkeys([sam.SamNam for sam in self.allSampleClasses]) # Dictionary for each sample holding all the samples muts.
            for sam in self.allSampleClasses:
                with open(Options.inDir + sam.SamNam + ".variantfilter.ann.vcf", 'r') as inFile:
                    lines = inFile.readlines()
                for line in lines:
                    if line.startswith('#') == False:
                        allMuts.append(line)
            trinucs = GetTrinucs(allMuts, Options) # Gets trinucleotides for every mutation...

            annMuts = {}
            annMutReads = {} # Contains the forward and reverse reads
            for sam in self.allSampleClasses:
                with open(Options.inDir + sam.SamNam + ".variantfilter.ann.vcf", 'r') as inFile:
                    lines = [line.replace('\n','') for line in inFile.readlines()]
                allLines = []
                for line in lines:
                    if line.startswith('#') == False:
                        allLines.append(line)
                annMuts.update({sam.SamNam:allLines})
                readCounts = GetForwardAndReverseReads(sam.bamFile, allLines, mapq=10) # Get the forward and reverse reads
                annMutReads.update({sam.SamNam:readCounts})

            outTableName = Options.inDir + Options.inDir.split('/')[len(Options.inDir.split('/'))-2] + '.strelka2.passed.variantfilter.digest.txt'
            header = "sample\tchr\tpos\tref_nt\tref_nt_pred\tmut_nt\treads\treads_fw\treads_rv\tdepth\tvaf\ttrinuc\tgene_name\timpact\tprotein_change\tnt_change\ttranscript_id\tgene_id\n"
            outLines = []
            for sam in self.allSampleClasses:
                sample = sam.SamNam
                muts = annMuts[sam.SamNam]
                keys = sam.rawSNVSHead[len(sam.rawSNVSHead) - 1].split('\t')

                for headline in sam.rawSNVSHead:
                    if headline.startswith("##reference="):
                        ref_fasta = headline.split("##reference=file://")[1]

                for mut in muts:
                    m = mut.split('\t')  # Get the mutation
                    m = dict(zip(keys, m))

                    info = {}
                    for item in m['INFO'].split(';'):
                        infopair = item.split('=')
                        if len(infopair) != 2:
                            info.update({infopair[0]:'.'})
                        else:
                            info.update({infopair[0]:infopair[1]})
                    info = info['ANN'].split(',')[0]
                    # Check for multiple events
                    if len(m['ALT'].split(',')) > 1:
                        pass
                        #TODO implement this functionality for multiple events...
            #             for z in range(0, len(m[4].split(','))):
            #                 chr = m[0]
            #                 pos = m[1]
            #                 ref_nt = m[3]
            #                 ref_nt_pred = m[3]
            #                 mut_nt = m[4].split(',')[z]
            #                 mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
            #                 reads = int(mvar['AD'].split(',')[z + 1])
            #                 reads_fw = int(mvar['F1R2'].split(',')[z + 1])
            #                 reads_rv = int(mvar['F2R1'].split(',')[z + 1])
            #                 depth = sum([int(reads) for reads in mvar['AD'].split(',')])
            #                 vaf = reads / float(depth)
            #                 if len(mut_nt) == 1 and len(ref_nt) == 1:
            #                     try:
            #                         trinuc = trinucs[str(chr) + ':' + str(pos)]
            #                     except KeyError:
            #                         trinuc = '.'
            #                         logging.warning("Unable to locate trinucleotide context...")
            #                 else:
            #                     trinuc = '.'
            #                 gene_name = info.split('|')[3]
            #                 impact = info.split('|')[1].split('&')[0]
            #                 protein_change = info.split('|')[10]
            #                 nt_change = info.split('|')[9]
            #                 transcript_id = info.split('|')[6]
            #                 gene_id = info.split('|')[4]
            #                 l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
            #                      vaf, trinuc, gene_name, impact, protein_change, nt_change, transcript_id,gene_id]
            #                 strLine = [str(item) for item in l]
            #                 theLine = '\t'.join(strLine)
            #                 outLines.append(theLine)
                    else:
                        chr = m['#CHROM']
                        pos = m['POS']
                        ref_nt = GetReferenceBase(ref_fasta,chr,pos,int(pos)+len(m['REF'])-1)
                        ref_nt_pred = m['REF']
                        mut_nt = m['ALT']
                        formatKey = m['FORMAT'].split(':')
                        tumorVals = dict(zip(formatKey, m['TUMOR'].split(':')))
                        reads_fw = annMutReads[sam.SamNam]['\t'.join(mut.split('\t')[0:2])][0]
                        reads_rv = annMutReads[sam.SamNam]['\t'.join(mut.split('\t')[0:2])][1]
                        reads = int(tumorVals[base[mut_nt]].split(',')[1])
                        depth = int(tumorVals['DP'])-int(tumorVals['FDP'])
                        vaf = int(tumorVals[base[mut_nt]].split(',')[0])/float(depth)
                        if len(mut_nt) == 1 and len(ref_nt) == 1:
                            try:
                                trinuc = trinucs[str(chr) + ':' + str(pos)]
                            except KeyError:
                                trinuc = '.'
                                logging.warning("Unable to locate trinucleotide context...")
                        else:
                            trinuc = '.'
                        gene_name = info.split('|')[3]
                        impact = info.split('|')[1].split('&')[0]
                        protein_change = info.split('|')[10]
                        nt_change = info.split('|')[9]
                        transcript_id = info.split('|')[6]
                        gene_id = info.split('|')[4]
                        l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
                             vaf, trinuc, gene_name, impact, protein_change, nt_change, transcript_id, gene_id]
                        strLine = [str(item) for item in l]
                        theLine = '\t'.join(strLine)
                        outLines.append(theLine)
                with open(outTableName, 'w') as outFile:
                    outFile.write(header)
                    outFile.write('\n'.join(outLines) + '\n')
        logging.info("Table Completed...")

class DragenSample:

    def __init__(self, Options, theFlags, SamNam):
        self.SamNam = SamNam
        self.bamFile = None
        self.rawSNVS = []
        self.rawSNVSHead = []
        self.SNVSHeadout = ''
        self.rawSNVSflags = []
        self.rawIndels = []
        self.rawIndelsHead = []
        self.rawIndelsflags = []
        self.rawPassIndels = []

        # Step 1: Read in all corresponding files for GermMuts and Somatic SNVs
        self._readFile()
        self._newHeader(Options)

        # Step 1a: Get positions for cluster events
        self.chromsPos = PrepareClusteredEventCheck(self.rawSNVS)

        # Step 2: Perform Filtering
        logging.info("Filtering calls for %s..." % (self.SamNam))
        logging.info("SOMATIC filtering...")
        self._filterRawSomaticSNVs(Options, theFlags)
        # self._filterRawIndelSNVs(Options, theFlags)

        # if Options.usesamplemap:
        #     self.germjointsampleref = self._constructjointrefdictGerm()
        #     self.snvjointsampleref = self._constructjointrefdictSomatic()

    def _readFile(self):
        somSNVS = self.SamNam

        logging.info("Reading somatic SNVs...")
        with gzip.open(somSNVS, 'rb') as inSom:
            SomLines = [line.decode('UTF-8').replace('\n', '') for line in inSom.readlines()]
        for line in SomLines:
            if line.startswith('#'):
                self.rawSNVSHead.append(line)
            else:
                self.rawSNVS.append(line)

    def _newHeader(self, Options):
        argparsevars = {'minnr': '-nr', 'minvr': '-vr', 'minvaf': '-af', 'clusterFlank': '-l', 'table': '-t',
                        'inDir': '-i', 'caller': '-c', 'usesamplemap': '--usesamplemap', 'samplemap': '--samplemap',
                        'uminvr': '-uvr', 'uminvaf': '-uaf', 'verbose': '-v', 'maxevents': '-e', 'varsamdepth': '-vd',
                        'ann': '-a', 'ref': '-ref', 'jarfile': '-jar', 'refgenome': '-reffasta',
                        'ignoreLowEVS': '-ievs'}
        commandLine = 'FilterVariants.py ' + ' '.join(
            [' '.join([argparsevars[arg], str(getattr(Options, arg))]) for arg in vars(Options)])
        outline = '##FilterVariants=<ID=FilterVariants,CommandLine="%s">' % (commandLine)
        newHeader = ''
        for line in self.rawSNVSHead:
            if line.startswith('#'):
                if line.startswith('##reference='):
                    newHeader+=outline+'\n'
                    newHeader += line + '\n'
                else:
                    newHeader+=line+'\n'
        self.SNVSHeadout = newHeader

    def _filterRawSomaticSNVs(self, Options, theFlags):
        vafs = []
        candidates = 0
        checks = {'^': 0, '`': 0, '&': 0, '!': 0, '*':0, '@':0, ';':0}
        base = {'A':'AU','C':'CU','G':'GU','T':'TU'}
        self.rawSNVSflags = ['' for i in range(0, len(self.rawSNVS))]
        keys = self.rawSNVSHead[len(self.rawSNVSHead)-1].split('\t')
        for i, line in enumerate(self.rawSNVS):
            line = line.split('\t')

            lineDat = dict(zip(keys, line))
            alt = lineDat['ALT']
            ref = lineDat['REF']

            # Have to add the tumor normal pairs as these are not present in these vcfs in an easy to handle way.
            # Get Column 1 Normal
            # Get Column 2 Tumor
            lineDat.update({'NORMAL': line[9]})
            lineDat.update({'TUMOR': line[10]})

            # Filter multiple events at same loci
            if len(lineDat['ALT'].split(',')) > Options.maxevents:
                self.rawSNVSflags[i] += '^'
                checks['^'] += 1
            formatKey = lineDat['FORMAT'].split(':')
            normVals = dict(zip(formatKey, lineDat['NORMAL'].split(':')))
            tumorVals = dict(zip(formatKey, lineDat['TUMOR'].split(':')))

            # Filter depth less than opt in normal
            if int(normVals['DP']) < Options.minnr:
                self.rawSNVSflags[i] += '*'
                checks['*'] += 1
            # Filter variant reads in normal check
            if int(normVals['AD'].split(',')[1]) > 0:
                self.rawSNVSflags[i] += '@'
                checks['@'] += 1

            # Filter Tumor Total Depth
            if int(tumorVals['DP']) < Options.varsamdepth:
                self.rawSNVSflags[i] += '`'
                checks['`'] += 1

            # Filter Variant Reads Depth
            if int(tumorVals['AD'].split(',')[1]) < Options.minvr:
                self.rawSNVSflags[i] += '&'
                checks['&'] += 1

            # check clustered event Options.clusterFlank
            if ClusterMutCount(int(lineDat['POS']), self.chromsPos[lineDat['#CHROM']], Options) > 0:
                self.rawSNVSflags[i] += ';'
                checks[';'] += 1

            # Filter for VAF check
            if float(tumorVals['DP']) > 0:
                try:
                    vaf = int(tumorVals['AD'].split(',')[1])/float(int(tumorVals['AD'].split(',')[1])+int(tumorVals['AD'].split(',')[0]))
                except ZeroDivisionError:
                    vaf=0.0
                if vaf < Options.minvaf:
                    self.rawSNVSflags[i] += '!'
                    checks['!'] += 1

                # Total candidates
                if self.rawSNVSflags[i] == '':
                    candidates += 1

                    if vaf>1.0:
                        print(vaf)
                        print(tumorVals)
                        print(line)
                        print(self.rawSNVSflags[i])
                        logging.error("Appears to be a problem with FilterVariant.py...see source code or use a different method.")
                        sys.exit()

                    vafs.append(vaf)

        if Options.verbose:
            for item in checks:
                if checks[item] > 0:
                    logging.info(
                        '%s: %s (%.1f%%)' % (theFlags[item], checks[item], checks[item] / float(len(self.rawSNVS)) * 100))

            logging.info('Mean VAF: %.2f ;  Min VAF: %.2f ; Max VAF: %.2f' % (
            sum(vafs) / float(len(vafs)), min(vafs), max(vafs)))

            logging.info("Total candidate somatic SNVs: %s" % (candidates))

    def _filterRawIndelSNVs(self, Options, theFlags):
        vafs = []
        candidates = 0
        checks = {'^': 0, '`': 0, '&': 0, '!': 0, '*':0, '@':0, ';':0}
        base = {'A':'AU','C':'CU','G':'GU','T':'TU'}
        self.rawIndelsflags = ['' for i in range(0, len(self.rawIndels))]
        keys = self.rawIndelsHead[len(self.rawIndelsHead)-1].split('\t')
        for i, line in enumerate(self.rawIndels):
            if line.split('\t')[6] == "PASS":
                self.rawPassIndels.append(line)
                #TODO implement way of filtering Strelka2 Indels...

    def _constructjointrefdictGerm(self):
        '''
        Constructs a dictionary with the following information {mutpos: [index,flags]}
        :return: a nested dictionary
        '''
        theRef = {}
        for i,m in enumerate(self.filteredGermMuts):
            key = m.split('\t')[0] + ':' + m.split('\t')[1] + '_' + m.split('\t')[3] + '/' + m.split('\t')[4]
            theRef.update({key:[i,'']})
        return(theRef)

    def _constructjointrefdictSomatic(self):
        '''
        Constructs a dictionary with the following information {mutpos: [index,flags]}
        :return: a nested dictionary
        '''
        theRef = {}
        for i,m in enumerate(self.rawSNVS):
            key = m.split('\t')[0] + ':' + m.split('\t')[1] + '_' + m.split('\t')[3] + '/' + m.split('\t')[4]
            theRef.update({key:[i,self.rawSNVSflags[i]]})
        return(theRef)

class DragenMain:

    def __init__(self, Options, theFlags):
        self.samplemap = []
        self.germlineDirs = {}
        self.somaticDirs = {}
        self.theSamples = []
        self.allSampleClasses = []
        self.germlineFilters = {}

        if Options.usesamplemap == True:
            self._sampleMapping(Options)

        # Step 1: Find all of the samples and set the above three items
        self._getSamples(Options)

        # Step 2: Read in all of the information for the files and setup classes
        self._setupSample(Options, theFlags)

        # # Step 3: Perform Further Germline filtering...
        # self._germlineFilterPart2(Options, theFlags)
        #
        # # Step 4: Append flags with Germline Evidence Flag...
        # self._AddGermLineFlag(Options, theFlags)
        #
        # # Step 5: Update flags based on joint filtering...
        # if Options.usesamplemap:
        #     self._PerformJointFiltering(Options, theFlags)
        # else:
        #     pass

        # Step 6: Build files with headers....
        self._Strelka2WriteVCFFile(Options)

        # if Options.ann:
        #     logging.info("Annotating vcf files using snpEff...")
        #     for sam in self.allSampleClasses:
        #         self._strelka2_AnnotateVCFFiles(Options.inDir, sam.SamNam, Options)
        #
        # if Options.table:
        #     self._BuildStrelka2Table(Options)

    def _sampleMapping(self, Options):
        with open(Options.samplemap, 'r') as samMap:
            self.samplemap = [sam.replace('\n', '').split(',') for sam in samMap.readlines()]

    def _getSamples(self, Options):
        '''
        Sets up the class variables for which directories house which files
        :param Options: Execution options.
        :return: None.
        '''
        infiles = glob.glob(Options.inDir + '*.vcf.gz')
        if Options.verbose:
            logging.info("A total of %s vcf files have been found..." % (len(infiles)))
        logging.info("...")
        logging.info("...")

        self.theSamples = list(set([os.path.abspath(thefile) for thefile in infiles]))

    def _setupSample(self, Options, theFlags):
        for sam in self.theSamples:
            samClass = DragenSample(Options, theFlags, sam)
            self._Strelka2WriteVCFFile(Options, samClass)
            samClass = None

    def _germlineFilterPart2(self, Options, theFlags):
        if Options.usesamplemap:
            logging.info('...')
            logging.info('...')
            logging.info("Performing joint filtering on germline variants by updating filter flags...")
            for sam in self.samplemap:
                if len(sam) > 1:
                    sampleIDs = []
                    for indsam in sam:
                        for i, sample in enumerate(self.allSampleClasses):
                            if indsam in sample.SamNam:
                                sampleIDs.append(i)

                    # The mutation must be present in all samples for a germline filter to be applied...
                    # ~~~~~ Get Candidate Sites ~~~~#
                    logging.info("Gathering shared germline variant loci within specimen...")
                    # Step 1: loop over each sample set of positions w/ mutations and find the matching ones
                    matchingPositions = {}
                    for theid in sampleIDs:
                        for pos in self.allSampleClasses[theid].germjointsampleref:
                            try:
                                matchingPositions[pos] += 1
                            except KeyError:
                                matchingPositions.update({pos: 1})
                    # Step 2: Trim to shared positions w/ mutations to examine
                    candidateSites = []
                    for pos in matchingPositions:
                        if matchingPositions[pos] == len(sam):
                            candidateSites.append(pos)

                    logging.info("A total of %s candidate positions for joint filtering on %s..."%(len(candidateSites),','.join(sam)))

                    finalGermlinePos = []
                    for pos in candidateSites:
                        flagsets = []
                        for theid in sampleIDs:
                            try:
                                # Must do it this way to avoid errors later on. This means that >2 samples in a specimen and mutation shared between < total samples
                                a = theid
                                b = self.allSampleClasses[theid].SamNam
                                c = self.allSampleClasses[theid].germjointsampleref[pos]
                                d = self.allSampleClasses[theid].filteredGermMuts[self.allSampleClasses[theid].germjointsampleref[pos][0]]
                                flagsets.append(a)  # id of the sample inside of self.allMutectSamples
                                flagsets.append(b)  # Get Identifier
                                flagsets.append(c)  # Get Index and Flagset
                                flagsets.append(d)  # Get full mutation information
                            except KeyError:
                                pass

                        gpos = []
                        for z in range(0,len(flagsets),4):
                            indpos = '\t'.join(flagsets[z+3].split('\t')[0:5])

                            gpos.append(indpos)

                        gpos = list(set(gpos))
                        if len(gpos)==1:
                            finalGermlinePos.append(gpos[0])
                        else:
                            logging.error("GERMLINE MUTATIONS NOT UNIQUE. DO NOT TRUST RESULTS...")
                            sys.exit()
                    for theSam in sam:
                        self.germlineFilters.update({theSam:finalGermlinePos}) # Holds the final germline filters to be applied...
                else: # No paired for joint filtering
                    # TODO need to look at this section in greater detail when mixed matched and unmatched samples
                    finalGermlinePos = []
                    for i, sample in enumerate(self.allSampleClasses):
                        if sam[0] in sample.SamNam:
                            for m in self.allSampleClasses[i].filteredGermMuts:
                                finalGermlinePos.append('\t'.join(m.split('\t')[0:5]))
                            self.germlineFilters.update({sam[0]: finalGermlinePos})

        else:
            # Use all germline filters without flags, no need to further filter these.
            logging.info('...')
            logging.info('...')
            logging.info("Preparing germline variant extraction for filtering...")
            for i in range(0,len(self.allSampleClasses)):
                finalGermlinePos = []
                for m in self.allSampleClasses[i].filteredGermMuts:
                    finalGermlinePos.append('\t'.join(m.split('\t')[0:5]))
                self.germlineFilters.update({self.allSampleClasses[i].filteredGermMuts.SamNam:finalGermlinePos})

    def _AddGermLineFlag(self, Options, theFlags):
        # Updates the flags with variant called as germline flag
        logging.info("Appending variant flags with germline flag...")
        for i in range(0,len(self.allSampleClasses)):
            count = 0
            for z,m in enumerate(self.allSampleClasses[i].rawSNVS):
                cond = '\t'.join(m.split('\t')[0:5])
                if cond in self.germlineFilters[self.allSampleClasses[i].SamNam]:
                    count += 1
                    self.allSampleClasses[i].rawSNVSflags[z]+="G"
            logging.info("%s variants with germline evidence in %s..."%(count, self.allSampleClasses[i].SamNam))

    def _PerformJointFiltering(self, Options, theFlags):
        base = {'A': 'AU', 'C': 'CU', 'G': 'GU', 'T': 'TU'}
        logging.info('...')
        logging.info('...')
        logging.info("Performing joint filtering by updating filter flags...")
        count = 0
        # Responsible for modifying the flags for each of the variant filtered sets for the final output.
        for sam in self.samplemap:
            if len(sam) > 1:
                sampleIDs = []
                for indsam in sam:
                    for i, sample in enumerate(self.allSampleClasses):
                        if indsam in sample.SamNam:
                            sampleIDs.append(i)

                # Perform joint filter flag adjustments
                # ~~~~~ Get Candidate Sites ~~~~#
                # logging.info("Gathering shared mutation loci within specimen...")
                # Step 1: loop over each sample set of positions w/ mutations and find the matching ones
                matchingPositions = {}
                for theid in sampleIDs:
                    for pos in self.allSampleClasses[theid].snvjointsampleref:
                        try:
                            matchingPositions[pos] += 1
                        except KeyError:
                            matchingPositions.update({pos: 1})
                # Step 2: Trim to shared positions w/ muations to examine
                candidateSites = []
                for pos in matchingPositions:
                    if matchingPositions[pos] > 1:
                        candidateSites.append(pos)

                for pos in candidateSites:
                    flagsets = []
                    for theid in sampleIDs:
                        try:
                            # Must do it this way to avoid errors later on. This means that >2 samples in a specimen and mutation shared between < total samples
                            a = theid
                            b = self.allSampleClasses[theid].SamNam
                            c = self.allSampleClasses[theid].snvjointsampleref[pos]
                            d = self.allSampleClasses[theid].rawSNVS[self.allSampleClasses[theid].snvjointsampleref[pos][0]]
                            e = self.allSampleClasses[theid].rawSNVSHead[len(self.allSampleClasses[theid].rawSNVSHead)-1].split('\t')
                            flagsets.append(a)  # id of the sample inside of self.allMutectSamples
                            flagsets.append(b)  # Get Identifier
                            flagsets.append(c)  # Get Index and Flagset
                            flagsets.append(d)  # Get full mutation information
                            flagsets.append(e)  # Header information
                        except KeyError:
                            pass

                    # Check if a flagset needs to be jointly filtered
                    checkcount = 0
                    flagspresent = []
                    for z in range(0, len(flagsets), 5):
                        if flagsets[z + 2][1] != '':
                            checkcount += 1
                            flagspresent.append(flagsets[z + 2][1])

                    if checkcount > 0:
                        # Check if all flags are the same within the first (cond1 and cond2) or truly different...then adjust flags
                        if (len(flagspresent) == len(flagsets) / 4 and len(list(set(flagspresent))) != 1) or len(flagspresent) == 1:
                            if any(i in '.'.join(flagspresent) for i in ['!', '&']):
                                # ~~~~ Adjust flags ~~~~#
                                allVafs = []
                                allnvr = []

                                for z in range(0, len(flagsets), 5):
                                    m = flagsets[z + 3].split('\t')
                                    m = dict(zip( flagsets[z+4] , m ))

                                    alt = m['ALT']
                                    ref = m['REF']

                                    mvar = dict(zip(m['FORMAT'].split(':'), m['TUMOR'].split(':')))

                                    allnvr.append( int(mvar[base[alt]].split(',')[1]) )
                                    allVafs.append(int(mvar[base[alt]].split(',')[0])/float(int(mvar['DP'])-int(mvar['FDP'])))

                                # Step 1: Adjust for variant supporting reads flags (&)
                                if all(vr >= Options.uminvr for vr in allnvr) and any(vr >= Options.minvr for vr in allnvr):
                                    for z in range(0, len(flagsets), 5):
                                        self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]] = self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]].replace('&', '')
                                        count +=1
                                # Step 2: Adjust for variant allele frequencies flags (!)
                                if all(vaf >= Options.uminvaf for vaf in allVafs) and any(vaf >= Options.minvaf for vaf in allVafs):
                                    for z in range(0, len(flagsets), 5):
                                        self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]] = self.allSampleClasses[flagsets[z]].rawSNVSflags[flagsets[z + 2][0]].replace('!', '')
                                        count +=1
                            else:
                                pass  # Don't do any flag adjustments as none of the joint filtering options apply

        logging.info("A total of %s filter flags have been updated..."%(count))
        logging.info("Joint filtering completed...")
        logging.info('...')
        logging.info('...')

    def _Strelka2WriteVCFFile(self, Options, dragenSample):
        fullFail = ['@', '*', '&', '!', '^', '`', ';']  # All except '?'
        filtered = 0
        accepted = 0
        germflag = 0
        lowevsAccepted = 0
        finalaccepted = 0
        finalMutsToWrite = []
        for k,m in enumerate(dragenSample.rawSNVS):
            flagIDs = dragenSample.rawSNVSflags[k]
            if any(f in flagIDs for f in fullFail) == False:
                if 'G' in flagIDs:
                    germflag +=1
                else:
                    accepted += 1
                    if 'LowEVS' in m.split('\t')[6]:
                        lowevsAccepted += 1
                        if Options.ignoreLowEVS:
                            finalMutsToWrite.append(m)
                    else:
                        finalMutsToWrite.append(m)
            else:
                filtered +=1
        logging.info('...')
        logging.info("?:?%s total sites rejected by variant filter: %s" % (dragenSample.SamNam, filtered))
        logging.info("?:?%s total sites accepted by variant filter: %s" % (dragenSample.SamNam, accepted))
        logging.info("?:?%s total sites with germline evidence in accepted: %s" % (dragenSample.SamNam, germflag))
        logging.info("?:?%s total sites with LowEVS flag in accepted: %s"  % (dragenSample.SamNam, lowevsAccepted))
        logging.info("?:?%s total SNVs/Indels: %s" % (dragenSample.SamNam, len(finalMutsToWrite)))
        logging.info('...')

        outputName = dragenSample.SamNam.replace('.vcf.gz','.variantfilter.vcf')
        with open(outputName, 'w') as outFile:
            outFile.write(dragenSample.SNVSHeadout)
            for m in finalMutsToWrite:
                outFile.write(m + '\n')

    def _strelka2_AnnotateVCFFiles(self, inputpath, samname, Options):
        inputName = Options.inDir + samname + '.variantfilter.vcf'
        outputName = Options.inDir + samname + '.variantfilter.ann.vcf'
        cmd = 'java -Xmx8G -jar %s ann -noStats -canon %s %s > %s' % (
        Options.jarfile, Options.ref, inputName, outputName)
        os.system(cmd)
        if os.path.exists(outputName) and os.path.exists(inputName):
            os.remove(inputName)

    def _BuildStrelka2Table(self, Options):
        logging.info("Constructing digested table for samples...")
        logging.info("Extracting forward and reverse reads from bam files...")
        logging.info("Extracting reference base...")
        base = {'A':'AU','C':'CU','G':'GU','T':'TU'}
        if Options.ann==False:
            pass
        #     fullFail = ['@','*','&','!','^','`',';']  # All except '?'
        #
        #     allMuts = []
        #     for sam in self.allMutectSamples:
        #         with open(sam.samplePath + sam.sampleidentifier.replace("passed.vcf","passed.variantfilter.vcf"), 'r') as inFile:
        #             lines = inFile.readlines()
        #         for line in lines:
        #             if line.startswith('#')==False:
        #                 allMuts.append(line)
        #     trinucs = GetTrinucs(allMuts, Options)
        #
        #     outTableName = self.allMutectSamples[0].samplePath + 'Mutect2.passed.variantfilter.digest.txt'
        #     header = "sample\tchr\tpos\tref_nt\tref_nt_pred\tmut_nt\treads\treads_fw\treads_rv\tdepth\tvaf\ttrinuc\tgene_name\timpact\tprotein_change\tnt_change\tprotein_id\n"
        #     outLines = []
        #     for sam in self.allMutectSamples:
        #         sample = sam.sampleidentifier.replace('.passed.vcf','')
        #         for i, flagIDs in enumerate(sam.flags):
        #             if any(f in flagIDs for f in fullFail) == False:
        #                 m = sam.rawmuts[i].split('\t') # Get the mutation
        #                 # Check for multiple events
        #                 if len(m[4].split(',')) > 1:
        #                     for z in range(0,len(m[4].split(','))):
        #                         chr = m[0]
        #                         pos = m[1]
        #                         ref_nt = m[3]
        #                         ref_nt_pred = m[3]
        #                         mut_nt = m[4].split(',')[z]
        #                         mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
        #                         reads = int(mvar['AD'].split(',')[z+1])
        #                         reads_fw = int(mvar['F1R2'].split(',')[z+1])
        #                         reads_rv = int(mvar['F2R1'].split(',')[z+1])
        #                         depth = sum([int(reads) for reads in mvar['AD'].split(',')])
        #                         vaf = reads / float(depth)
        #                         if len(mut_nt) == 1 and len(ref_nt) == 1:
        #                             try:
        #                                 trinuc = trinucs[str(chr) + ':' + str(pos)]
        #                             except KeyError:
        #                                 trinuc = '.'
        #                                 logging.warning("Unable to locate trinucleotide context...")
        #                         else:
        #                             trinuc = '.'
        #                         gene_name = '.'
        #                         impact = '.'
        #                         protein_change = '.'
        #                         nt_change = '.'
        #                         protein_id = '.'
        #                         l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
        #                              vaf, trinuc, gene_name, impact, protein_change, nt_change, protein_id]
        #                         strLine = [str(item) for item in l]
        #                         theLine = '\t'.join(strLine)
        #                         outLines.append(theLine)
        #                 else:
        #                     chr = m[0]
        #                     pos = m[1]
        #                     ref_nt = m[3]
        #                     ref_nt_pred = m[3]
        #                     mut_nt = m[4]
        #                     mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
        #                     reads = int(mvar['AD'].split(',')[1])
        #                     reads_fw = int(mvar['F1R2'].split(',')[1])
        #                     reads_rv = int(mvar['F2R1'].split(',')[1])
        #                     depth = sum([int(reads) for reads in mvar['AD'].split(',')])
        #                     vaf = reads/float(depth)
        #                     if len(mut_nt)==1 and len(ref_nt)==1:
        #                         try:
        #                             trinuc = trinucs[str(chr)+':'+str(pos)]
        #                         except KeyError:
        #                             trinuc = '.'
        #                             logging.warning("Unable to locate trinucleotide context...")
        #                     else:
        #                         trinuc = '.'
        #                     gene_name = '.'
        #                     impact = '.'
        #                     protein_change = '.'
        #                     nt_change = '.'
        #                     protein_id = '.'
        #                     l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
        #                          vaf, trinuc, gene_name, impact, protein_change, nt_change, protein_id]
        #                     strLine = [str(item) for item in l]
        #                     theLine = '\t'.join(strLine)
        #                     outLines.append(theLine)
        #     with open(outTableName, 'w') as outFile:
        #         outFile.write(header)
        #         outFile.write('\n'.join(outLines) + '\n')
        else: # If annotations are done...
            # Get trinucleotide context information
            allMuts = []
            # mutsBySam = dict.fromkeys([sam.SamNam for sam in self.allSampleClasses]) # Dictionary for each sample holding all the samples muts.
            for sam in self.allSampleClasses:
                with open(Options.inDir + sam.SamNam + ".variantfilter.ann.vcf", 'r') as inFile:
                    lines = inFile.readlines()
                for line in lines:
                    if line.startswith('#') == False:
                        allMuts.append(line)
            trinucs = GetTrinucs(allMuts, Options) # Gets trinucleotides for every mutation...

            annMuts = {}
            annMutReads = {} # Contains the forward and reverse reads
            for sam in self.allSampleClasses:
                with open(Options.inDir + sam.SamNam + ".variantfilter.ann.vcf", 'r') as inFile:
                    lines = [line.replace('\n','') for line in inFile.readlines()]
                allLines = []
                for line in lines:
                    if line.startswith('#') == False:
                        allLines.append(line)
                annMuts.update({sam.SamNam:allLines})
                readCounts = GetForwardAndReverseReads(sam.bamFile, allLines, mapq=10) # Get the forward and reverse reads
                annMutReads.update({sam.SamNam:readCounts})

            outTableName = Options.inDir + Options.inDir.split('/')[len(Options.inDir.split('/'))-2] + '.strelka2.passed.variantfilter.digest.txt'
            header = "sample\tchr\tpos\tref_nt\tref_nt_pred\tmut_nt\treads\treads_fw\treads_rv\tdepth\tvaf\ttrinuc\tgene_name\timpact\tprotein_change\tnt_change\ttranscript_id\tgene_id\n"
            outLines = []
            for sam in self.allSampleClasses:
                sample = sam.SamNam
                muts = annMuts[sam.SamNam]
                keys = sam.rawSNVSHead[len(sam.rawSNVSHead) - 1].split('\t')

                for headline in sam.rawSNVSHead:
                    if headline.startswith("##reference="):
                        ref_fasta = headline.split("##reference=file://")[1]

                for mut in muts:
                    m = mut.split('\t')  # Get the mutation
                    m = dict(zip(keys, m))

                    info = {}
                    for item in m['INFO'].split(';'):
                        infopair = item.split('=')
                        if len(infopair) != 2:
                            info.update({infopair[0]:'.'})
                        else:
                            info.update({infopair[0]:infopair[1]})
                    info = info['ANN'].split(',')[0]
                    # Check for multiple events
                    if len(m['ALT'].split(',')) > 1:
                        pass
                        #TODO implement this functionality for multiple events...
            #             for z in range(0, len(m[4].split(','))):
            #                 chr = m[0]
            #                 pos = m[1]
            #                 ref_nt = m[3]
            #                 ref_nt_pred = m[3]
            #                 mut_nt = m[4].split(',')[z]
            #                 mvar = dict(zip(m[8].split(':'), m[sam.tumor_col].split(':')))
            #                 reads = int(mvar['AD'].split(',')[z + 1])
            #                 reads_fw = int(mvar['F1R2'].split(',')[z + 1])
            #                 reads_rv = int(mvar['F2R1'].split(',')[z + 1])
            #                 depth = sum([int(reads) for reads in mvar['AD'].split(',')])
            #                 vaf = reads / float(depth)
            #                 if len(mut_nt) == 1 and len(ref_nt) == 1:
            #                     try:
            #                         trinuc = trinucs[str(chr) + ':' + str(pos)]
            #                     except KeyError:
            #                         trinuc = '.'
            #                         logging.warning("Unable to locate trinucleotide context...")
            #                 else:
            #                     trinuc = '.'
            #                 gene_name = info.split('|')[3]
            #                 impact = info.split('|')[1].split('&')[0]
            #                 protein_change = info.split('|')[10]
            #                 nt_change = info.split('|')[9]
            #                 transcript_id = info.split('|')[6]
            #                 gene_id = info.split('|')[4]
            #                 l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
            #                      vaf, trinuc, gene_name, impact, protein_change, nt_change, transcript_id,gene_id]
            #                 strLine = [str(item) for item in l]
            #                 theLine = '\t'.join(strLine)
            #                 outLines.append(theLine)
                    else:
                        chr = m['#CHROM']
                        pos = m['POS']
                        ref_nt = GetReferenceBase(ref_fasta,chr,pos,int(pos)+len(m['REF'])-1)
                        ref_nt_pred = m['REF']
                        mut_nt = m['ALT']
                        formatKey = m['FORMAT'].split(':')
                        tumorVals = dict(zip(formatKey, m['TUMOR'].split(':')))
                        reads_fw = annMutReads[sam.SamNam]['\t'.join(mut.split('\t')[0:2])][0]
                        reads_rv = annMutReads[sam.SamNam]['\t'.join(mut.split('\t')[0:2])][1]
                        reads = int(tumorVals[base[mut_nt]].split(',')[1])
                        depth = int(tumorVals['DP'])-int(tumorVals['FDP'])
                        vaf = int(tumorVals[base[mut_nt]].split(',')[0])/float(depth)
                        if len(mut_nt) == 1 and len(ref_nt) == 1:
                            try:
                                trinuc = trinucs[str(chr) + ':' + str(pos)]
                            except KeyError:
                                trinuc = '.'
                                logging.warning("Unable to locate trinucleotide context...")
                        else:
                            trinuc = '.'
                        gene_name = info.split('|')[3]
                        impact = info.split('|')[1].split('&')[0]
                        protein_change = info.split('|')[10]
                        nt_change = info.split('|')[9]
                        transcript_id = info.split('|')[6]
                        gene_id = info.split('|')[4]
                        l = [sample, chr, pos, ref_nt, ref_nt_pred, mut_nt, reads, reads_fw, reads_rv, depth,
                             vaf, trinuc, gene_name, impact, protein_change, nt_change, transcript_id, gene_id]
                        strLine = [str(item) for item in l]
                        theLine = '\t'.join(strLine)
                        outLines.append(theLine)
                with open(outTableName, 'w') as outFile:
                    outFile.write(header)
                    outFile.write('\n'.join(outLines) + '\n')
        logging.info("Table Completed...")

def GetForwardAndReverseReads(bam, muts, mapq=30):
    '''
    Obtains the forward and reverse read counts covering a position of interest given a mapping quality of 30 or other provided by user. Now using pysam rather than subprocess.
    Old commands are:
    fwcmd = ['samtools','view','-F','16','--threads','8','-L','tmp.bedlike.bed',bam]
    rvcmd = ['samtools','view','-f','16','--threads','8','-L','tmp.bedlike.bed',bam]
    :param bam: Bam file to extract read information from
    :param muts: List of Muts
    :return: Dictionary of {position:[fwd,rv]}
    '''
    logging.info("Extracting read counts...")

    bamfile = pysam.AlignmentFile(bam, 'rb')

    reads = {}
    for mut in muts:
        pos = '\t'.join(mut.split('\t')[0:2]) + '\t' + str(int(mut.split('\t')[1])+len(mut.split('\t')[4]))
        chrom = pos.split('\t')[0]
        start = int(pos.split('\t')[1])
        end = int(pos.split('\t')[2])

        bamiter = bamfile.fetch(chrom, start, end)
        r1 = 0
        r2 = 0
        for seq in bamiter:
            if seq.mapping_quality > mapq:
                if seq.is_read1:
                    r1 += 1
                else:
                    r2 += 1
        reads.update({'\t'.join(mut.split('\t')[0:2]):[r1,r2]})


    return(reads)

def GetReferenceBase(RefFasta, chrom, start, end):
    '''
    Obtains the base(s) for positions from the reference sequence for a variant within a VCF file.

    :param RefFasta: Reference fasta used for
    :param chrom: Chromosome in question
    :param start: Pos in vcf
    :param end: Pos plus length of the reference - 1
    :return: The base(s) from the reference fasta
    '''
    cmd = ["samtools","faidx",RefFasta,'%s:%s-%s'%(chrom,start,end)]
    val = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    val.wait()
    base = subprocess.Popen([ "grep","-v",'">"'], stdin=val.stdout, stdout=subprocess.PIPE)
    base.wait()
    theBase = base.communicate()[0].decode(encoding="utf-8").split('\n')[1]
    return(theBase)

def GetTrinucs(allMuts, Options):
    '''
    Takes in all of the passed and filtered variants and the program options and returns a dictionary of positions with
    the trinucleotide context. Note: requires bedtools

    :param allMuts: Passed and filtered variants.
    :param Options: Program options
    :return: Dictionary (Genomic position : trinucleotide context)
    '''
    # Step 1: Prepare unique bed
    # allMuts = [':'.join(m) for m in allMuts]
    logging.info("Creating temporary bed file...")
    with open('__tmpbed.bed','w') as tmpbed:
        for m in allMuts:
            if ',' not in m.split('\t')[4] and len(m.split('\t')[4]) == 1:
                chr, pos = m.split('\t')[0], int(m.split('\t')[1])
                bedline = '%s\t%s\t%s\n'%(chr, pos-2, pos+1)
                tmpbed.write(bedline)

    # Step 2: Sort bed and make uniqe
    cmd = ['bedtools','sort','-i','__tmpbed.bed']
    sort = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    uniq = subprocess.Popen(['uniq'], stdin=sort.stdout, stdout=subprocess.PIPE)
    with open('__tmpbed.sorted.uniquified.bed','wb') as tmpbed:
        tmpbed.write(uniq.communicate()[0])

    logging.info("Extracting regions for trinucleotide contexts...")
    cmd = "bedtools getfasta -tab -fi %s -bed __tmpbed.sorted.uniquified.bed -fo __tri.tsv"%(Options.refgenome)
    extract = subprocess.Popen(cmd.split(' '))
    extract.wait()

    with open('__tri.tsv', 'r') as trifile:
        lines = [t.replace('\n','').split('\t') for t in trifile.readlines()]
    tridict = {}
    for line in lines:
        key = line[0].split(':')[0] + ':' + str(int(line[0].split(':')[1].split('-')[0])+2)
        val = line[1]
        tridict.update({key:val})

    # Clean up temporary files
    os.system('rm -f __tmpbed.bed')
    os.system('rm -f __tri.tsv')
    os.system('rm -f __tmpbed.sorted.uniquified.bed')

    return(tridict)

def PrepareClusteredEventCheck(rawmuts):
    chroms = list(set([m.split('\t')[0] for m in rawmuts]))
    chromPos = dict.fromkeys(chroms)
    for item in chromPos:
        chromPos[item] = []
    for m in rawmuts:
        chr = m.split('\t')[0]
        chromPos[chr].append(int(m.split('\t')[1]))
    return(chromPos)

def ClusterMutCount(num, arr, Options):
    '''
    Helper function for Checking clustered mutational events. It returns the count of clustered events around a single mutation.

    :param num: A loci position from within a chrom
    :param arr: List of positions on that chromosome
    :return: A count of mutations +/- the num
    '''
    count = 0
    for index in range (len (arr)):
        if abs(num - arr[index]) < Options.clusterFlank and abs(num - arr[index])!=0:
            count+=1
    return(count)

def GatherMutectSamples(Options):
    '''
    Obtains the samples within the specified directory and maps to a sample map if one is provided.

    :param Options: Parser arguments
    :return: A list of files (length = 1 if no sample mapping, list = 2 if sample mapping [0] is a list [1] of the samples
    '''
    infiles = glob.glob(Options.inDir+'*.passed.vcf')
    if Options.verbose:
        logging.info("A total of %s vcf files have been found..."%(len(infiles)))
    logging.info("...")
    logging.info("...")
    if Options.usesamplemap==False:
        return([os.path.abspath(thefile) for thefile in infiles])
    elif Options.usesamplemap==True:
        with open(Options.samplemap, 'r') as samMap:
            samMap = [sam.replace('\n','').split(',') for sam in samMap.readlines()]
        return([samMap, [os.path.abspath(thefile) for thefile in infiles]])

def GatherDragenSamples(Options):
    '''
        Obtains the samples within the specified directory and maps to a sample map if one is provided.

        :param Options: Parser arguments
        :return: A list of files (length = 1 if no sample mapping, list = 2 if sample mapping [0] is a list [1] of the samples
    '''
    infiles = glob.glob(Options.inDir + '*.vcf.gz')
    if Options.verbose:
        logging.info("A total of %s vcf files have been found..." % (len(infiles)))
    logging.info("...")
    logging.info("...")
    if Options.usesamplemap == False:
        return ([os.path.abspath(thefile) for thefile in infiles])
    elif Options.usesamplemap == True:
        with open(Options.samplemap, 'r') as samMap:
            samMap = [sam.replace('\n', '').split(',') for sam in samMap.readlines()]
        return ([samMap, [os.path.abspath(thefile) for thefile in infiles]])

def WriteFilteredVCFFiles(outpath, samname, header, rawmuts, flags): # Currently only for mutect2
    '''
    Final filtered (but not annotated vcf) file writer.
    :param outpath: Path to write the samples to
    :param samname: Sample identifier name
    :param header: Header of the vcf file as a string
    :param rawmuts: All of the mutations in a list unfiltered
    :param flags: All of the flags for each mutation
    :return: None
    '''
    fullFail = ['@','*','&','!','^','`',';'] # All except '?'
    outputName = outpath + samname.replace("passed.vcf","passed.variantfilter.vcf")
    with open(outputName, 'w') as outFile:
        outFile.write(header)
        for i, flagIDs in enumerate(flags):
            if any(f in flagIDs for f in fullFail) == False:
                outFile.write(rawmuts[i]+'\n')

def AnnotateVCFFiles(inputpath, samname, Options):
    inputName = inputpath + samname.replace("passed.vcf","passed.variantfilter.vcf")
    outputName = inputpath + samname.replace("passed.vcf","passed.variantfilter.ann.vcf")
    cmd = 'java -Xmx8G -jar %s ann -noStats -canon %s %s > %s' % (Options.jarfile, Options.ref, inputName, outputName)
    os.system(cmd)
    if os.path.exists(outputName) and os.path.exists(inputName):
        os.remove(inputName)

def main():
    SetupEnvironment()
    localpath = os.path.abspath(__file__).replace('BuildVariantScript.py', '')  # path to scripts working directory
    Options = Parser()

    theFlags = {'@': 'Failing variant reads in normal check', '?': 'VAF > 0.05 in normal', '*': 'Failing normal reads check',
     '&': 'Failing variant supporting reads check',
     '!': 'Failing VAF < cutoff', '^': "Failing events check", '`': 'Failing variant total depth check',
     ';': 'Failing clustered event check', 'G':"Variant called in germline"}

    if Options.caller==1:
        inFiles = GatherMutectSamples(Options)
        MuTect2Filter(Options, inFiles, theFlags)
    elif Options.caller==2:
        Strelka2Main(Options, theFlags)
    elif Options.caller==3:
        pass
    elif Options.caller==4:
        pass
    elif Options.caller==5:
        DragenMain(Options, theFlags)
    else:
        logging.error("Unspecified caller.")
        sys.exit()

    logging.info("Complete.")

if __name__=="__main__":
    main()
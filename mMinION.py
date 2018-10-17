from sqlalchemy.ext.declarative import declarative_base
import datetime
import subprocess
import vcf
from types import FunctionType
import logging
import json
import argparse
import numpy
import pandas
from scipy.stats import binned_statistic

with open('/Users/mattparker/PycharmProjects/mMinION/conf/config.json') as json_data_file:
    config = json.load(json_data_file)

def decorator(f,):
    """
    decorates functions with a logger so that this process is automated
    :param f: function
    :return: wrapper
    """
    def wrapper(self, *args, **kwargs):
        self.logger.info("Running '{0}' with args: '{1}'".format(f.__name__.upper(),args))
        start = datetime.datetime.now()
        try:
            result = f(self, *args, **kwargs)
            end =datetime.datetime.now()
            duration = end - start
            self.logger.info("Finnished '{0}' in {1}".format(f.__name__.upper(),str(duration)))
            return result
        except Exception as e:
            self.logger.warning("ERROR in {0}: {1}".format(f.__name__,e))
            exit()


    return wrapper


class MetaClass(type):
    def __new__(meta, classname, bases, classDict):
        """
        decorates all functions in the class with the decorator defined above

        :param classname:
        :param bases:
        :param classDict:
        :return:
        """
        newClassDict = {}
        for attributeName, attribute in classDict.items():
                if isinstance(attribute, FunctionType):
                    if attributeName != "__init__" and attributeName != "run_command":
                        # replace it with a wrapped version
                        print attributeName
                        attribute = decorator(attribute)
                newClassDict[attributeName] = attribute
        return type.__new__(meta, classname, bases, newClassDict)


class mMinION():
    __metaclass__ = MetaClass
    base = declarative_base()

    def __init__(self,samplename,fastq,output,sangerresults):
        self.output = output
        self.samplename = samplename
        self.fastq = fastq
        self.sangerresults = sangerresults
        # set up logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # complete setup of logger
        handler = logging.FileHandler(self.output + "/pipeline.log", 'w')
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)

        self.logger.addHandler(handler)
        self.logger.addHandler(console)

    def run_command(self, cmd):
        self.logger.info("Command: " + cmd)
        ##todo chnage to popen
        try:
            result = subprocess.call(cmd, shell=True)
            if result != 0:
                exit(result)
        except subprocess.CalledProcessError as e:
            print('Error executing command: ' + str(e.returncode))
            exit(1)

    def correction_trimming(self):

        #correct reads
        command = [config["software"]["canu"],"-correct","-p",self.samplename,"-d",self.output,"genomeSize=16.5k overlapper=mhap utgReAlign=true","-nanopore-raw",self.fastq]
        self.run_command(" ".join(command))

        #trim reads
        corrected = self.output + "/" + self.samplename + ".correctedReads.fasta.gz"
        command = [config["software"]["canu"], "-trim","-p",self.samplename,"-d",self.output,"genomeSize=16.5k overlapper=mhap utgReAlign=true","-nanopore-corrected",corrected]
        self.run_command(" ".join(command))

        trimmed = self.output + "/" + self.samplename + ".trimmedReads.fasta.gz"
        return trimmed

    def mapping(self,fastq):

        rg = r'"@RG\tID:'+self.samplename+r'\tSM:mMinION"'

        outfile = self.output + "/" + self.samplename + ".sorted.bam"
        print config["software"]
        command = [config["software"]["minimap2"], "-a", '-R',rg, config["reference"]["genome"], fastq,"|",config["software"]["samtools"],"sort","-","|",config["software"]["samtools"],"view","-bh","-",">",outfile]
        self.run_command(" ".join(command))

        command = [config["software"]["samtools"],"index",outfile]
        self.run_command(" ".join(command))

        return outfile

    def coverage(self,bam):
        """
        this uses samtools, sambamba fails ti deal with gapped reads very well
        """
        outfile = self.output + "/" + self.samplename + ".sorted.bam.depth"
        command = [config["software"]["samtools"],"depth","-a", "-Q 20",bam,">",outfile]
        self.run_command(" ".join(command))



        return outfile

    def transform_coverage(self,depth):
        out = []
        position = []
        coverage = []
        with open(depth,"rb") as f:
            for line in f:

                chrom,pos,cov = line.rstrip().split("\t")
                position.append(int(pos))
                coverage.append(int(cov))
                new_line = [chrom,str(pos),str(pos),str(cov)]
                out.append("\t".join(new_line))


        df = pandas.DataFrame({'position':position,'coverage':coverage})
        print coverage

        self.max_coverage = max(coverage)+100
        print binned_statistic(coverage, coverage, statistic='mean', bins=1650, range=None)



        circos_file = self.output + "/" + self.samplename + ".sorted.bam.depth.circos"
        circos_out = open(circos_file,"w")

        circos_out.write("\n".join(out))
        circos_out.close()

        return circos_file

    def sv(self,bam):

        #run nanosv

        outfile = self.output + "/" + self.samplename + ".sorted.bam.vcf"
        command = [config["software"]["nanosv"],"-o",outfile,"-s",config["software"]["samtools"],"-b",config["reference"]["bed"],"-t","6",bam]
        self.run_command(" ".join(command))

        #filter calls which overlap the uncovered region around 0kb

        return outfile

    def transform_sv(self,sv):
        vcf_reader = vcf.Reader(open(sv, 'r'))
        out = []
        for record in vcf_reader:
            print record
            print record.ALT[0]
            print type(record.ALT[0])

            if str(record.ALT[0]) == "<INS>":
                print record.INFO["END"]
                alt=str(int(record.POS)+int(record.INFO["END"]))
            else:
                alt = str(''.join(i for i in str(record.ALT[0]) if i.isdigit()))
            if record.QUAL < 200:
                color="lgrey"
            else:
                color="red"
            if "LowQual" in record.FILTER:
                color="vvlgrey"
            line = [record.CHROM,str(record.POS),str(record.POS),record.CHROM,alt,alt,"color="+color+",thickness=5p"]
            out.append("\t".join(line))

        links_file = self.output + "/" + self.samplename + ".sorted.bam.sv.circos"
        circos_out = open(links_file, "w")

        circos_out.write("\n".join(out))
        circos_out.close()

        return links_file


    def parse_sanger_svs(self):
        sanger_file = self.output + "/" + self.samplename + ".labels.circos"
        sanger_out = open(sanger_file, "w")

        print self.sangerresults

        with open ("/Users/mattparker/PycharmProjects/mMinION/conf/labels.txt", "rb") as l:
            for line in l:
                sanger_out.write(line)

        with open(self.sangerresults, "rb") as f:
            for line in f:
                sample, miseqid, barcode, deletion_status, deletion_start, deletion_end, mutation, mutation_frequency = line.rstrip().split(
                    "\t")
                if self.samplename in line:
                    if deletion_start != "NA":
                        sanger_out.write(str("MT ")+str(deletion_start)+str(" ")+str(deletion_start)+str(" SANGER_1 link_color=blue\n"))
                        sanger_out.write(str("MT ") + str(deletion_end) + str(" ") + str(deletion_end) + str(" SANGER_2 link_color=blue\n"))

        sanger_out.close
        return sanger_file

    def mpileup(self,bam):
        outfile = self.output + "/" + self.samplename + ".mpileup"
        command = [config["software"]["samtools"], "mpileup", bam, ">", outfile]
        self.run_command(" ".join(command))

        return outfile


    def call_variants(self,mpileup):
        snp_outfile = self.output + "/" + self.samplename + ".snp.vcf"
        command = ["java","-jar",config["software"]["varscan"],"mpileup2snp","mpileup --output-vcf 1 --strand-filter 0",">",snp_outfile]
        self.run_command(" ".join(command))

        #do indels not follow VCF spec? process them so they do
        indel_outfile = self.output + "/" + self.samplename + ".indel.vcf"
        command = ["java", "-jar", config["software"]["varscan"], "mpileup2indel",
                   "mpileup --output-vcf 1 --strand-filter 0", ">", indel_outfile]
        self.run_command(" ".join(command))

    def annotate_variants(self,vcf):
        outfile = self.output + "/" + self.samplename + ".annotated.vcf"
        command = [config["software"]["vcfanno"],config["reference"]["vcfanno"],vcf,">",outfile]
        self.run_command(" ".join(command))

        return outfile

    def prepare_ciros(self,links,coverage,sanger):
        #replace lines in conf with actual locations
        print self.max_coverage
        conf_file = self.output + "/" + self.samplename + ".circos.conf"
        circos_out = open(conf_file, "w")

        with open("/Users/mattparker/PycharmProjects/mMinION/conf/circos.conf", "rb") as f:
            for line in f:
                line = line.replace("{{LINKS}}",links)
                line = line.replace("{{COVERAGE}}",coverage)
                line = line.replace("{{MAX}}", str(self.max_coverage))
                line = line.replace("{{LABELS}}",sanger)
                circos_out.write(line)

        circos_out.close()
        command = ["cp","/Users/mattparker/PycharmProjects/mMinION/conf/ticks.conf",self.output+"/ticks.conf"]
        self.run_command(" ".join(command))
        command = ["cp","/Users/mattparker/PycharmProjects/mMinION/conf/ideogram.conf",self.output+"/ideogram.conf"]
        self.run_command(" ".join(command))
        command = ["cp", "/Users/mattparker/PycharmProjects/mMinION/conf/karyotype.txt", self.output + "/karyotype.txt"]
        self.run_command(" ".join(command))
        return conf_file


    def circos(self,conf):

        circos_file = self.samplename+".circos"
        command = [config["software"]["circos"], "-conf", conf,"-outputdir",self.output,"-outputfile",circos_file,"-debug_group text"]
        self.run_command(" ".join(command))




def main():

    parser = argparse.ArgumentParser(description='mMinION - Run fastq file for a sample through mitochondraial analysis pipeline')
    parser.add_argument('--samplename',help='The name of the sample being processed - files will take this name')
    parser.add_argument('--fastq',help='Full path to the Fastq file')
    parser.add_argument('--outdir',help='The output directory for all files generated')
    parser.add_argument('--bam',help='Provide bam file to skip mapping')
    parser.add_argument('--coverage',help='Provide samtools depth output to skip this step')
    parser.add_argument('--sv',help='Provide SV calls to skip SV calling')
    parser.add_argument('--sangerresults')

    parser.add_argument('--correction', action='store_true', help='Turn on read correction & trimming with canu. WARNING! This is computationally intensive')

    args = parser.parse_args()



    m = mMinION(args.samplename,args.fastq,args.outdir,args.sangerresults)
    if args.bam:
        bam = args.bam
    else:
        if args.correction:
            trimmed_fastq = m.correction_trimming()
            bam = m.mapping(trimmed_fastq)
        else:
            bam = m.mapping(args.fastq)

    if args.coverage:
        coverage = args.coverage
    else:
        coverage = m.coverage(bam)

    #coverage = "/Users/mattparker/Documents/Projects/Wood/DATA/MINION/PORETOOLS_FASTQ/test/test.sorted.bam.depth"
    circos_coverage = m.transform_coverage(coverage)
    # sv = m.sv(bam)

    if args.sv:
        sv = args.sv
    else:
        sv = sv = m.sv(bam)


    #sv = "/Users/mattparker/Documents/Projects/Wood/DATA/MINION/PORETOOLS_FASTQ/test/test.sorted.bam.vcf"


    links = m.transform_sv(sv)
    sanger = m.parse_sanger_svs()

    circos_conf = m.prepare_ciros(links, circos_coverage,sanger)
    circos = m.circos(circos_conf)
    print circos_conf

    #preprocessing steps
    #convert to fastq
    #poretools fastq /Volumes/Untitled/20180206_1708_Multiplex_Exp/uploaded/fast5/pass/ > ~/Documents/Projects/Wood/DATA/MINION/PORETOOLS_FASTQ/20180206_1708_Multiplex_Exp.fastq
    #porechop -i PORETOOLS_FASTQ/20180206_1708_Multiplex_Exp.fastq -b ~/Documents/Projects/Wood/DATA/MINION/PORETOOLS_FASTQ/ --barcode_threshold 85 -t 16

if __name__ == '__main__':
    main()
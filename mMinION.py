from sqlalchemy.ext.declarative import declarative_base
import vcf
import logging
import json
import argparse
from lib.common import Common
import os
from lib.common import MetaClass

class mMinION(Common):
    __metaclass__ = MetaClass
    base = declarative_base()

    def __init__(self,samplename,fastq,output,sangerresults):
        self.output = output
        self.samplename = samplename
        self.fastq = fastq
        self.sangerresults = sangerresults
        with open('/Users/mattparker/PycharmProjects/mMinION/conf/config.json') as json_data_file:
            self.config = json.load(json_data_file)

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


    def correction_trimming(self,corrected_error_rate,min_read_length):

        #correct reads
        command = [self.config["software"]["canu"],"-correct","-p",self.samplename,"-d",self.output,"genomeSize=16.5k overlapper=mhap utgReAlign=true stopOnReadQuality=false -nanopore-raw",self.fastq]
        self.run_command(" ".join(command))

        #trim reads
        corrected = self.output + "/" + self.samplename + ".correctedReads.fasta.gz"
        command = [self.config["software"]["canu"], "-trim","-p",self.samplename,"-d",self.output,"genomeSize=16.5k overlapper=mhap utgReAlign=true stopOnReadQuality=false -nanopore-corrected",corrected]
        self.run_command(" ".join(command))

        trimmed = self.output + "/" + self.samplename + ".trimmedReads.fasta.gz"
        return trimmed

    def trim_adapters(self,fastqs):
        out_files = []
        fastq = fastqs.split(",")
        for i in fastq:
            output = self.output + "/" + os.path.basename(i).replace("fastq","trimmed.fastq.gz")
            out_files.append(output)
            command = [self.config["software"]["porechop"], "-i", i, "-o",output]
            self.run_command(" ".join(command))
        return out_files


    def mapping(self,fastq):

        rg = r'"@RG\tID:'+self.samplename+r'\tSM:mMinION"'

        outfile = self.output + "/" + self.samplename + ".sorted.bam"
        print self.config["software"]
        command = [self.config["software"]["minimap2"], "--MD -ax map-ont", '-R',rg, self.config["reference"]["genome"], fastq,"|",self.config["software"]["samtools"],"sort","-","|",self.config["software"]["samtools"],"view","-bh","-",">",outfile]
        self.run_command(" ".join(command))

        command = [self.config["software"]["samtools"],"index",outfile]
        self.run_command(" ".join(command))

        return outfile

    def filter(self,bam):
        #samtools view -h N4.trimmed.sorted.bam | awk 'length($10) > 1000 || $1 ~ /^@/' | samtools view -bS > N4.trimmed.sorted.filtered.bam

        outfile = self.output + "/" + self.samplename + ".trimed.sorted.filtered.bam"
        command = [self.config["software"]["samtools"],"view","-h",bam,"| awk 'length($10) > 200 || $1 ~ /^@/' |",self.config["software"]["samtools"],"view -bS >",outfile]
        self.run_command(" ".join(command))

        command = [self.config["software"]["samtools"], "index", outfile]
        self.run_command(" ".join(command))

        return outfile

    def sv(self,bam):

        #run nanosv

        outfile = self.output + "/" + self.samplename + ".sorted.bam.vcf"
        command = [self.config["software"]["nanosv"],"-o",outfile,"-s",self.config["software"]["samtools"],"-b",self.config["reference"]["bed"],"-t","6",bam]
        self.run_command(" ".join(command))

        #filter calls which overlap the uncovered region around 0kb

        return outfile

    def transform_sv(self, sv):
        vcf_reader = vcf.Reader(open(sv, 'r'))
        out = []
        for record in vcf_reader:
            print record
            print record.ALT[0]
            print type(record.ALT[0])

            if str(record.ALT[0]) == "<INS>":
                print record.INFO["END"]
                alt = str(int(record.POS) + int(record.INFO["END"]))
            else:
                alt = str(''.join(i for i in str(record.ALT[0]) if i.isdigit()))
            if record.QUAL < 200:
                color = "lgrey"
            else:
                color = "red"
            if "LowQual" in record.FILTER:
                color = "vvlgrey"
            line = [record.CHROM, str(record.POS), str(record.POS), record.CHROM, alt, alt,
                    "color=" + color + ",thickness=5p"]
            out.append("\t".join(line))

        links_file = self.output + "/" + self.samplename + ".sorted.bam.sv.circos"
        circos_out = open(links_file, "w")

        circos_out.write("\n".join(out))
        circos_out.close()

        return links_file

    def sniffles_sv(self,bam):

        outfile = self.samplename + ".sorted.bam.sniffles.vcf"
        bam_path = os.path.dirname(bam)
        bam = os.path.basename(bam)

        self.docker_command(image=self.config["software"]["sniffles_image"],executable=self.config["software"]["sniffles_path"],data_dir=bam_path,output_dir=self.output,arguments="-q 30 --allelefreq 0.01 -v /output/"+outfile+" -m /data/"+bam)

        return self.output + "/" + self.samplename + ".sorted.bam.sniffles.vcf"

    def transform_sniffles_sv(self, sv):
        vcf_reader = vcf.Reader(open(sv, 'r'))
        red_out = []
        grey_out = []
        for record in vcf_reader:

            if str(record.ALT[0]) == '<DEL>':

                if "IMPRECISE" in record.INFO:
                    if record.INFO["IMPRECISE"] == True:
                        color = "vvlgrey"

                if "PRECISE" in record.INFO:
                    if record.INFO["PRECISE"] == True:
                        color = "red"

                alt = str(record.INFO["END"])
                # print record.FORMAT
                # #reference reads
                # reference_reads = record.samples[0]["DR"]
                # #variant reads
                # variant_reads = record.samples[0]["DV"]
                # total_reads = reference_reads+variant_reads
                #
                #
                # print record.samples[0]["GT"]
                #
                # if reference_reads == 0:
                #     color="vvlgrey"
                # else:
                #     color = "red"
                # print total_reads
                # print self.median_coverage
                # print (total_reads / float(self.median_coverage))
                # if (total_reads / float(self.median_coverage)) < 0.1:
                #     print "low_cov"
                #     color="vvlgrey"
                #
                line = [record.CHROM, str(record.POS), str(record.POS), record.CHROM, alt, alt,
                        "color=" + color + ",thickness=5p"]
                if color == "vvlgrey":
                    grey_out.append("\t".join(line))
                elif color == "red":
                    red_out.append("\t".join(line))

        links_file = self.output + "/" + self.samplename + ".sorted.bam.sv.circos"
        circos_out = open(links_file, "w")

        circos_out.write("\n".join(grey_out))
        circos_out.write("\n")
        circos_out.write("\n".join(red_out))
        circos_out.write("\n")
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


    def call_variants(self,mpileup):
        snp_outfile = mpileup + ".snp.vcf"
        command = ["java","-jar",self.config["software"]["varscan"],"mpileup2snp",mpileup," --min-var-freq 0.01 --output-vcf 1 --strand-filter 0",">",snp_outfile]
        self.run_command(" ".join(command))

        #do indels not follow VCF spec? process them so they do
        indel_outfile = mpileup + ".indel.vcf"
        command = ["java", "-jar", self.config["software"]["varscan"], "mpileup2indel",
                   mpileup," --output-vcf 1 --strand-filter 0", ">", indel_outfile]
        self.run_command(" ".join(command))

        bgzip_snps = self.bgzip_tabix(snp_outfile)
        bgzip_indels = self.bgzip_tabix(indel_outfile)

        combined = mpileup + ".snp.indel.vcf"
        command = [self.config["software"]["bcftools"],"merge","--force-samples",bgzip_snps,bgzip_indels,">",combined]
        self.run_command(" ".join(command))

        return combined

    def annotate_variants(self,vcf):
        outfile = self.output + "/" + self.samplename + ".sorted.bam.mpileup.snp.indel.annotated.vcf"
        command = [self.config["software"]["vcfanno"],self.config["reference"]["vcfanno"],vcf,">",outfile]
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




def main():

    parser = argparse.ArgumentParser(description='mMinION - Run fastq file for a sample through mitochondraial analysis pipeline')
    parser.add_argument('--samplename',help='The name of the sample being processed - files will take this name')
    parser.add_argument('--fastq',help='Full path to the Fastq file')
    parser.add_argument('--trimmed_fastq', help='Full path to the Fastq file')
    parser.add_argument('--outdir',help='The output directory for all files generated')
    parser.add_argument('--bam',help='Provide bam file to skip mapping')
    parser.add_argument('--coverage',help='Provide samtools depth output to skip this step')
    parser.add_argument('--sv',help='Provide SV calls to skip SV calling')
    parser.add_argument('--variants', help='Provide SNV/INDEL calls to skip calling')
    parser.add_argument('--correctedErrorRate')
    parser.add_argument('--minReadLength')
    parser.add_argument('--sangerresults')

    parser.add_argument('--correction', action='store_true', help='Turn on read correction & trimming with canu. WARNING! This is computationally intensive')

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)


    m = mMinION(args.samplename,args.fastq,args.outdir,args.sangerresults)
    if args.bam:
        bam = m.filter(args.bam)
    else:
        if args.correction:
            trimmed_fastq = m.correction_trimming(corrected_error_rate=args.correctedErrorRate,min_read_length=args.minReadLength)
            bam = m.mapping(trimmed_fastq)
        else:
            if args.trimmed_fastq:
                trimmed_fastq = args.trimmed_fastq.split(",")
            else:
                trimmed_fastq = m.trim_adapters(args.fastq)
            fastq = " ".join(trimmed_fastq)
            temp_bam = m.mapping(fastq)
            bam = m.filter(temp_bam)


    if args.coverage:
        coverage = args.coverage
    else:
        coverage = m.coverage(bam)


    circos_coverage = m.transform_coverage(coverage)

    if args.sv:
        sv = args.sv
    else:
        sv = m.sv(bam)

    if args.variants:
        variants = args.variants
    else:
        pileup = m.mpileup(bam)
        raw_variants = m.call_variants(pileup)
        variants = m.annotate_variants(raw_variants)



    links = m.transform_sniffles_sv(sv)
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
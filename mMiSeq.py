from sqlalchemy.ext.declarative import declarative_base
import argparse
import os
from lib.common import MetaClass,Common
import logging
import json


class mMiSeq(Common):
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

def main():

    parser = argparse.ArgumentParser(description='mMinION - Run fastq file for a sample through mitochondraial analysis pipeline')
    parser.add_argument('--samplename',help='The name of the sample being processed - files will take this name')
    parser.add_argument('--fastq',help='Full path to the Fastq file')
    parser.add_argument('--outdir',help='The output directory for all files generated')
    parser.add_argument('--bam',help='Provide bam file to skip mapping')
    parser.add_argument('--coverage',help='Provide samtools depth output to skip this step')
    parser.add_argument('--sv',help='Provide SV calls to skip SV calling')
    parser.add_argument('--variants', help='Provide SNV/INDEL calls to skip calling')
    parser.add_argument('--correctedErrorRate')
    parser.add_argument('--minReadLength')
    parser.add_argument('--sangerresults')



    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)


    m = mMiSeq(args.samplename,args.fastq,args.outdir,args.sangerresults)
    if args.bam:
        bam = args.bam
    else:
        bam = m.mapping(args.fastq)

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



    links = m.transform_sv(sv)
    sanger = m.parse_sanger_svs()

    circos_conf = m.prepare_ciros(links, circos_coverage,sanger)
    circos = m.circos(circos_conf)
    print circos_conf

    # preprocessing steps
    # convert to fastq
    # poretools fastq /Volumes/Untitled/20180206_1708_Multiplex_Exp/uploaded/fast5/pass/ > ~/Documents/Projects/Wood/DATA/MINION/PORETOOLS_FASTQ/20180206_1708_Multiplex_Exp.fastq
    # porechop -i PORETOOLS_FASTQ/20180206_1708_Multiplex_Exp.fastq -b ~/Documents/Projects/Wood/DATA/MINION/PORETOOLS_FASTQ/ --barcode_threshold 85 -t 16

if __name__ == '__main__':
    main()
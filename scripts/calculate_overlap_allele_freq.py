import vcf
import argparse
from collections import defaultdict
import json
import seaborn as sns
import pandas
def main():

    parser = argparse.ArgumentParser(description='produce overlap table for use in plotting')
    parser.add_argument('--miseq',help='miseq vcf file')
    parser.add_argument('--minion',help='minion vcf file')
    parser.add_argument('--sample', help='samplename')
    parser.add_argument('--output')

    args = parser.parse_args()

    miseq_vcf = vcf.Reader(open(args.miseq, 'r'))

    #result is structured like so]
    # pos>ref>alt>platform>{genotype:,frequency}
    result = {}

    for record in miseq_vcf:

        GT = record.samples[0]['GT']
        GQ = record.samples[0]['GQ']
        ABQ = record.samples[0]['ABQ']
        FREQ = float(record.samples[0]['FREQ'].replace("%",""))
        if record.POS not in result:
            result[record.POS] = {}
            result[record.POS][record.REF] = {}
            result[record.POS][record.REF][str(record.ALT[0])] = {}
            result[record.POS][record.REF][str(record.ALT[0])]["miseq"] = { "genotype":GT,"frequency":FREQ, "genotype_quality": GQ, "average_base_quality": ABQ }


    minion_vcf = vcf.Reader(open(args.minion, 'r'))

    for record in minion_vcf:

        GT = record.samples[0]['GT']
        GQ = record.samples[0]['GQ']
        ABQ = record.samples[0]['ABQ']
        FREQ = float(record.samples[0]['FREQ'].replace("%", ""))
        if record.POS not in result:
            result[record.POS] = {}
            result[record.POS][record.REF] = {}
            result[record.POS][record.REF][str(record.ALT[0])] = {}
            result[record.POS][record.REF][str(record.ALT[0])]["minion"] = {"genotype": GT, "frequency": FREQ, "genotype_quality": GQ, "average_base_quality": ABQ}
        else:
            if record.REF not in result[record.POS]:
                result[record.POS][record.REF] = {}
                result[record.POS][record.REF][str(record.ALT[0])] = {}
                result[record.POS][record.REF][str(record.ALT[0])]["minion"] = {"genotype": GT, "frequency": FREQ, "genotype_quality": GQ, "average_base_quality": ABQ}
            else:
                if str(record.ALT[0]) not in result[record.POS][record.REF]:
                    result[record.POS][record.REF][str(record.ALT[0])] = {}
                    result[record.POS][record.REF][str(record.ALT[0])]["minion"] = {"genotype": GT, "frequency": FREQ, "genotype_quality": GQ, "average_base_quality": ABQ}
                else:
                    result[record.POS][record.REF][str(record.ALT[0])]["minion"] = {"genotype": GT, "frequency": FREQ, "genotype_quality": GQ, "average_base_quality": ABQ}


    output = ["\t".join(["Sample","POS","REF","ALT","miseq_GT","miseq_GQ","miseq_ABQ","miseq_FREQ","minion_GT","minion_GQ","minion_ABQ","minion_FREQ"])]
    for pos in result:
        line = []
        for ref in result[pos]:
            for alt in result[pos][ref]:
                line.append(args.sample)
                line.append(str(pos))
                line.append(ref)
                line.append(alt)

                if "miseq" in result[pos][ref][alt]:
                    line.append(result[pos][ref][alt]["miseq"]["genotype"])
                    line.append(str(result[pos][ref][alt]["miseq"]["genotype_quality"]))
                    line.append(str(result[pos][ref][alt]["miseq"]["average_base_quality"]))
                    line.append(str(result[pos][ref][alt]["miseq"]["frequency"]))
                else:
                    line.append(".")
                    line.append("0")
                    line.append("0")
                    line.append("0")
                if "minion" in result[pos][ref][alt]:
                    line.append(result[pos][ref][alt]["minion"]["genotype"])
                    line.append(str(result[pos][ref][alt]["minion"]["genotype_quality"]))
                    line.append(str(result[pos][ref][alt]["minion"]["average_base_quality"]))
                    line.append(str(result[pos][ref][alt]["minion"]["frequency"]))
                else:
                    line.append(".")
                    line.append("0")
                    line.append("0")
                    line.append("0")


        output.append("\t".join(line))

    f = open(args.output,"w")
    f.write("\n".join(output))
    f.close()


if __name__ == '__main__':
    main()
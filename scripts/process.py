import sys
import argparse

parser = argparse.ArgumentParser(description="Processing variant file to identifie the passed variant", usage="%(prog)s [options]",
formatter_class=argparse.ArgumentDefaultsHelpFormatter, argument_default=argparse.SUPPRESS)

parser.add_argument("input", help="Input file from clariovant", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument("output", help="The output file from clariovant", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
parser.add_argument("-f", "--filter", help="Minimum threshold for variant to be passed (default: %(default)s)", type=int, default=200 )

args = parser.parse_args()

myFile = args.input
dataOut = args.output
threshold = args.filter


# myFile = sys.argv[1]
# with open(myFile, "r") as dataIn, open(myFile+"_filter.vcf", 'w') as dataOut:

for line in myFile:
        lineSplit = line.split()
        if line.startswith("#"):
            dataOut.write(line)
        elif lineSplit[4].startswith("<"):
            pass
        else:
            if int(float(lineSplit[5])) >= threshold:
                lineSplit[6] = 'PASS'
                lineSplit[5] = str(int(float(lineSplit[5])))
                dataOut.write("{}\n".format("\t".join(lineSplit)))

myFile.close()
dataOut.close()

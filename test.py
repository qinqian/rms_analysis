import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--foo', help='foo of the %(prog)s program')
args = parser.parse_args()

print(args)




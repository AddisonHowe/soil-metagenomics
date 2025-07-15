"""Main script.

"""

import sys
import argparse

def parse_args(args):
    parser = argparse.ArgumentParser()
    return parser.parse_args(args)

def main(args):
    print("Hello world")
    return


#######################
##  Main Entrypoint  ##
#######################

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)

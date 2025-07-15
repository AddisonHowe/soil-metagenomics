import sys

def main():
    from mgsa.main import parse_args, main
    args = parse_args(sys.argv[1:])
    main(args)

if __name__ == "__main__":
    main()

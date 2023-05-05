from . import __version__
import argparse
import sys
from .cellatlas_build import setup_build_args, validate_build_args

# Steps to add new subcommands
# Create cellatlas_subcommand.py (create setup_subcmd_args, validate_subcmd_args, run_subcmd in that file)
# (in this file) from cellatlas_subcmd import setup_subcmd_args, validate_subcmd_args
# Add setup_subcmd_args to command_to_parser along with its key==str(subcmd)
# Add validate_subcmd_args to COMMAND_TO_FUNCTION along with its key==str(subcmd)


def main():
    # setup parsers
    parser = argparse.ArgumentParser(
        description=f"cellatlas {__version__}: Format sequence specification files"
    )
    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<CMD>",
    )

    # Setup the arguments for all subcommands
    command_to_parser = {
        "build": setup_build_args(subparsers),
    }

    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Setup validator and runner for all subcommands (validate and run if valid)
    COMMAND_TO_FUNCTION = {
        "build": validate_build_args,
    }
    COMMAND_TO_FUNCTION[sys.argv[1]](parser, args)


if __name__ == "__main__":
    main()

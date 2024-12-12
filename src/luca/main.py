import click

from . import __version__
from .count import count
from .merge import merge


@click.group(invoke_without_command=False)
@click.version_option(__version__)
def main():
    pass


main.add_command(count)
main.add_command(merge)

import click
from mutation_filtering.analysis import run_analysis

@click.command()
@click.argument("ivar_dir", type=click.Path(exists=True, file_okay=False))
@click.argument("output_dir", type=click.Path())
@click.option("--workers", default=4, help="Number of parallel workers")
@click.option("--log-level", default="INFO", help="Logging level")
@click.option("--strict", is_flag=True, default=False, help="Strict error handling for missing data")
@click.option("--suspicious-file", type=click.Path(exists=True), help="CSV of suspicious mutations")
def cli(ivar_dir, output_dir, workers, log_level, strict, suspicious_file):
    """Main entry point to run mutation filtering pipeline."""
    run_analysis(ivar_dir, output_dir, workers, log_level, strict, suspicious_file)

if __name__ == "__main__":
    cli()

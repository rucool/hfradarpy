#!/usr/bin/env python

"""Tests for `hfradarpy` package."""

import pytest

from click.testing import CliRunner

import hfradarpy
from hfradarpy import cli
from pathlib import Path

data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
output_path = (Path(__file__).parent.with_name("examples") / "output").resolve()


def test_command_line_interface():
    """Test the CLI."""
    radial_dir = data_path / "radials" / "ruv" / "SEAB" 
    path_save = output_path / "cli_output" 
    point1 = "-73.63, 40.29"
    point2 = "-73.63, 40.31"
    runner = CliRunner()
    result = runner.invoke(
        cli.extract_timeseries, [
            "--point", point1, 
            "--point", point2,
            "--path_data", str(radial_dir),
            "--path_save", str(path_save),
            "--type", 'csv']
    )
    assert result.exit_code == 0
    assert "hfradarpy.cli.extract_timeseries" in result.output
    help_result = runner.invoke(cli.main, ["--help"])
    assert help_result.exit_code == 0
    assert "--help  Show this message and exit." in help_result.output

"""L0: the command line lists every command and loads the example configuration."""
import subprocess
import sys
from pathlib import Path

import pytest

from plbuild.cli import COMMANDS, main

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "bianconi2012" / "plbuild.yaml"


def test_help_lists_all_commands():
    out = subprocess.run([sys.executable, "-m", "plbuild", "--help"], capture_output=True, text=True, check=True).stdout
    for name in ("search", "fetch", "build", "validate", "explain"):
        assert name in out
    assert set(COMMANDS) == {"search", "fetch", "build", "validate", "explain"}


@pytest.mark.parametrize("cmd", ["fetch", "build", "validate"])
def test_stub_commands_load_config_then_report_not_implemented(cmd, capsys):
    assert main([cmd, "--config", str(EXAMPLE)]) == 2
    assert "not implemented" in capsys.readouterr().err


def test_bad_config_exits_1_with_file_and_line(tmp_path, capsys):
    cfg = tmp_path / "plbuild.yaml"
    cfg.write_text("vocabulary: v.yaml\ncompartments: {cell: {}}\nmodels: [BIOMD0000000427]\nbogus: 1\n")
    (tmp_path / "v.yaml").write_text("X:\n  class: molecule\n  unit: item\n")
    assert main(["build", "--config", str(cfg)]) == 1
    assert f"{cfg}:4:" in capsys.readouterr().err

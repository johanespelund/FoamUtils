"""Tests for FoamUtils.file_utils helper functions."""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from FoamUtils.file_utils import get_sorted_times, load


# ── get_sorted_times ──────────────────────────────────────────────────────────

class TestGetSortedTimes:
    """Tests for get_sorted_times – returns time-directory names sorted numerically."""

    def test_returns_sorted_list(self, tmp_path):
        for d in ("0.5", "0", "1", "0.1", "10"):
            (tmp_path / d).mkdir()
        result = get_sorted_times(str(tmp_path))
        assert result == ["0", "0.1", "0.5", "1", "10"]

    def test_single_entry(self, tmp_path):
        (tmp_path / "0").mkdir()
        result = get_sorted_times(str(tmp_path))
        assert result == ["0"]

    def test_integer_time_dirs(self, tmp_path):
        for d in ("3", "1", "2"):
            (tmp_path / d).mkdir()
        result = get_sorted_times(str(tmp_path))
        assert result == ["1", "2", "3"]


# ── load ──────────────────────────────────────────────────────────────────────

class TestLoad:
    """Tests for load – reads whitespace-delimited files with leading # comments."""

    def _write_dat(self, path, content):
        """Write content to a file and return its path as a string."""
        path.write_text(content)
        return str(path)

    def test_reads_two_column_file(self, tmp_path):
        content = (
            "# Time value\n"
            "# Time value\n"
            "0.1 1.0\n"
            "0.2 2.0\n"
        )
        filepath = self._write_dat(tmp_path / "data.dat", content)
        df = load(filepath)
        assert list(df.columns) == ["Time", "value"]
        assert len(df) == 2
        assert df["Time"].tolist() == pytest.approx([0.1, 0.2])
        assert df["value"].tolist() == pytest.approx([1.0, 2.0])

    def test_multiple_comment_lines(self, tmp_path):
        content = (
            "# comment 1\n"
            "# comment 2\n"
            "# Time a b\n"
            "1.0 10.0 20.0\n"
            "2.0 11.0 21.0\n"
        )
        filepath = self._write_dat(tmp_path / "data.dat", content)
        df = load(filepath)
        assert list(df.columns) == ["Time", "a", "b"]
        assert len(df) == 2

    def test_returns_dataframe(self, tmp_path):
        content = (
            "# T v\n"
            "0 0\n"
        )
        filepath = self._write_dat(tmp_path / "data.dat", content)
        result = load(filepath)
        assert isinstance(result, pd.DataFrame)

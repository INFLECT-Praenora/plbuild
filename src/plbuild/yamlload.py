"""YAML loading that keeps line numbers, so every input error can name its file and line."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


class ConfigError(Exception):
    """An input file is missing or malformed. The message starts with `file:line:`."""

    def __init__(self, path: Path | str, line: int | None, message: str):
        self.path = str(path)
        self.line = line
        self.message = message
        where = f"{self.path}:{line}" if line is not None else self.path
        super().__init__(f"{where}: {message}")


class LineDict(dict):
    """A dict that remembers the line (1-based) of each key and of itself."""

    line: int = 0
    key_lines: dict

    def line_of(self, key: Any) -> int:
        return self.key_lines.get(key, self.line)


class LineList(list):
    line: int = 0
    item_lines: list

    def line_of(self, i: int) -> int:
        return self.item_lines[i] if i < len(self.item_lines) else self.line


class _Loader(yaml.SafeLoader):
    pass


def _construct_mapping(loader: _Loader, node: yaml.MappingNode) -> LineDict:
    loader.flatten_mapping(node)
    d = LineDict()
    d.line = node.start_mark.line + 1
    d.key_lines = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=True)
        try:
            hash(key)
        except TypeError:
            raise ConfigError(loader.name, key_node.start_mark.line + 1, "a mapping key must be a scalar")
        if key in d:
            raise ConfigError(loader.name, key_node.start_mark.line + 1, f"duplicate key '{key}'")
        d[key] = loader.construct_object(value_node, deep=True)
        d.key_lines[key] = key_node.start_mark.line + 1
    return d


def _construct_sequence(loader: _Loader, node: yaml.SequenceNode) -> LineList:
    seq = LineList(loader.construct_object(child, deep=True) for child in node.value)
    seq.line = node.start_mark.line + 1
    seq.item_lines = [child.start_mark.line + 1 for child in node.value]
    return seq


_Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)
_Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_SEQUENCE_TAG, _construct_sequence)


def load_yaml(path: Path | str) -> Any:
    """Load a YAML file. Syntax errors become ConfigError with the file and line."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise ConfigError(path, None, "file not found")
    loader = _Loader(text)
    loader.name = str(path)
    try:
        data = loader.get_single_data()
    except yaml.MarkedYAMLError as e:
        mark = e.problem_mark or e.context_mark
        line = mark.line + 1 if mark is not None else None
        raise ConfigError(path, line, f"YAML syntax error: {e.problem or e.context}")
    finally:
        loader.dispose()
    return data

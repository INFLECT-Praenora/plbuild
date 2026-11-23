"""Output files (spec section 11): the library, the full candidate pool, the report and the manifest.

`curated.ant` is the library EnsGen loads: the drawn templates and curated processes, each with an explanatory comment
block so people can read and edit it. It holds no build date, so the same inputs give a byte-identical file; the date
is in manifest.json.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import html
import json
from pathlib import Path

from . import __version__
from .extract import ModelResult
from .process import Process, entry_text
from .select import Replacement


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _library_text(title: str, entries: list[Process], header: list[str]) -> str:
    templates = sorted((p for p in entries if p.kind != "curated"), key=lambda p: p.pid)
    curated = sorted((p for p in entries if p.kind == "curated"), key=lambda p: p.pid)
    bar = "// " + "=" * 76
    out = [bar, f"// {title}", f"// Written by plbuild {__version__}.", "//"]
    out += [f"// {h}" if h else "//" for h in header]
    out += ["//",
            "// How to edit: change values and bounds, delete entries or add new ones by hand. Keep the",
            "// `// @...` lines and the statements below them; the comment blocks above each entry are",
            "// only explanation. Templates use placeholder species (S, P, X, Y, A, ...) that EnsGen binds",
            "// to real species; curated processes use real species named <entity>_<compartment>.",
            bar, ""]
    out += [bar, f"// Templates ({len(templates)})", bar, ""]
    for p in templates:
        out += [entry_text(p)]
    out += [bar, f"// Curated processes ({len(curated)})", bar, ""]
    for p in curated:
        out += [entry_text(p)]
    return "\n".join(out)


def write_outputs(out_dir: Path, *, results: list[ModelResult], pool: list[Process], drawn: list[Process],
                  gaps: list[Replacement], size: int | None, seed: int, time_unit: str,
                  input_hashes: dict[str, str]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs_hash = sha256(json.dumps(input_hashes, sort_keys=True).encode())[:16]
    n_cur = sum(p.kind == "curated" for p in pool)
    header = [
        f"Library time unit: {time_unit}",
        f"Inputs hash: {inputs_hash} (input files and SBML; details in manifest.json)",
        f"Curated processes from: {', '.join(r.model_id for r in results)}",
        f"Drawn: {len(drawn)} of {len(pool)} candidates ({n_cur} curated, {len(pool) - n_cur} templates), seed {seed}"
        if size is not None else f"All {len(pool)} candidates ({n_cur} curated, {len(pool) - n_cur} templates)",
    ]
    if gaps:
        header += ["", "Curated processes left out by the draw, and the drawn templates that can stand in for them:"]
        for g in gaps:
            header.append(f"  {g.pid}  {g.title}: " + (g.replaced_by or "NO REPLACEMENT in this draw"
                                                      + (f" (would need {g.candidates[0]})" if g.candidates else "")))
    library = _library_text("EnsGen process library", drawn, header)
    pool_text = _library_text("EnsGen candidate pool (every process the library is drawn from; EnsGen does not load "
                              "this file)", pool, header[:3])
    paths = {"library": out_dir / "curated.ant", "pool": out_dir / "pool.ant",
             "report_json": out_dir / "report.json", "report_html": out_dir / "report.html",
             "manifest": out_dir / "manifest.json"}
    paths["library"].write_bytes(library.encode("utf-8"))
    paths["pool"].write_bytes(pool_text.encode("utf-8"))

    report = _report(results, pool, drawn, gaps, size, seed)
    paths["report_json"].write_text(json.dumps(report, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    paths["report_html"].write_text(_report_html(report), encoding="utf-8")

    manifest = {
        "tool": "plbuild", "version": __version__,
        "build_date": _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat(),
        "models": [{"id": r.model_id, "name": r.name, "sbml": r.sbml_path.name,
                    "sbml_sha256": input_hashes.get(f"sbml:{r.model_id}")} for r in results],
        "inputs_sha256": input_hashes, "inputs_hash": inputs_hash,
        "selection": {"size": size, "seed": seed, "pool": len(pool), "drawn": len(drawn)},
        "outputs_sha256": {"curated.ant": sha256(library.encode("utf-8")), "pool.ant": sha256(pool_text.encode("utf-8"))},
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return paths


def _report(results, pool, drawn, gaps, size, seed) -> dict:
    drawn_ids = {p.pid for p in drawn}
    models = []
    for r in results:
        models.append({
            "id": r.model_id, "name": r.name,
            "reactions_found": r.n_reactions, "kept": len(r.processes), "excluded": len(r.excluded),
            "in_library": sum(p.pid in drawn_ids for p in r.processes),
            "excluded_reactions": [{"reaction": e.reaction, "code": e.code, "detail": e.detail} for e in r.excluded],
            "unresolved": r.unresolved,
            "conversions": r.conversions,
        })
    models.sort(key=lambda m: -m["excluded"])
    return {
        "totals": {"models": len(results), "reactions_found": sum(m["reactions_found"] for m in models),
                   "kept": sum(m["kept"] for m in models), "excluded": sum(m["excluded"] for m in models),
                   "pool": len(pool), "in_library": len(drawn),
                   "dropped_curated": len(gaps), "gaps": sum(g.replaced_by is None for g in gaps)},
        "selection": {
            "size": size, "seed": seed,
            "library": [{"id": p.pid, "kind": p.kind, "title": p.title} for p in sorted(drawn, key=lambda p: p.pid)],
            "dropped_curated": [{"id": g.pid, "title": g.title, "source": g.source, "replaced_by": g.replaced_by,
                                 "possible_replacements": g.candidates} for g in gaps],
        },
        "models": models,
    }


def _table(headers: list[str], rows: list[list]) -> str:
    h = "".join(f"<th>{html.escape(x)}</th>" for x in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table><tr>{h}</tr>{body}</table>"


def _report_html(rep: dict) -> str:
    t = rep["totals"]
    parts = ["<!DOCTYPE html><html><head><meta charset='utf-8'><title>plbuild report</title><style>",
             "body{font-family:sans-serif;margin:2em}table{border-collapse:collapse;margin:1em 0}",
             "td,th{border:1px solid #ccc;padding:4px 8px;text-align:left}th{background:#eee}.gap{color:#b00}",
             "</style></head><body><h1>plbuild report</h1>",
             _table(list(t.keys()), [list(t.values())]),
             "<h2>Curated processes left out of the library</h2>",
             _table(["id", "process", "source", "replaced by drawn templates"],
                    [[g["id"], g["title"], g["source"], g["replaced_by"] or "NO REPLACEMENT"]
                     for g in rep["selection"]["dropped_curated"]]),
             "<h2>Library entries</h2>",
             _table(["id", "kind", "title"], [[e["id"], e["kind"], e["title"]] for e in rep["selection"]["library"]])]
    for m in rep["models"]:
        parts.append(f"<h2>{html.escape(m['id'])}: {html.escape(m['name'])}</h2>")
        parts.append(_table(["found", "kept", "excluded", "in library"],
                            [[m["reactions_found"], m["kept"], m["excluded"], m["in_library"]]]))
        if m["excluded_reactions"]:
            parts.append(_table(["reaction", "reason", "detail"],
                                [[e["reaction"], e["code"], e["detail"]] for e in m["excluded_reactions"]]))
        if m["unresolved"]:
            parts.append(_table(["kind", "id", "name", "reason"],
                                [[u["kind"], u["id"], u["name"], u["reason"]] for u in m["unresolved"]]))
    parts.append("</body></html>\n")
    return "".join(parts)

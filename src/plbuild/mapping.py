"""Species and compartment mapping (spec sections 6 and 7): override, then annotation, then normalized name."""
from __future__ import annotations

import re

import libsbml

from .config import CompartmentSpec, Entity, ModelOverrides

GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "κ": "kappa", "λ": "lambda",
         "μ": "mu", "ν": "nu", "σ": "sigma", "τ": "tau", "ω": "omega"}


def normalize(text: str) -> str:
    """Lower case; spaces, hyphens and underscores removed; Greek letters spelled out."""
    t = text.lower()
    for g, name in GREEK.items():
        t = t.replace(g, name)
    return re.sub(r"[\s\-_]", "", t)


def annotations(sbase: libsbml.SBase) -> dict[str, list[str]]:
    """Biological qualifier -> resource URIs, e.g. {'is': [...], 'isVersionOf': [...], 'hasPart': [...]}."""
    out: dict[str, list[str]] = {}
    for cv in sbase.getCVTerms() or []:
        if cv.getQualifierType() != libsbml.BIOLOGICAL_QUALIFIER:
            continue
        q = libsbml.BiolQualifierType_toString(cv.getBiologicalQualifierType())
        out.setdefault(q, []).extend(cv.getResourceURI(i) for i in range(cv.getNumResources()))
    return out


def _uri_key(uri: str) -> str:
    """identifiers.org URIs come in several spellings; compare on the last two path parts."""
    u = uri.rstrip("/").replace("urn:miriam:", "").replace(":", "/")
    parts = [p for p in u.split("/") if p]
    return "/".join(parts[-2:]).lower()


class Mapper:
    def __init__(self, vocabulary: dict[str, Entity], compartments: dict[str, CompartmentSpec]):
        self.vocab = vocabulary
        self.comps = compartments
        self._ent_by_uri: dict[str, set[str]] = {}
        self._ent_by_name: dict[str, set[str]] = {}
        for e in vocabulary.values():
            for uri in e.identifiers:
                self._ent_by_uri.setdefault(_uri_key(uri), set()).add(e.name)
            for n in [e.name, *e.synonyms]:
                self._ent_by_name.setdefault(normalize(n), set()).add(e.name)
        self._comp_by_uri: dict[str, set[str]] = {}
        self._comp_by_name: dict[str, set[str]] = {}
        for c in compartments.values():
            for uri in c.identifiers:
                self._comp_by_uri.setdefault(_uri_key(uri), set()).add(c.name)
            for n in [c.name, *c.synonyms]:
                self._comp_by_name.setdefault(normalize(n), set()).add(c.name)

    # ------------------------------------------------------------------ species

    def species(self, sp: libsbml.Species, ov: ModelOverrides | None) -> tuple[str | None, str]:
        """Return (entity name or None, how it was decided). None means unresolved; the reason says why."""
        sid = sp.getId()
        if ov and sid in ov.species:
            return ov.species[sid], "override"
        ann = annotations(sp)
        is_hits = self._hits(self._ent_by_uri, ann.get("is", []))
        if len(is_hits) == 1:
            return is_hits.pop(), "annotation (is)"
        if len(is_hits) > 1:
            return None, f"ambiguous annotation: {sorted(is_hits)}"
        ver_hits = self._hits(self._ent_by_uri, ann.get("isVersionOf", []) + ann.get("hasVersion", []))
        if len(ver_hits) == 1:
            return ver_hits.pop(), "annotation (isVersionOf/hasVersion)"
        parts = ann.get("hasPart", [])
        if len(parts) >= 2:
            names = []
            for uri in parts:
                hit = self._hits(self._ent_by_uri, [uri])
                if len(hit) != 1:
                    names = []
                    break
                names.append(hit.pop())
            if names:
                cname = "_".join(sorted(names))
                if cname in self.vocab:
                    return cname, "complex (hasPart)"
        for text in (sp.getName(), sid):
            if not text:
                continue
            hits = self._ent_by_name.get(normalize(text), set())
            if len(hits) == 1:
                return next(iter(hits)), f"name '{text}'"
            if len(hits) > 1:
                return None, f"ambiguous name '{text}': {sorted(hits)}"
        return None, "no vocabulary entry matches"

    # ------------------------------------------------------------------ compartments

    def compartment(self, comp: libsbml.Compartment, ov: ModelOverrides | None, n_compartments: int,
                    default: str | None) -> tuple[str | None, str]:
        cid = comp.getId()
        if ov and cid in ov.compartments:
            return ov.compartments[cid], "override"
        ann = annotations(comp)
        hits = self._hits(self._comp_by_uri, ann.get("is", []) + ann.get("isVersionOf", []))
        if len(hits) == 1:
            return hits.pop(), "annotation"
        for text in (comp.getName(), cid):
            if text:
                hits = self._comp_by_name.get(normalize(text), set())
                if len(hits) == 1:
                    return next(iter(hits)), f"name '{text}'"
        if n_compartments == 1 and default:
            return default, "model default_compartment"
        return None, "no configured compartment matches"

    @staticmethod
    def _hits(index: dict[str, set[str]], uris: list[str]) -> set[str]:
        out: set[str] = set()
        for uri in uris:
            out |= index.get(_uri_key(uri), set())
        return out

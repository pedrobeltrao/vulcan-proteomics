"""Cross-linker definitions as data (lit/04 A1: mass, specificity, cleavage stubs, mono-/loop-link masses).

Masses are monoisotopic, Da, as added to the unmodified peptide(s):
  crosslink:  M = m(alpha) + m(beta) + mass
  loop-link:  M = m(peptide) + mass
  mono-link:  M = m(peptide) + mass + quench (hydrolysed: H2O; amidated: NH3)
  cleavable:  after MS2 cleavage each peptide carries one stub; stubs[0] + stubs_complement = mass
DSSO stubs (Kao et al. MCP 2011): alkene A = +54.01056 (C3H2O), unsaturated thiol T = +85.98264 (C3H2OS),
sulfenic acid S = +103.99320 (= T + H2O); A + S = 158.00376 = intact DSSO; T - A = 31.97208 (doublet spacing).
"""
from __future__ import annotations

from dataclasses import dataclass, field

H2O = 18.010565
NH3 = 17.026549


@dataclass(frozen=True)
class Linker:
    name: str
    mass: float                        # added mass of the intact cross-link
    sites: str = "KSTY^"               # linkable residues; '^' = peptide N-terminus
    cleavable: bool = False
    stubs: tuple = ()                  # stub masses observed on each released peptide (cleavable only)
    doublet: tuple = ()                # (light stub, heavy stub) forming the MS2 doublet
    mono_quench: tuple = (H2O, NH3)    # mono-link end groups

    @property
    def mono_masses(self) -> tuple:
        return tuple(self.mass + q for q in self.mono_quench)

    @property
    def loop_mass(self) -> float:
        return self.mass

    @property
    def doublet_delta(self) -> float:
        return self.doublet[1] - self.doublet[0] if self.doublet else 0.0


DSS = Linker("DSS", 138.068080, sites="KSTY^")
BS3 = Linker("BS3", 138.068080, sites="KSTY^")
DSSO = Linker("DSSO", 158.003765, sites="KSTY^", cleavable=True,
              stubs=(54.010565, 85.982635, 103.993200), doublet=(54.010565, 85.982635))

LINKERS = {x.name: x for x in (DSS, BS3, DSSO)}

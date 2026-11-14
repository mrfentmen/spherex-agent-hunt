"""spherex_hunt: a fully agentic hunter for moving objects in NASA SPHEREx data.

NASA Space Apps Challenge 2026, Challenge 12 (Planet X and SPHEREx).
Team: Abdelrahman El Sayad (solo), New York local event.

The pipeline is built around one rule: nothing is called a discovery until
it survives every proof gate. Candidates that fail a gate are killed.
Candidates with an unverified gate stay candidates, never discoveries.
"""

__version__ = "0.1.0"

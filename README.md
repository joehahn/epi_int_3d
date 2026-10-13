# epi_int_3d

*Work in progress — currently in the planning stage.*

epi_int_3d will extend the 2D epi_int_lite streamline N-body integrator for gravitating
planetary rings into the vertical direction, and use it to simulate spiral bending waves.
The project also serves as a demo of using an AI coding agent (Claude Code) to accelerate
research: reviewing papers, deriving equations, porting and extending code, testing,
running simulations, and drafting the journal paper.

- Plan: [PLAN.md](PLAN.md)
- Prompts that drive the project: [prompts/](prompts/)
- Equation verification status: [notes/equations.md](notes/equations.md)
- Reference papers: [papers/](papers/)

Built on:
- epi_int_lite (Python): https://github.com/joehahn/epi_int_lite — Hahn et al. 2025, ApJ, 986, 220
- epi_int (IDL): https://github.com/joehahn/epi_int — Hahn & Spitale 2013, ApJ, 772, 122

Joe Hahn, jmh.datasciences@gmail.com

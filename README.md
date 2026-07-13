# STIRAP Interactive Learning Tool

An interactive, in-browser visualization of a generalized fractional Stimulated Raman Adiabatic Passage (STIRAP) scheme for controlling a coherent superposition of quantum states. Built as a learning tool for exploring coherent population transfer through live, adjustable plots rather than static figures.

**[Live demo](http://stirap.ngolubev.com/)**

## What it covers

The interface is organized into two main sections:

- **Control Panel** - adjust the energy levels of the three-level system, the initial/final populations, the Pump and Stokes laser pulse parameters (detuning, duration, timing), and the time-propagation window.
- **Plot** - updates in real time as parameters change, showing:
  - **Envelope and Pulse Functions** - the envelopes and electric fields of the Pump and Stokes pulses.
  - **Angle** - the time evolution of the mixing angle that governs the adiabatic population transfer.
  - **Energy-Level Diagram** - the three-level system and its laser couplings.
  - **Populations** - the time-dependent populations of the three levels, compared against the exact numerical solution alongside the rotating-wave approximation (RWA).

Each panel and plot has an accompanying info icon with a plain-language explanation of what's being shown.

## Theory

A detailed description of the generalized fractional STIRAP protocol, its theoretical foundations, and its application to coherent control of quantum dynamics in atomic systems can be found in:

- [M. A. Alarcón, K. Hauser, and N. V. Golubev; Phys. Rev. A **113**, 013112 (2026)](https://doi.org/10.1103/2k75-h2pm)

## Running locally

This is a static site with no build step. Serve the directory with any static file server and open it in a browser, for example:

```bash
python3 -m http.server
```

Then visit `http://localhost:8000`.

## Tech stack

- Vanilla HTML/CSS/JS (no framework)
- [Plotly.js](https://plotly.com/javascript/) for plotting
- [MathJax](https://www.mathjax.org/) for math rendering
- [math.js](https://mathjs.org/) and a bundled `numeric.js` for the numerical time propagation

## Developers

Brayton Bosuku, Miguel Alarcón, and Nikolay Golubev

## License

[GPL-3.0](LICENSE)

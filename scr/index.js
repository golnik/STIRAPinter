// Conversion of units from SI to atomic units
const fsperau = (2.488843e-17)/(1.0e-15);
const auI = 3.50944e16;
const evperAU = 27.2114079527e0;


// Shared spacing so every plot's whitespace around the drawing area matches.
// Top margin leaves room for the horizontal legend row above the traces;
// left margin leaves room for the wider y-axis title gap below.
const PLOT_MARGIN = { l: 70, r: 15, b: 55, t: 52, pad: 10 };

// Draws the legend as a horizontal row above the plotting area, in its
// reserved top margin, instead of overlapping the data.
const LEGEND_TOP = { orientation: 'h', x: 0.5, xanchor: 'center', y: 1, yanchor: 'bottom', font: { size: 18 } };

// Fixed gap between each axis and its title, so the gap doesn't vary with
// tick label width (Plotly's automatic standoff differs per plot otherwise).
// The y-axis gets a larger gap since its tick labels sit between the axis
// line and the title.
const AXIS_TITLE_STANDOFF_X = 15;
const AXIS_TITLE_STANDOFF_Y = 30;

const defaultValues = {
    C0: 0.25,
    CF: 0,
    t_s: 12,
    detuning0: 0,
    duration: 12,
    time0: 0,
    timef: 100,
    toggle_curves: true,
    E1: 0,
    E2: 5.266919,
    E3: 1.306423
};


function syncPositionSlider(durationValue) {
    const durationVal = parseFloat(durationValue);
    const min = -2 * durationVal;
    const max = 2 * durationVal;


    positionSlider.min = min;
    positionSlider.max = max;


    const newDelay = durationVal;
    positionSlider.value = newDelay;
    positionOutput.textContent = newDelay.toFixed(1);
}


//Returns all the values back to their defaults
document.getElementById("default").onclick = function resetToDefaults() {
    Object.entries(defaultValues).forEach(([key, value]) => {
        const element = document.getElementById(key);
        if (!element) return;


        if (element.type === "checkbox") {
            element.checked = value;
        } else {
            element.value = value;
        }
    });


    document.getElementById("time0").value = defaultValues.time0;
    document.getElementById("timef").value = defaultValues.timef;
    document.getElementById("detuning0_Value").textContent = defaultValues.detuning0.toFixed(1);
    document.getElementById("duration_Value").textContent = defaultValues.duration.toFixed(1);
    document.getElementById("C0_Value").textContent = defaultValues.C0.toFixed(2);
    document.getElementById("CF_Value").textContent = defaultValues.CF.toFixed(2);
    document.getElementById("t_s_Value").textContent = defaultValues.t_s.toFixed(1);


    syncPositionSlider(defaultValues.duration);
    updatePlots();
};


// UI Elements mapping
const sliders = {
    C0: document.getElementById("C0"),
    CF: document.getElementById("CF"),
    t_s: document.getElementById("t_s"),
    detuning0: document.getElementById("detuning0"),
    duration: document.getElementById("duration")
};


const outputs = {
    C0: document.getElementById("C0_Value"),
    CF: document.getElementById("CF_Value"),
    t_s: document.getElementById("t_s_Value"),
    detuning0: document.getElementById("detuning0_Value"),
    duration: document.getElementById("duration_Value")
};


const TWO_DECIMAL_KEYS = new Set(["C0", "CF"]);
const ONE_DECIMAL_KEYS = new Set(["detuning0", "duration", "t_s"]);

// Update displayed values dynamically
Object.keys(sliders).forEach(key => {
    sliders[key].addEventListener("input", function() {
        if (TWO_DECIMAL_KEYS.has(key)) {
            outputs[key].textContent = parseFloat(this.value).toFixed(2);
        } else if (ONE_DECIMAL_KEYS.has(key)) {
            outputs[key].textContent = parseFloat(this.value).toFixed(1);
        } else {
            outputs[key].textContent = this.value;
        }
    });
});


// Debounce function to prevent the heavy ODE solver from freezing the browser during slider drag
function debounce(func, wait) {
    let timeout;
    return function(...args) {
        clearTimeout(timeout);
        timeout = setTimeout(() => func.apply(this, args), wait);
    };
}
const durationSlider = document.getElementById('duration');
const positionSlider = document.getElementById("t_s");
const positionOutput = document.getElementById("t_s_Value");
durationSlider.addEventListener("input", function() {
    syncPositionSlider(this.value);
});


// Keeps t0 < tf so the ODE integrator (numeric.dopri only steps forward)
// never gets called with an empty or inverted time window, which otherwise
// silently collapses every plot to a single degenerate point.
const time0Input = document.getElementById("time0");
const timefInput = document.getElementById("timef");
function syncTimeWindow(changedId) {
    const t0 = parseFloat(time0Input.value);
    const tf = parseFloat(timefInput.value);
    if (!isFinite(t0) || !isFinite(tf) || t0 < tf) return;

    if (changedId === "time0") {
        timefInput.value = t0 + 1;
    } else {
        time0Input.value = Math.max(0, tf - 1);
    }
}
time0Input.addEventListener("input", () => syncTimeWindow("time0"));
timefInput.addEventListener("input", () => syncTimeWindow("timef"));


// Core calculation and plotting function
function updatePlots() {

    //Creates all the consts from the user input, and treats them as floats
    const time0 = parseFloat(document.getElementById("time0").value);
    const timef = parseFloat(document.getElementById("timef").value);
    const delta_0 = parseFloat(document.getElementById('detuning0').value);
    const gboth = parseFloat(document.getElementById('duration').value);
   
    const tp_inp = 50.0;
    const ts_inp = tp_inp + parseFloat(document.getElementById("t_s").value);
    const C_0 = parseFloat(document.getElementById("C0").value);
    const C_F = parseFloat(document.getElementById("CF").value);
    const showCurves = document.getElementById("toggle_curves").checked;
    const E1_inp = parseFloat(document.getElementById("E1").value);
    const E2_inp = parseFloat(document.getElementById("E2").value);
    const E3_inp = parseFloat(document.getElementById("E3").value);

    let test = Gaussion_Values(tp_inp, ts_inp, C_0, C_F, delta_0, gboth, E1_inp, E2_inp, E3_inp);
    let Delta = delta_creation(test.E1, test.E2, test.E3, test.wp, test.ws);


    sliders.t_s.min = -2 * test.gp * fsperau;
    sliders.t_s.max = 2 * test.gp * fsperau;
    test.delta = Delta;


    test.t0 = time0 / fsperau;
    test.tf = timef / fsperau;


    const t_values = [];
    const stagent_p = [];
    const stagent_s = [];
    const no_frequency_p = [];
    const no_frequency_s = [];
    const angle_va = [];
    const beta_va = [];
    const alpha_va = [];

    const steps = 2000;
    const dt = (test.tf - test.t0) / steps;


    // Loop through the wanted time from & #of steps to get x & y values
    for (let i = 0; i <= steps; i++) {
        const t = test.t0 + i * dt;
        t_values.push(t * fsperau);


        let stagent_p0 = (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta));
        stagent_p.push(stagent_p0 * Math.cos(test.wp * t));
        no_frequency_p.push(stagent_p0);


        let stagent_s0 = (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta));
        stagent_s.push(stagent_s0 * Math.cos(test.ws * t));
        no_frequency_s.push(stagent_s0);


        angle_va.push(Math.atan(-stagent_p0 / stagent_s0) / Math.PI);
        alpha_va.push(Math.acos(Math.sqrt(C_0))/Math.PI);
        beta_va.push(Math.acos(Math.sqrt(C_F))/Math.PI);
    }

    // Plot I Configurations
    const Gau_p = {
        x: t_values, y: stagent_p, name: 'Ω<sub>P</sub>(t)cos(ω<sub>P</sub>t)',
        mode: 'lines', line: {color: 'blue', width: 1.5},
        visible: showCurves
    };
    const Gau_s = {
        x: t_values, y: stagent_s, name: 'Ω<sub>S</sub>(t)cos(ω<sub>S</sub>t)',
        mode: 'lines', line: {color: 'red', width: 1.5},
        visible: showCurves
    };
    const Gau_P_Stagnent = {
        x: t_values, y: no_frequency_p, name: 'Ω<sub>P</sub>(t)',
        mode: 'lines', line: {color: 'blue', width: 2.5}
    };
    const Gau_S_Stagnent = {
        x: t_values, y: no_frequency_s, name: 'Ω<sub>S</sub>(t)',
        mode: 'lines', line: {color: 'red', width: 2.5}
    };


    const layout = {
        font: { family: "'STIX Two Text', serif",
            size: 14},
        xaxis: {title: {text: 'Time (fs)', standoff: AXIS_TITLE_STANDOFF_X}},
        yaxis: {title: {text: 'Amplitude', standoff: AXIS_TITLE_STANDOFF_Y}, range: [-2*test.Op,2*test.Op]},
        showlegend: true,
        legend: LEGEND_TOP,
        margin: PLOT_MARGIN
    };
    Plotly.react('plot', [Gau_p, Gau_s, Gau_P_Stagnent, Gau_S_Stagnent], layout, {responsive: true, displayModeBar: false});


    //Plot 2 Configurations
    const angles = {
        x: t_values,
        y: angle_va,
        name: 'θ(t)',
        mode: 'lines',
        line: {color: 'green', width: 1.5}
    }


    const Beta_Line = {
        x: t_values,
        y: beta_va,
        name: 'β',
        mode: 'lines',
        line: {color: 'grey', width: 1.5, dash: 'dot'}
    }
    const Alpha_Line = {
        x: t_values,
        y: alpha_va,
        name: 'α',
        mode: 'lines',
        line: {color: 'grey', width: 1.5, dash: 'dot'}
    }
   
    // NEW: Updated layoutII to include annotations for the Greek letters
    const layoutII = {
        font: { family: "'STIX Two Text', serif",
            size: 14},
        xaxis: {title: {text: 'Time (fs)', standoff: AXIS_TITLE_STANDOFF_X}},
        yaxis: {title: {text: 'Angle (rad)', standoff: AXIS_TITLE_STANDOFF_Y}, range: [0.0, 0.55]},
        showlegend: false,
        annotations: [
            {
                x: time0 + (timef - time0) * 0.05,         // Near the initial time
                y: Math.acos(Math.sqrt(C_0))/Math.PI,        // Use the current alpha numerical value
                text: '<b>α</b>', // Greek letter alpha
                showarrow: false,
                yshift: 10,    // Shifts the text slightly above the line
                font: {size: 16, color: 'black'}
            },
            {
                x: timef - (timef - time0) * 0.05,         // Near the final time
                y: Math.acos(Math.sqrt(C_F))/Math.PI,        // Use the current beta numerical value
                text: '<b>β</b>', // Greek letter beta
                showarrow: false,
                yshift: 10,    // Shifts the text slightly above the line
                font: {size: 16, color: 'black'}
            }
        ],
        margin: { ...PLOT_MARGIN, t: 15 }
    };


    Plotly.newPlot('plotII', [angles, Beta_Line, Alpha_Line], layoutII, {responsive: true, displayModeBar: false});


    // Plot III Configurations
    // Does the calculations for the populations transfer graphs using data from the system of equations established below
    let res = population_calculations_NRW(test);
    let resII = population_calculations_RWA(test);


    const tt = [], c1 = [], c3 = [], tt_RWA = [], c1_RWA = [], c3_RWA = [];

    // res (NRW) and resII (RWA) are two independent adaptive-step integrations with different number of steps
    for (let i = 0; i < res[0].length; i++) {
        tt.push(res[0][i] * fsperau);
        c1.push(res[1][0][i]**2 + res[1][1][i]**2);
        c3.push(res[1][4][i]**2 + res[1][5][i]**2);
    }
    for (let i = 0; i < resII[0].length; i++) {
        tt_RWA.push(resII[0][i] * fsperau);
        c1_RWA.push(resII[1][0][i]**2 + resII[1][1][i]**2);
        c3_RWA.push(resII[1][4][i]**2 + resII[1][5][i]**2);
    }

       const C1 = {
        x: tt, y: c1, name: '|C<sub>1</sub>(t)|<sup>2</sup>',
        mode: 'lines', line: {color: 'blue', width: 2},
        visible: showCurves
    };
    const C3 = {
        x: tt, y: c3, name: '|C<sub>3</sub>(t)|<sup>2</sup>',
        mode: 'lines', line: {color: 'red', width: 2},
        visible: showCurves
    };
    const C1_RWAs = {
        x: tt_RWA, y: c1_RWA, name: '|C<sub>1</sub><sup>RWA</sup>(t)|<sup>2</sup>',
        mode: 'lines', line: {color: 'orange', width: 2, dash: 'dot'}
    };
    const C3_RWAs = {
        x: tt_RWA, y: c3_RWA, name: '|C<sub>3</sub><sup>RWA</sup>(t)|<sup>2</sup>',
        mode: 'lines', line: {color: 'purple', width: 2, dash: 'dot'}
    };
    // Target reference markers (the desired initial/final populations set by
    // the Population Controls sliders), each in its own well-distinguished
    // color so they read as targets rather than competing with the computed
    // dynamics (which use blue/red/orange/purple).
    const INITIAL_COLOR = '#16a34a';  // Green
    const FINAL_COLOR = '#db2777';    // Magenta/pink
    const t_init = tt[0];
    const t_final = tt[tt.length - 1];
    // cliponaxis:false keeps these fully visible even though they sit exactly
    // on the plot's tight time bounds, where they'd otherwise get half-clipped.
    const C_init = {
        x: [t_init], y: [C_0], name: 'Initial',
        mode: 'markers', marker: {color: INITIAL_COLOR, size: 9, symbol: 'circle-open', line: {width: 2}},
        cliponaxis: false
    };
    const C_final = {
        x: [t_final], y: [C_F], name: 'Desired Final',
        mode: 'markers', marker: {color: FINAL_COLOR, size: 9, symbol: 'circle-open', line: {width: 2}},
        cliponaxis: false
    };


    const layoutIII = {
        font: {family: "'STIX Two Text', serif",
            size: 14},
        xaxis: {title: {text: 'Time (fs)', standoff: AXIS_TITLE_STANDOFF_X}, range: [t_init, t_final]},
        yaxis: {title: {text: 'Probability', standoff: AXIS_TITLE_STANDOFF_Y}, range: [0, 1.2]},
        showlegend: true,
        legend: LEGEND_TOP,
        margin: PLOT_MARGIN
    };
    Plotly.react('plotIII', [C1, C3, C1_RWAs, C3_RWAs, C_init, C_final], layoutIII, {responsive: true, displayModeBar: false });
    
    //Plot IIII Configurations
    const upperY = test.E2 * evperAU;
    const lowerY1 = test.E1* evperAU;
    const lowerY2 = test.E3* evperAU;


    const freqp = test.wp * evperAU;
    const freqs = test.ws * evperAU;


    const bwp = 4*Math.sqrt(2*Math.log(2.0))/test.gp * evperAU;
    const bws = 4*Math.sqrt(2*Math.log(2.0))/test.gs * evperAU;


    const LEVEL_COLOR = '#4a1420'; // Dark maroon energy-level bars
    const PUMP_COLOR = '#2255dd';  // Blue: |1⟩ -> |3⟩ (matches Pump color elsewhere)
    const STOKES_COLOR = '#e74c3c'; // Red: |2⟩ -> |3⟩ (matches Stokes color elsewhere)

    const trace1 = {
        x: [1, 3],
        y: [lowerY1, lowerY1],
        mode: 'lines+text',
        text: ['', '|1⟩'], // Label on the right side
        textposition: 'middle right',
        textfont: { size: 22 },
        line: { color: LEVEL_COLOR, width: 5 },
        hoverinfo: 'none',
        showlegend: false
    };


    const trace2 = {
        x: [7, 9],
        y: [lowerY2, lowerY2],
        mode: 'lines+text',
        text: ['', '|3⟩'],
        textposition: 'middle right',
        textfont: { size: 22 },
        line: { color: LEVEL_COLOR, width: 5 },
        hoverinfo: 'none',
        showlegend: false
    };


    const trace3 = {
        x: [4, 6],
        y: [upperY, upperY],
        mode: 'lines+text',
        text: ['', '|2⟩'],
        textposition: 'middle right',
        textfont: { size: 22 },
        line: { color: LEVEL_COLOR, width: 5 },
        hoverinfo: 'none',
        showlegend: false
    };


    // Bounds for the arrow-style energy axis drawn to the left of the diagram.
    // axisX matches the xaxis range's left edge, where Plotly's side:'left'
    // ticks are anchored, so the custom arrow lines up with them.
    // Padding is based on the level spacing, not the laser bandwidth (bwp) -
    // otherwise the Pulse Duration slider would rescale the whole diagram
    // even though it doesn't actually move the energy levels.
    const axisX = -0.5;
    const levelMin = Math.min(lowerY1, lowerY2, upperY);
    const levelMax = Math.max(lowerY1, lowerY2, upperY);
    const levelSpan = Math.max(levelMax - levelMin, 0.5);
    const axisYMin = levelMin - 0.3*levelSpan;
    const axisYMax = levelMax + 0.45*levelSpan;

    // 2. Define the Layout (Arrows and removing axes)
    const layoutIV = {
        font: { family: "'STIX Two Text', serif", size: 14 },
        // This diagram has no legend row and no visible x-axis, so it doesn't need
        // PLOT_MARGIN's t/b space reserved for those - shrinking them lets the
        // diagram fill the card down to the bottom instead of leaving a gap.
        margin: { ...PLOT_MARGIN, l: 60, pad: 0, t: 36, b: 10 },
        xaxis: {
            visible: false, // Hide the X axis entirely - it's just layout spacing, not physical
            range: [axisX, 10]  // Set a fixed internal coordinate system, with room for the energy axis
        },
        yaxis: {
            visible: true,
            side: 'left',
            showline: false, // The arrow annotation below draws the axis line instead
            zeroline: false,
            showgrid: false, // Default gridlines span the full plot; using custom rule lines below instead
            ticks: 'outside',
            tickwidth: 3,
            tickcolor: 'black',
            ticklen: 8,
            tickfont: { size: 16 },
            tickvals: [lowerY1, lowerY2, upperY], // Exact energies, so ticks line up with the level bars
            tickformat: '.2f', // Display with 2 digits after the decimal point
            range: [axisYMin, axisYMax]
        },
        shapes: [
            // Rule line from the axis to where the |1⟩ level bar starts
            {
                type: 'line',
                x0: axisX, x1: 1, y0: lowerY1, y1: lowerY1,
                line: { color: '#999999', width: 1, dash: 'dot' }
            },
            // Rule line from the axis to where the |3⟩ level bar starts
            {
                type: 'line',
                x0: axisX, x1: 4, y0: upperY, y1: upperY,
                line: { color: '#999999', width: 1, dash: 'dot' }
            },
            // Rule line from the axis to where the |2⟩ level bar starts
            {
                type: 'line',
                x0: axisX, x1: 7, y0: lowerY2, y1: lowerY2,
                line: { color: '#999999', width: 1, dash: 'dot' }
            },
            // Pump Laser Bandwidth - Tail (State 1)
            {
                type: 'rect',
                x0: 1.5, x1: 2.5,      // Centers the box around the tail at X=2
                y0: lowerY1-bwp/2, y1: lowerY1+bwp/2,   // The "height" or bandwidth thickness
                fillcolor: PUMP_COLOR,
                opacity: 0.2,          // Makes it lightly shaded
                line: { width: 0 }     // Removes the hard border
            },
            // Pump Laser Bandwidth - Head (State 3)
            {
                type: 'rect',
                x0: 4.1, x1: 4.9,      // Centers around the head at X=4.5
                y0: lowerY1+freqp-bwp/2, y1: lowerY1+freqp+bwp/2,
                fillcolor: PUMP_COLOR,
                opacity: 0.2,
                line: { width: 0 }
            },
            // Stokes Laser Bandwidth - Tail (State 2)
            {
                type: 'rect',
                x0: 7.5, x1: 8.5,      // Centers around the tail at X=8
                y0: lowerY2-bwp/2, y1: lowerY2+bwp/2,
                fillcolor: STOKES_COLOR,
                opacity: 0.2,
                line: { width: 0 }
            },
            // Stokes Laser Bandwidth - Head (State 3)
            {
                type: 'rect',
                x0: 5.1, x1: 5.9,      // Centers around the head at X=5.5
                y0: lowerY2+freqs-bwp/2, y1: lowerY2+freqs+bwp/2,
                fillcolor: STOKES_COLOR,
                opacity: 0.2,
                line: { width: 0 }
            },
            // Pump Arrow connecting line (the arrowheads are drawn as
            // annotations below - Plotly only guarantees the x/y end of an
            // annotation arrow lands exactly on target, not the ax/ay end,
            // so a single double-headed annotation can't be exact at both
            // ends. Drawing the line as a shape plus two head-only
            // annotations keeps both tips exact.)
            {
                type: 'line',
                x0: 2, y0: lowerY1, x1: 4.5, y1: lowerY1+freqp,
                line: { color: PUMP_COLOR, width: 3 }
            },
            // Stokes Arrow connecting line
            {
                type: 'line',
                x0: 8, y0: lowerY2, x1: 5.5, y1: lowerY2+freqs,
                line: { color: STOKES_COLOR, width: 3 }
            }
        ],
        annotations: [
            // Pump Arrow head (near |3⟩), exact at the bandwidth box center
            {
                ax: 3.25, ay: lowerY1 + freqp/2, // Midpoint - direction reference only
                axref: 'x',
                ayref: 'y',
                x: 4.5,
                y: lowerY1+freqp,
                xref: 'x',
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1,
                arrowwidth: 3,
                arrowcolor: PUMP_COLOR
            },
            // Pump Arrow tail (near |1⟩), exact at the bandwidth box center
            {
                ax: 3.25, ay: lowerY1 + freqp/2,
                axref: 'x',
                ayref: 'y',
                x: 2,
                y: lowerY1,
                xref: 'x',
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1,
                arrowwidth: 3,
                arrowcolor: PUMP_COLOR
            },
            // Stokes Arrow head (near |3⟩), exact at the bandwidth box center
            {
                ax: 6.75, ay: lowerY2 + freqs/2,
                axref: 'x',
                ayref: 'y',
                x: 5.5,
                y: lowerY2+freqs,
                xref: 'x',
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1,
                arrowwidth: 3,
                arrowcolor: STOKES_COLOR
            },
            // Stokes Arrow tail (near |2⟩), exact at the bandwidth box center
            {
                ax: 6.75, ay: lowerY2 + freqs/2,
                axref: 'x',
                ayref: 'y',
                x: 8,
                y: lowerY2,
                xref: 'x',
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1,
                arrowwidth: 3,
                arrowcolor: STOKES_COLOR
            },
            // Pump transition energy label
            {
                x: 3.25, y: lowerY1 + freqp/2,
                xref: 'x', yref: 'y',
                text: freqp.toFixed(2) + ' eV',
                showarrow: false,
                xshift: -40,
                font: { size: 16, color: PUMP_COLOR }
            },
            // Stokes transition energy label
            {
                x: 6.75, y: lowerY2 + freqs/2,
                xref: 'x', yref: 'y',
                text: freqs.toFixed(2) + ' eV',
                showarrow: false,
                xshift: 40,
                font: { size: 16, color: STOKES_COLOR }
            },
            // Arrow-style energy axis: vertical line with an arrowhead at the top
            {
                ax: axisX, ay: axisYMin,
                axref: 'x', ayref: 'y',
                x: axisX, y: axisYMax,
                xref: 'x', yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 0.9,
                arrowwidth: 3,
                arrowcolor: 'black'
            },
            // Axis label centered above the arrow tip
            {
                x: axisX, y: axisYMax,
                xref: 'x', yref: 'y',
                text: 'E [eV]',
                showarrow: false,
                xanchor: 'center',
                yanchor: 'bottom',
                yshift: 8,
                font: { size: 20, color: 'black' }
            }
        ],
        plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background
        paper_bgcolor: 'rgba(0,0,0,0)'
    };


    // 3. Render the Plot
    Plotly.newPlot('energyPlot', [trace1, trace2, trace3], layoutIV, {staticPlot: true, responsive: true});
}




// Math/Physics Helpers
// E1, E2, E3 are the user-adjustable level positions (eV). wp/ws (the pump/Stokes
// laser frequencies) stay resonant with the E2-E1 and E2-E3 spacing plus the
// user's detuning, so moving a level automatically retunes the lasers to match.
function Gaussion_Values(tp, ts, C0, CF, det, gboth, E1, E2, E3) {
    return {
        Os: Math.sqrt((5.803548e11)*(4.33**2)/auI),
        ts: ts/fsperau,
        gs: gboth/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        ws: (E2-E3+det)/evperAU,
        Op: Math.sqrt((5.803548e11)*(4.33**2)/auI),
        tp: tp/fsperau,
        gp: gboth/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        wp: (E2-E1+det)/evperAU,
        E1: E1/evperAU,
        E2: E2/evperAU,
        E3: E3/evperAU,
        mu12: 1,
        mu23: -1,
        alpha: Math.acos(Math.sqrt(C0)),
        beta: Math.acos(Math.sqrt(CF))
    };
}


function envelope(t, t0, alpha, as, ap, mu) {
    let exs = Math.exp(alpha*(t-t0));
    return (1/mu)*alpha*(as-ap)*exs/((1+exs)*Math.sqrt((1-as+(1-ap)*exs)*(as+ap*exs)));
}


function delta_creation(E1, E2, E3, wp, ws) {
    let delta_12 = E1 - E2 + wp;
    let delta_23 = E3 - E2 + ws;
    if (Math.abs(delta_12).toFixed(5) !== Math.abs(delta_23).toFixed(5)) {
        console.warn("The two deltas are not equal, please check the values of the energies and frequencies");
    }
    return delta_12;
}


function Gaussion_Creation_S(t, ts, gs) {
    return Math.exp(-((t-ts)**2/gs**2));
}


function Gaussion_Creation_P(t, tp, gp) {
    return Math.exp(-((t-tp)**2/gp**2));
}


function population_calculations_NRW(test) {
    function F(x, OmegaS, OmegaP) {
        return [
            test.E1 * x[1] - (OmegaP * x[3]),
            -(test.E1 * x[0]) + (OmegaP * x[2]),
            test.E2 * x[3] - (OmegaP * x[1]) - (OmegaS * x[5]),
            -test.E2 * x[2] + (OmegaP * x[0]) + (OmegaS * x[4]),
            test.E3 * x[5] - (OmegaS * x[3]),
            -(test.E3 * x[4]) + (OmegaS * x[2])
        ];
    }
    function f(t, x) {
        let OmegaS = -(test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta)) * Math.cos(test.ws * (t-test.t0));
        let OmegaP = -(test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta)) * Math.cos(test.wp * (t-test.t0));
        return F(x, OmegaS, OmegaP);
    }
    let x0 = [Math.cos(test.alpha), 0, 0, 0, Math.sin(test.alpha), 0];
    let sol = numeric.dopri(test.t0, test.tf, x0, f, 1e-8, 20000);
    return [sol.x, numeric.transpose(sol.y)];
}


//System of Equations created from the Hamiltonion
function population_calculations_RWA(test) {
    function F(x, OmegaS, OmegaP) {
        return [
            test.delta * x[1] - (OmegaP * x[3])/2,
            -(test.delta * x[0]) + (OmegaP * x[2])/2,
            -(OmegaP * x[1])/2 - (OmegaS * x[5])/2,
            (OmegaP * x[0])/2 + (OmegaS * x[4])/2,
            test.delta * x[5] - (OmegaS * x[3])/2,
            -(test.delta * x[4]) + (OmegaS * x[2])/2
        ];
    }
    function f(t, x) {
        let OmegaS = (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta));
        let OmegaP = (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta));
        return F(x, OmegaS, OmegaP);
    }
    let x0 = [Math.cos(test.alpha), 0, 0, 0, Math.sin(test.alpha), 0];
    let sol = numeric.dopri(test.t0, test.tf, x0, f, 1e-8, 2000);
    return [sol.x, numeric.transpose(sol.y)];
}


// Event Listeners for completely dynamic updating
const inputsToWatch = ["C0", "CF", "t_s", "detuning0", "duration", "time0", "timef", "toggle_curves", "E1", "E2", "E3"];


// Apply a 300ms debounce to prevent the UI from freezing when rapidly sliding
const debouncedUpdate = debounce(updatePlots, 300);


inputsToWatch.forEach(id => {
    document.getElementById(id).addEventListener('input', debouncedUpdate);
    document.getElementById(id).addEventListener('change', debouncedUpdate);
});


// Initial Plot Generation on Page Load
window.onload = updatePlots;

// Does the Info Box toggle
function show_info(info) {
    const info_box = document.getElementById(info);
    if (info_box) {
        info_box.style.display = 'block';
    }
}

function hide_info(info) {
    const info_box = document.getElementById(info);
    if (info_box) {
        info_box.style.display = 'none';
    }
}


function toggleHelp(){
    const overlay = document.getElementById('about-modal-overlay');
    overlay.style.display = overlay.style.display === 'flex' ? 'none' : 'flex';
}

// Keeps every plot's pixel size in sync with its (flexbox-driven) card height,
// so charts resize to fit rather than getting clipped by the card.
const plotIds = ['plot', 'plotII', 'plotIII', 'energyPlot'];
function resizeAllPlots() {
    plotIds.forEach(id => {
        const el = document.getElementById(id);
        if (el) Plotly.Plots.resize(el);
    });
}
window.addEventListener('resize', debounce(resizeAllPlots, 150));

// Web fonts finish loading after onload and can shift card heights slightly;
// resize once more when that settles so plots aren't left mis-sized.
if (document.fonts && document.fonts.ready) {
    document.fonts.ready.then(resizeAllPlots);
}

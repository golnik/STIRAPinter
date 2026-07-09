// Conversion of units from SI to atomic units
const fsperau = (2.488843e-17)/(1.0e-15);
const auI = 3.50944e16;
const evperAU = 27.2114079527e0;

// UI Elements mapping
const sliders = {
    C0: document.getElementById("C0"),
    CF: document.getElementById("CF"),
    t_s: document.getElementById("t_s"),
    total_time: document.getElementById("total_time"),
    detuning0: document.getElementById("detuning0"),
    duration: document.getElementById("duration")
};

const outputs = {
    C0: document.getElementById("C0_Value"),
    CF: document.getElementById("CF_Value"),
    t_s: document.getElementById("t_s_Value"),
    total_time: document.getElementById("total_time_Value"),
    detuning0: document.getElementById("detuning0_Value"),
    duration: document.getElementById("duration_Value")
};

// Update displayed values dynamically
Object.keys(sliders).forEach(key => {
    sliders[key].addEventListener("input", function() {
        outputs[key].textContent = this.value;
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
durationSlider = document.getElementById('duration');
const positionSlider = document.getElementById("t_s");
const positionOutput = document.getElementById("t_s_Value");
durationSlider.addEventListener("input", function() {
    const durationVal = parseFloat(this.value);

    // 1. Update the ranges of the position slider 
    // (You can adjust the math here if you want different min/max bounds)
    positionSlider.min = 50.0 - 2 * durationVal;
    positionSlider.max = 50.0 + 2 * durationVal;

    // 2. Force the position value to be exactly 50 + duration
    const newPosition = 50.0 + durationVal;
    
    // 3. Apply the new value to the slider handle and text output
    positionSlider.value = newPosition;
    positionOutput.textContent = newPosition;
});

// Core calculation and plotting function
function updatePlots() {
    const time0 = parseFloat(document.getElementById("time0").value);
    const timef = parseFloat(document.getElementById("timef").value);
    const delta_0 = parseFloat(document.getElementById('detuning0').value);
    const gboth = parseFloat(document.getElementById('duration').value);
    
    const tp_inp = 50.0;
    const ts_inp = parseFloat(document.getElementById("t_s").value);
    const C_0 = parseFloat(document.getElementById("C0").value);
    const C_F = parseFloat(document.getElementById("CF").value);
    const showCurves = document.getElementById("toggle_curves").checked;
    
    let test = Gaussion_Values(tp_inp, ts_inp, C_0, C_F, delta_0, gboth);
    let Delta = delta_creation(test.E1, test.E2, test.E3, test.wp, test.ws);

    sliders.t_s.min = 50.0 - 2 * test.gp * fsperau;
    sliders.t_s.max = 50.0 + 2 * test.gp * fsperau;
    test.delta = Delta;

    test.t0 = time0 / fsperau;
    test.tf = timef / fsperau;

    const t_values = [];
    const t_fq = [];
    const t_lq = [];
    const stagent_p = [];
    const stagent_s = [];
    const no_frequency_p = [];
    const no_frequency_s = [];
    const angle_va = [];
    const beta_va = [];
    const alpha_va = [];
    const init_pop = [];
    const final_pop = [];
    
    const steps = 2000;
    const steps_fq = parseInt(steps/4);
    const steps_lq = parseInt(3/4 * steps);
    const dt = (test.tf - test.t0) / steps;

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

        if(i<=steps_fq){
            t_fq.push(t*fsperau);
            init_pop.push(C_0);
        }

        if(i>steps_lq){
            t_lq.push(t*fsperau);
            final_pop.push(C_F);
        }
    } 

    // Plot I Configurations
    const Gau_p = {
        x: t_values, y: stagent_p, name: 'Pump Frequency Pulse',
        mode: 'lines', line: {color: 'blue', width: 1.5},
        visible: showCurves
    };
    const Gau_s = {
        x: t_values, y: stagent_s, name: 'Stokes Frequency Pulse',
        mode: 'lines', line: {color: 'red', width: 1.5},
        visible: showCurves
    };
    const Gau_P_Stagnent = {
        x: t_values, y: no_frequency_p, name: 'Pump Stagnent Pulse',
        mode: 'lines', line: {color: 'blue', width: 2.5}
    };
    const Gau_S_Stagnent = {
        x: t_values, y: no_frequency_s, name: 'Stokes Stagnent Pulse',
        mode: 'lines', line: {color: 'red', width: 2.5}
    };

    const layout = {
        title: 'Envelope and Pulse Functions', 
        font: { family: "'STIX Two Text', serif",
            size: 14},  
        xaxis: {title: {text:'Time (fs)'}},
        yaxis: {title: {text: 'Amplitude'},range: [-2*test.Op,2*test.Op]},
        margin: { t: 40, l: 50, r: 20, b: 40 }
    };  
    Plotly.react('plot', [Gau_p, Gau_s, Gau_P_Stagnent, Gau_S_Stagnent], layout, {responsive: true}); 

    const angles = {
        x: t_values,
        y: angle_va,
        name: 'Angle',
        mode: 'lines',
        line: {color: 'green', width: 1.5}
    }

    const Beta_Line = {
        x: t_values,
        y: beta_va,
        name: 'Beta',
        mode: 'lines',
        line: {color: 'grey', width: 1.5, dash: 'dot'}
    }
    const Alpha_Line = {
        x: t_values,
        y: alpha_va,
        name: 'Alpha',
        mode: 'lines',
        line: {color: 'grey', width: 1.5, dash: 'dot'}
    }
    
    // NEW: Updated layoutII to include annotations for the Greek letters
    const layoutII = {
        title: {text: 'Adiabatic mixing angle.'},
        font: { family: "'STIX Two Text', serif",
            size: 14},
        xaxis: {title: {text: 'Time (fs)'}},
        yaxis: {title: {text: 'Angle (rad)'}, range: [0.0, 0.55]},
        showlegend: true,
        annotations: [
            {
                x: 25,         // Set x-axis position to 25fs
                y: Math.acos(Math.sqrt(C_0))/Math.PI,        // Use the current alpha numerical value
                text: '<b>α</b>', // Greek letter alpha
                showarrow: false, 
                yshift: 10,    // Shifts the text slightly above the line
                font: {size: 16, color: 'black'}
            },
            {
                x: 75,         // Set x-axis position to 75fs
                y: Math.acos(Math.sqrt(C_F))/Math.PI,        // Use the current beta numerical value
                text: '<b>β</b>', // Greek letter beta
                showarrow: false,
                yshift: 10,    // Shifts the text slightly above the line
                font: {size: 16, color: 'black'}
            }
        ]
    };

    Plotly.newPlot('plotII', [angles, Beta_Line, Alpha_Line], layoutII);

    // Plot III Configurations
    let res = population_calculations_NRW(test);
    let resII = population_calculations_RWA(test);

    const tt = [], c1 = [], c3 = [], tt_RWA = [], c1_RWA = [], c3_RWA = [];
    let ii = res[0].length;

    for(let i = 0; i < ii; i++) {
        tt.push(res[0][i] * fsperau);
        c1.push(res[1][0][i]**2 + res[1][1][i]**2);
        c3.push(res[1][4][i]**2 + res[1][5][i]**2);
        tt_RWA.push(resII[0][i] * fsperau);
        c1_RWA.push(resII[1][0][i]**2 + resII[1][1][i]**2);
        c3_RWA.push(resII[1][4][i]**2 + resII[1][5][i]**2);
    }

    const C1 = {
        x: tt, y: c1, name: 'C1',
        mode: 'lines', line: {color: 'blue', width: 2}
    };
    const C3 = {
        x: tt, y: c3, name: 'C3',
        mode: 'lines', line: {color: 'red', width: 2}
    };
    const C1_RWAs = {
        x: tt_RWA, y: c1_RWA, name: 'C1 (RWA)',
        mode: 'lines', line: {color: 'orange', width: 2, dash: 'dot'},
        visible: showCurves
    };
    const C3_RWAs = {
        x: tt_RWA, y: c3_RWA, name: 'C3 (RWA)',
        mode: 'lines', line: {color: 'purple', width: 2, dash: 'dot'},
        visible: showCurves
    };
    const C_init = {
        x: t_fq, y: init_pop, name: 'C_1 initial',
        mode:'lines', line:{color: 'black', width: 2, dash: 'dashed'}
    };
    const C_final = {
        x: t_lq, y: final_pop, name: 'C_1 final',
        mode:'lines', line:{color: 'black', width:2, dash:'dashed'}
    };

    const layoutIII = {
        title: 'Populations',
        font: {family: "'STIX Two Text', serif",
            size: 14},
        xaxis: {title: {text:'Time (fs)'}},
        yaxis: {title: {text:'Probability'}, range: [0, 1.2]},
        annotations: [
            {
                x: 25,         // Set x-axis position to 25fs
                y: C_0,        // Use the current alpha numerical value
                text: '<b>Initial</b>', // Greek letter alpha
                showarrow: false, 
                yshift: 10,    // Shifts the text slightly above the line
                font: {size: 16, color: 'black'}
            },
            {
                x: 75,         // Set x-axis position to 75fs
                y: C_F,        // Use the current beta numerical value
                text: '<b>Final</b>', // Greek letter beta
                showarrow: false,
                yshift: 10,    // Shifts the text slightly above the line
                font: {size: 16, color: 'black'}
            }
        ],
        margin: { t: 40, l: 50, r: 20, b: 40 }
    };   
    Plotly.react('plotIII', [C1, C3, C1_RWAs, C3_RWAs,C_init,C_final], layoutIII, {responsive: true});


    // Define Y-coordinates for the energy levels
    const upperY = test.E2 * evperAU;
    const lowerY1 = test.E1* evperAU;
    const lowerY2 = test.E3* evperAU;

    const freqp = test.wp * evperAU;
    const freqs = test.ws * evperAU;

    const bwp = 4*Math.sqrt(2*Math.log(2.0))/test.gp * evperAU;
    const bws = 4*Math.sqrt(2*Math.log(2.0))/test.gs * evperAU;

    const trace1 = {
        x: [1, 3], 
        y: [lowerY1, lowerY1],
        mode: 'lines+text',
        text: ['', '|1⟩'], // Label on the right side
        textposition: 'middle right',
        textfont: { size: 18 },
        line: { color: '#34495e', width: 4 },
        hoverinfo: 'none',
        showlegend: false
    };

    const trace2 = {
        x: [7, 9], 
        y: [lowerY2, lowerY2],
        mode: 'lines+text',
        text: ['', '|2⟩'],
        textposition: 'middle right',
        textfont: { size: 18 },
        line: { color: '#34495e', width: 4 },
        hoverinfo: 'none',
        showlegend: false
    };

    const trace3 = {
        x: [4, 6], 
        y: [upperY, upperY],
        mode: 'lines+text',
        text: ['', '|3⟩'],
        textposition: 'middle right',
        textfont: { size: 18 },
        line: { color: '#34495e', width: 4 },
        hoverinfo: 'none',
        showlegend: false
    };

    // 2. Define the Layout (Arrows and removing axes)
    const layoutIV = {
        // Remove margins to maximize the drawing space
        margin: { t: 10, b: 10, l: 10, r: 10 }, 
        xaxis: {
            visible: false, // Hide the X axis entirely
            range: [0, 10]  // Set a fixed internal coordinate system
        },
        yaxis: {
            visible: false, // Hide the Y axis entirely
            range: [lowerY1-2*bwp, upperY+3*bwp] // Padding above and below the states
        },
        shapes: [
            // Pump Laser Bandwidth - Tail (State 1)
            {
                type: 'rect',
                x0: 1.5, x1: 2.5,      // Centers the box around the tail at X=2
                y0: lowerY1-bwp/2, y1: lowerY1+bwp/2,   // The "height" or bandwidth thickness
                fillcolor: '#e74c3c',
                opacity: 0.2,          // Makes it lightly shaded
                line: { width: 0 }     // Removes the hard border
            },
            // Pump Laser Bandwidth - Head (State 3)
            {
                type: 'rect',
                x0: 4.1, x1: 4.9,      // Centers around the head at X=4.5
                y0: lowerY1+freqp-bwp/2, y1: lowerY1+freqp+bwp/2,
                fillcolor: '#e74c3c',
                opacity: 0.2,
                line: { width: 0 }
            },
            // Stokes Laser Bandwidth - Tail (State 2)
            {
                type: 'rect',
                x0: 7.5, x1: 8.5,      // Centers around the tail at X=8
                y0: lowerY2-bwp/2, y1: lowerY2+bwp/2,
                fillcolor: '#3498db',
                opacity: 0.2,
                line: { width: 0 }
            },
            // Stokes Laser Bandwidth - Head (State 3)
            {
                type: 'rect',
                x0: 5.1, x1: 5.9,      // Centers around the head at X=5.5
                y0: lowerY2+freqs-bwp/2, y1: lowerY2+freqs+bwp/2,
                fillcolor: '#3498db',
                opacity: 0.2,
                line: { width: 0 }
            }
        ],
        annotations: [
            // Pump Arrow (|1⟩ to |3⟩)
            {
                ax: 2,         // Tail X
                ay: lowerY1,      // Tail Y (slightly above lower state)
                axref: 'x', 
                ayref: 'y',
                x: 4.5,        // Head X
                y: lowerY1+freqp,       // Head Y (slightly below upper state)
                xref: 'x', 
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1.5,
                arrowwidth: 2,
                arrowcolor: '#e74c3c' // Red
            },
            // Stokes Arrow (|2⟩ to |3⟩)
            {
                ax: 8,         // Tail X
                ay: lowerY2,      // Tail Y
                axref: 'x', 
                ayref: 'y',
                x: 5.5,        // Head X
                y: lowerY2+freqs,       // Head Y
                xref: 'x', 
                yref: 'y',
                showarrow: true,
                arrowhead: 2,
                arrowsize: 1.5,
                arrowwidth: 2,
                arrowcolor: '#3498db' // Blue
            }
        ],
        plot_bgcolor: 'rgba(0,0,0,0)', // Transparent background
        paper_bgcolor: 'rgba(0,0,0,0)'
    };

    // 3. Render the Plot
    Plotly.newPlot('energyPlot', [trace1, trace2, trace3], layoutIV, {staticPlot: true});
}

// Math/Physics Helpers
function Gaussion_Values(tp, ts, C0, CF,det,gboth) {
    return {
        Os: Math.sqrt((5.803548e11)*(4.33**2)/auI),
        ts: ts/fsperau,
        gs: gboth/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        ws: (3.960496+det)/evperAU,
        Op: Math.sqrt((5.803548e11)*(4.33**2)/auI),
        tp: tp/fsperau,
        gp: gboth/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        wp: (5.266919+det)/evperAU,
        E1: 0,
        E2: 5.266919/evperAU,
        E3: 1.306423/evperAU,
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
        let OmegaS = -(test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta)) * Math.cos(test.ws * t);
        let OmegaP = -(test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta)) * Math.cos(test.wp * t);
        return F(x, OmegaS, OmegaP);
    }
    let x0 = [Math.cos(test.alpha), 0, 0, 0, Math.sin(test.alpha), 0];
    let sol = numeric.dopri(test.t0, test.tf, x0, f, 1e-8, 20000);
    return [sol.x, numeric.transpose(sol.y)];
}

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
const inputsToWatch = ["C0", "CF", "t_s", "total_time", "detuning0", "duration", "time0", "timef", "toggle_curves"];

// Apply a 300ms debounce to prevent the UI from freezing when rapidly sliding
const debouncedUpdate = debounce(updatePlots, 300);

inputsToWatch.forEach(id => {
    document.getElementById(id).addEventListener('input', debouncedUpdate);
    document.getElementById(id).addEventListener('change', debouncedUpdate); 
});

// Initial Plot Generation on Page Load
window.onload = updatePlots;
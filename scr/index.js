// Convertion of units from SI to atomic units
const fsperau = (2.488843e-17)/(1.0e-15);
const auI = 3.50944e16;
const evperAU = 27.2114079527e0;

// Gets the Users Input Values & Treats them as floats
document.getElementById("calculate").onclick = function d(){
    const time0 = parseFloat(document.getElementById("time0").value);
    const timef = parseFloat(document.getElementById("timef").value);
    const tp_inp = parseFloat(document.getElementById("t_p").value);
    const ts_inp = parseFloat(document.getElementById("t_s").value);
    const alph = parseFloat(document.getElementById("alpha").value);
    const bet = parseFloat(document.getElementById("beta").value);

    // Puts the User input into the Gaussion_Values function to create the constants for the Gaussians
    var test = Gaussion_Values(tp_inp, ts_inp, alph, bet);
    var Delta = delta_creation(test.E1, test.E2, test.E3, test.wp, test.ws);
    test.delta = Delta;

    // Converts the time from femtoseconds to atomic units
    test.t0 = time0/fsperau;
    test.tf = timef/fsperau;

    const tp = test.tp;
    const ts = test.ts;
    const alpha = test.alpha;
    const beta = test.beta; 

    // Graphing the Envelope and Pulse Functions
    
    // Creates the arrays for the values of the x & y(X: time, Y: Power)
    const t_values = [];
    const envope_values_p = [];
    const envope_values_s = [];
    const steps = 2000
    const angle_va = []
    const beta_va = []
    const alpha_va = []

    const dt = (test.tf - test.t0)/steps;

    // Iterates through the # of steps to get a time value, then uses that time value to get the corresponding envelope value for both the Pump and Stokes light
    for(let i = 0; i<= steps; i++){
        const t = test.t0 + i*dt
        t_values.push(t*fsperau)

        var envelop_p = (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta)) 
        envope_values_p.push(envelop_p * Math.cos(test.wp * t))

        var envelop_s = (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta)) 
        envope_values_s.push(envelop_s * Math.cos(test.ws * t))

        angle_va.push(Math.atan(-envelop_p/envelop_s)/Math.PI)

        beta_va.push(bet)

        alpha_va.push(alph)
    }

    // Creates the two axis values for the pupms & graphs them
    const Gau_p = {
        x: t_values,
        y: envope_values_p,
        name: 'Pump Pulse',
        mode: 'lines',
        line: {color: 'blue', width: 1.5}
    }

    const Gau_s = {
        x: t_values,
        y: envope_values_s,
        name: 'Stokes Pulse',
        mode: 'lines',
        line: {color: 'red', width: 1.5}
    }

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
        line: {color: 'red', width: 1.5}
    }

    const Alpha_Line = {
        x: t_values,
        y: alpha_va,
        name: 'Alpha',
        mode: 'lines',
        line: {color: 'blue', width: 1.5}
    }
    const layout = {
        title: 'Envelope and Pulse Functions',
        xaxis: {title: 'Time (fs)'},
        yaxis: {title: 'Amplitude'},
        showlegend: true 
        };  
        
    
    const layoutIII = {
        title: 'Angle',
        xaxis: {title: 'Time (fs)'},
        yaxis: {title: 'Angle (rad)', range: [-1/2, 1/2]},
        showlegend: true
        };

    Plotly.newPlot('plot', [Gau_p, Gau_s], layout); 
    Plotly.newPlot('plotIII', [angles, Beta_Line, Alpha_Line], layoutIII);

    res = population_calculations_NRW(test);

    var ii = res[0].length;
    const tt = [];
    const c1=[];
    const c3 = [];
    for(let i = 0; i<ii; i++){
        tt.push(res[0][i]*fsperau)
        c1.push(res[1][0][i]**2 + res[1][1][i]**2)
        c3.push(res[1][4][i]**2 + res[1][5][i]**2)
    }
    const C1 ={
        x: tt,
        y: c1,
        name: 'C1',
        mode: 'lines',
        line: {color: 'blue', width: 1.5}
    }

    const C3 ={
        x: tt,
        y: c3,
        name: 'C3',
        mode: 'lines',
        line: {color: 'red', width: 1.5}
    }

    const layoutII = {
        title: 'Populations',
        xaxis: {title: 'Time (fs)'},
        yaxis: {title: 'Population',range: [0, 1.2]},
        showlegend: true
        };   

    Plotly.newPlot('plotII', [C1, C3], layoutII);

    const angle = []
    const time_values = []

   
} 

// Creates the constants for the diffrent Gaussians depending on the light being used
// 0s/p - The strength of the light/Amplitude  ts/p - Distance from center(orgin)   gs/p - Duration of the entire Gaussian(FW@HM)   ws/p - Frequency of the laser 
function Gaussion_Values(tp, ts, alpha, beta){
    const args = {
        Os: Math.sqrt((5.803548e11)*(4.33**2)/auI),
        ts: ts/fsperau,
        gs: 12/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        ws: 3.960496/evperAU,

        Op: Math.sqrt((5.803548e11)*(4.33**2)/auI),
        tp: tp/fsperau,
        gp: 12/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        wp: 5.266919/evperAU,

        E1: 0,
        E2: 5.266919/evperAU,
        E3: 1.306423/evperAU,

        mu12: 1,
        mu23: -1,

        alpha: alpha * Math.PI,
        beta: beta * Math.PI
    }
    return args
}

// Creates the Envelope Function for the Pulse
function envelope(t, t0, alpha, as, ap, mu){
    var exs = Math.exp(alpha*(t-t0))
    return  (1/mu)*alpha*(as-ap)*exs/
            ((1+exs)*Math.sqrt((1-as+(1-ap)*exs)*(as+ap*exs)));
}

// Creates the two Deltas(Laser Detuning Value), and returns the one delta
function delta_creation(E1, E2, E3, wp, ws){
    var delta_12 = E1 - E2 + wp
    var delta_23 = E3 - E2 + ws
    var delta = delta_12
    if(Math.abs(delta_12) != Math.abs(delta_23)){
        RaiseError("The two deltas are not equal, please check the values of the energies and frequencies")
    }
    return delta
}

// Creates the Gaussion that correlates with Stokes Light
function Gaussion_Creation_S(t, ts, gs){
    return Math.exp(-((t-ts)**2/gs**2))
}

// Creates the Gaussion that correlates with Pump Light
function Gaussion_Creation_P(t, tp, gp){
    return Math.exp(-((t-tp)**2/gp**2))
}


// Calculates the Population transfers acroos the levels
function population_calculations_NRW(test){

    // Creates the System of Equations for the Matrix Multiplication by the WaveForm
    function F(x, OmegaS, OmegaP){
        return[
        test.E1 * x[1] - (OmegaP * x[3]),

        -(test.E1 * x[0]) + (OmegaP * x[2]),

        test.E2*x[3]-(OmegaP * x[1]) - (OmegaS * x[5]),

        -test.E2*x[2]+(OmegaP * x[0]) + (OmegaS * x[4]),
        
        test.E3 * x[5] - (OmegaS * x[3]),

        -(test.E3 * x[4]) + (OmegaS * x[2]) 
        
    ]
    }

    // Creates the Omegas, as well as inputing these into the system of equations to get the final values for the populations
    function f(t, x){
        var OmegaS = (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.
        beta)) * Math.cos(test.ws * (t-test.ts))
        
        var OmegaP = (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta)) * Math.cos(test.wp * (t-test.tp))

        return F(x, OmegaS, OmegaP)
            
    }

    
    var x0 = [Math.cos(test.alpha), 0, 0, 0, Math.sin(test.alpha), 0]

    // Uses the dopri function to get the final values for the populations, and returns these values
    var sol = numeric.dopri(test.t0, test.tf, x0, f, 1e-8, 100000)
    y_line = sol.y
    time = sol.x
    y_ret = numeric.transpose(y_line)
    return [time, numeric.transpose(y_line)]
}



// Calculates the Population transfers acroos the levels
function population_calculations(test){

    // Creates the System of Equations for the Matrix Multiplication by the WaveForm
    function F(x, OmegaS, OmegaP){
        return[
        test.delta * x[1] - (OmegaP * x[3])/2,
        -(test.delta * x[0]) + (OmegaP * x[2])/2,
        -(OmegaP * x[1])/2 - (OmegaS * x[5])/2,
        (OmegaP * x[0])/2 + (OmegaS * x[4])/2,
        test.delta * x[5] - (OmegaS * x[3])/2,
        -(test.delta * x[4]) - (OmegaS * x[2])/2 
        ]
    }

    // Creates the Omegas, as well as inputing these into the system of equations to get the final values for the populations
    function f(t, x){
        var OmegaS = (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta))
        var OmegaP = (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta))
        return F(x, OmegaS, OmegaP)
            
    }

    
    var x0 = [Math.cos(test.alpha), 0, 0, 0, Math.sin(test.alpha), 0]

    // Uses the dopri function to get the final values for the populations, and returns these values
    var sol = numeric.dopri(test.t0, test.tf, x0, f, 1e-8, 2000)
    y_line = sol.y
    time = sol.x
    y_ret = numeric.transpose(y_line)
    return [time, numeric.transpose(y_line)]
}


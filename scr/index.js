//Convertion of units from SI to atomic units
const fsperau = (2.488843e-17)/(1.0e-15);
const auI = 3.50944e16;
const evperAU = 27.2114079527e0;

// Gets the Users Input Values & Treats them as floats
document.getElementById("calculate").onclick = function d(){
    const time0 = parseFloat(document.getElementById("time0").value);
    const timef = parseFloat(document.getElementById("timef").value);
    const tp_inp = parseFloat(document.getElementById("tp").value);
    const alph = parseFloat(document.getElementById("alpha").value);
    const bet = parseFloat(document.getElementById("beta").value);

    var test = Gaussion_Values(tp_inp, alph, bet);
    var Delta = delta_creation(test.E1, test.E2, test.E3, test.wp, test.ws);

    const tp = test.tp;
    const alpha = test.alpha;
    const beta = test.beta;
    var Delta = delta_creation(test.E1, test.E2, test.E3, test.wp, test.ws);
    const time_initial = time0;
    const time_final = timef;
      
    // First Shot at Graphing the Pulses and Envolope Functions
    
    const t_values = [];
    const envope_values = [];
    const pulse_values = [];
    const steps = 2000
    const dt = (time_final - time_initial)/steps;
    const w0 = test.ws;
    const test_as = .8
    const test_ap = .2
    const test_mu = 1

// The problem lies in the math for as & ap, so look back at the Gaussian_Values math!

    for (let i = 0; i<= steps; i++){
        const t = time_initial + i*dt
        t_values.push(t)

        const t0_center = (time_initial + time_final)/2

        const env = envelope(t, t0_center, alpha, test_as, test_ap, test_mu)
        envope_values.push(env)

        const pul = env *Math.sin(5*t)
        pulse_values.push(pul)
    }

    const traceEnvolpe = {
        x: t_values,
        y: envope_values,
        mode: 'lines',
        name: 'Envelope Function',
        line: {color: 'blue', width: 2}
    };

    const tracePulse = {
        x: t_values,
        y: pulse_values,
        mode: 'lines',
        name: 'Pulse Function',
        line: {color: 'red', width: 1.5}
    };

    const layot = {
        title: 'Envelope and Pulse Functions',
        xaxis: {title: 'Time (fs)'},
        yaxis: {title: 'Amplitude'},
        template: 'plotly_dark',
        showlegend: true
    };
    Plotly.newPlot('plot', [traceEnvolpe, tracePulse], layot);

}


// Ensures the user input values are valid

// Creates the constants for the diffrent Gaussians depending on the light being used
// 0s/p - The strength of the light/Amplitude  ts/p - Distance from center(orgin)   gs/p - Duration of the entire Gaussian(FW@HM)   ws/p - Frequency of the laser 
function Gaussion_Values(tp, alpha, beta){
    const args = {
        Os: Math.sqrt((5.803548*(10**11)*(4.33**2))/auI),
        ts: 50/fsperau,
        gs: 12/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        ws: 9.960496/evperAU,

        Op: Math.sqrt((5.803548*(10**11)*(4.33**2))/auI),
        tp: tp/fsperau,
        gp: 12/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        wp: 11.266919/evperAU,

        E1: 0,
        E2: 11.266919/evperAU,
        E3: 1.306423/evperAU,

        mu12: 1,
        mu23: -1 * 1,

        alpha: alpha * Math.PI,
        beta: beta * Math.PI,
    }
    return args
}

// Creates the Envelope Function for the Pulse
function envelope(t, t0, alpha, as, ap, mu){
    var exs = Math.exp(alpha*(t-t0))
    return  (1/mu)*alpha*(as-ap)*exs/
            ((1+exs)*Math.sqrt((1-as+(1-ap)*exs)*(as+ap*exs)));
}

// Create the Pulse Function
function pulse(t, w0){
    return envelope(t)*Math.sin(w0*t)
}

// Creates the two Deltas(Laser Detuning Value), and returns the one delta
function delta_creation(E1, E2, E3, wp, ws){
    var delta_12 = E1 - E2 + wp
    var delta_23 = E3 - E2 + ws
    var delta = delta_12
    return delta
}

// Creates the Gaussion that correlates with Stokes Light
function Gaussion_Creation_S(t, ts, gs){
    return Math.exp((-(t-ts))**2/gs**2)
}

// Creates the Gaussion that correlates with Pump Light
function Gaussion_Creation_P(t, tp, gp){
    return Math.exp((-(t-tp))**2/gp**2)
}

// Calculates the Population transfers acroos the levels
function population_calculations(t, t0, tf, alpha, beta, E, ts, tp, gs, gp, delta){

    console.log(`t0: ${t0}, tf: ${tf}, alpha: ${alpha}, beta: ${beta}, E: ${E}, ts: ${ts}, tp: ${tp}, gs: ${gs}, gp: ${gp}, delta: ${delta}`)

    // Calculates the Omega_S coupler term
    function OmegaS_cal(t){
        return E * Math.sin(alpha) * Gaussion_Creation_S(t, ts, gs) +
               E * Math.sin(beta) * Gaussion_Creation_P(t, tp, gp)
    }
    console.log(`OmegaS: ${OmegaS_cal(t)}`)
    
    // Calculates the Omega_P coupler term
    function OmegaP_cal(t){
        return E * Math.cos(alpha) * Gaussion_Creation_S(t, ts, gs) +
               E * Math.cos(beta) * Gaussion_Creation_P(t, tp, gp)
    }
    console.log(`OmegaP: ${OmegaP_cal(t)}`)

    // Creates the System of Equations for the Matrix Multiplication by the WaveForm
    function F(t0, tf, x, OmegaS, OmegaP, delta){
        return[
        delta * x[1] - (OmegaP * x[3])/2,
        -(delta * x[0]) + (OmegaP * x[2])/2,
        -(OmegaP * x[1])/2 - (OmegaS * x[5])/2,
        (OmegaP * x[0])/2 + (OmegaS * x[4])/2,
        delta * x[5] - (OmegaS * x[3])/2,
        -(delta * x[4]) - (OmegaS * x[2])/2 
        ]
    }
    //console.log(`F: ${F(t0, tf, [Math.cos(alpha), 0, 0, 0, Math.sin(alpha), 0], OmegaS_cal(t), OmegaP_cal(t), delta)}`)

    function f(t, x){
        
        return F(t0, tf, x, OmegaS_cal(t), OmegaP_cal(t), delta)
            
    }
    //console.log(`f: ${f(t, [Math.cos(alpha), 0, 0, 0, Math.sin(alpha), 0])}`)   
    var x0 = [Math.cos(alpha), 0, 0, 0, Math.sin(alpha), 0]
    var sol = numeric.dopri(t0, tf, x0, f, 1e-8, 2000)
    //console.log(x0)
    //console.log(sol)
    y_line = sol.y
    time = sol.x
    //console.log(`time: ${time}`)
    //console.log(`y_line: ${y_line}`)

    //Make it print arrays and Time dependant
    //Next Step is to extract the time and the 6 dimensional table
}


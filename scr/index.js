//Convertion of units from SI to atomic units
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

    var test = Gaussion_Values(tp_inp, ts_inp, alph, bet);
    var Delta = delta_creation(test.E1, test.E2, test.E3, test.wp, test.ws);
    test.delta = Delta;
    test.t0 = time0/fsperau;
    test.tf = timef/fsperau;

    const tp = test.tp;
    const ts = test.ts;
    const alpha = test.alpha;
    const beta = test.beta; 

    //Graphing the Envelope and Pulse Functions
    
    const t_values = [];
    const envope_values_p = [];
    const envope_values_s = [];
    const steps = 2000
    time_final = timef / fsperau
    time_initial = time0 / fsperau

    const dt = (time_final - time_initial)/steps;

    for(let i = 0; i<= steps; i++){
        const t = time_initial + i*dt
        t_values.push(t)

        const envelop_p = (test.Op * test.mu12) * (Gaussion_Creation_P(t, tp, test.gp) * Math.sin(alpha) + Gaussion_Creation_S(t, ts, test.gs) * Math.sin(beta))

        envope_values_p.push(envelop_p)

        const envelop_s = (test.Os * test.mu23) * (Gaussion_Creation_P(t, tp, test.gp) * Math.cos(alpha) + Gaussion_Creation_S(t, ts, test.gs) * Math.cos(beta))
        envope_values_s.push(envelop_s)
    }

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

    const layout = {
        title: 'Envelope and Pulse Functions',
        xaxis: {title: 'Time (fs)'},
        yaxis: {title: 'Amplitude'},
        showlegend: true
        };   
    Plotly.newPlot('plot', [Gau_p, Gau_s], layout); 

    population_calculations(test);
}
// Creates the constants for the diffrent Gaussians depending on the light being used
// 0s/p - The strength of the light/Amplitude  ts/p - Distance from center(orgin)   gs/p - Duration of the entire Gaussian(FW@HM)   ws/p - Frequency of the laser 
function Gaussion_Values(tp, ts, alpha, beta){
    const args = {
        Os: Math.sqrt((5.803548*(10**11)*(4.33**2)))/auI,
        ts: ts/fsperau,
        gs: 12/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        ws: 9.960496/evperAU,

        Op: Math.sqrt((5.803548*(10**11)*(4.33**2)))/auI,
        tp: tp/fsperau,
        gp: 12/fsperau * 1/Math.sqrt(2*Math.log(2.0)),
        wp: 11.266919/evperAU,

        E1: 0,
        E2: 11.266919/evperAU,
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
function population_calculations(test){

    console.log(`t0: ${test.t0}, tf: ${test.tf}, alpha: ${test.alpha}, beta: ${test.beta}, E: ${test.E}, ts: ${test.ts}, tp: ${test.tp}, gs: ${test.gs}, gp: ${test.gp}, delta: ${test.delta}`)

    // Calculates the Omega_S coupler term
    function OmegaS_cal(t){
        return 
        (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta))
    }
    //console.log(`OmegaS: ${OmegaS_cal(t)}`)
    
    // Calculates the Omega_P coupler term
    function OmegaP_cal(t){
        return 
        (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta))
    }
    //    console.log(`OmegaP: ${OmegaP_cal(t)}`)

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

    function f(t, x){
        var OmegaS = (test.Os * test.mu23) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.cos(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.cos(test.beta))
        var OmegaP = (test.Op * test.mu12) * (Gaussion_Creation_P(t, test.tp, test.gp) * Math.sin(test.alpha) + Gaussion_Creation_S(t, test.ts, test.gs) * Math.sin(test.beta))
        //console.log(`t: ${t}, x: ${x}`)
        //console.log(`OmegaS: ${OmegaS}, OmegaP: ${OmegaP}`)
        return F(x, OmegaS, OmegaP)
            
    }

    
    var x0 = [Math.cos(test.alpha), 0, 0, 0, Math.sin(test.alpha), 0]

    //console.log(f(test.tp, x0)) 
    var sol = numeric.dopri(test.t0, test.tf, x0, f, 1e-8, 20000)
    y_line = sol.y
    time = sol.x
    console.log('Time steps: ', time)
    console.log(y_line.length)
    console.log(y_line[0].length)
    console.log('Done with one call to dopri')
    return [time, y_line]
}

//console.log(population_calculations(test.t0, test.tf, test.alpha, test.beta, test.E, test.ts, test.tp, test.gs, test.gp, test.delta))


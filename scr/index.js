//Convertion of units from SI to atomic units
const fsperau = (2.488843e-17)/(1.0e-15);
const auI = 3.50944e16;
const evperAU = 27.2114079527e0;

// Gets the Users Input Values
const time0 = document.getElementById("time0").value;
const timef = document.getElementById("timef").value;
const tp = document.getElementById("tp").value;
const alpha = document.getElementById("alpha").value;
const beta = document.getElementById("beta").value;


// Creates the constants for the diffrent Gaussians depending on the light being used
// 0s/p - The strength of the light/Amplitude  ts/p - Distance from center(orgin)   gs/p - Duration of the entire Gaussian(FW@HM)   ws/p - Frequency of the laser 
function Gaussian_Values(tp, alpha, beta){
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
var test = Gaussian_Values(tp, alpha, beta);
var test2 = Gaussion_Creation_S(time0, test.ts, test.gs);
console.log(`Test2: ${test2}`);
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
function population_calculations(t, t0, tf, alpha, beta, E, ts, tp, gs, gp, delta, pulse){

    // Calculates the Omega_S coupler term
    function OmegaS_cal(t){
        return E * Math.sin(alpha) * Gaussian_Creation_S(t, ts, gs) +
               E * Math.sin(beta) * Gaussion_Creation_P(t, tp, gp)
    }
    
    // Calculates the Omega_P coupler term
    function OmegaP_cal(t){
        return E * Math.cos(alpha) * Gaussian_Creation_S(t, ts, gs) +
               E * Math.cos(beta) * Gaussion_Creation_P(t, tp, gp)
    }

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

    function f(t, x){
        
        return F(t0, tf, x, OmegaS_cal(t), OmegaP_cal(t), delta)
            
    }   
    var x0 = [cos(alpha), 0, 0, 0, sin(alpha), 0]
    sol = numeric.dopri(t0, tf, x0, f, 1e-8, 2000)
    y_line = sol.y
    time = sol.x
}


/*
Pulse Function
      function pulse(t, envelope, w0, ...)
          return envelope(t)*Math.sin(w0*t)

 Envolope Function for Pulse
      function envolope(t, t0, alpha, as, ap mu)
          var exs = Math.exp(alpha*(t-t0))
          return  (1/mu)*alpha*(as-ap)*exs/
                  ((1+ex*Math.sqrt((1-as+(1-ap)*ex)(as+ap*ex)));

Gaussian Creation Function
function gaussian(t0, Os, Op, gs, gp, ts, tp, alpha, beta, ss, sp){
    return {
        Pump_Light: (Math.cos(alpha) * (Math.exp(-((t0 - ts)**2)/(gs)**2))),
        Stokes_Light: (Math.sin(beta) * (Math.exp(-((t0 - tp)**2)/(gp)**2)))

    };
}  


 Creates the Varibles to be used in the equations
function making_omegas(t, args){
     Makes Omega_Pump & Omega_Stoke
    var Omega_P = -args.mu12 * args.Op * (Math.exp(-((t-args.ts)/(args.gs)**2)) * Math.sin(args.alpha) + 
            Math.exp(-((t-args.tp)/(args.gp)**2)) * Math.sin(args.beta)) * Math.cos(args.wp * t);


    var Omega_S = -args.mu12 * args.Op * (Math.exp(-((t-args.ts)/(args.gs)**2)) * Math.cos(args.alpha) + 
            Math.exp(-((t-args.tp)/(args.gp)**2)) * Math.cos(args.beta)) * Math.cos(args.ws * t);

     Makes the Deltas that coupls energy level 1&2 and 2&3
     
}
*/



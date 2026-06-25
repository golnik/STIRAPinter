//Convertion of units from SI to atomic units
const fsperau = (2.488843e-17)/(1.0e-15);
const auI = 3.50944e16;
const evperAU = 27.2114079527e0;

// Creates the Varibles to be used in the equations
function making_omegas(t, args){
    // Makes Omega_Pump & Omega_Stoke
    var Omega_P = -args.mu12 * args.Op * (Math.exp(-((t-args.ts)/(args.gs)**2)) * Math.sin(args.alpha) + 
            Math.exp(-((t-args.tp)/(args.gp)**2)) * Math.sin(args.beta)) * Math.cos(args.wp * t);


    var Omega_S = -args.mu12 * args.Op * (Math.exp(-((t-args.ts)/(args.gs)**2)) * Math.cos(args.alpha) + 
            Math.exp(-((t-args.tp)/(args.gp)**2)) * Math.cos(args.beta)) * Math.cos(args.ws * t);

    // Makes the Deltas that coupls energy level 1&2 and 2&3
    var delta_12 = args.E1 - args.E2 + args.wp;
    var delta_23 = args.E3 - args.E2 + args.ws;
    //var delta = delta_12 = delta_23 
}
function population_calculations(alpha, beta, E, ts, tp, t0, tf, gs, gp, delta){

    function OmegaS_cal(t){
        return (E * Math.sin(alpha) * (Math.exp(-((t0 - ts)**2)/(gs)**2)))
               + (E * Math.sin(beta) * (Math.exp(-((t0 - tp)**2)/(gp)**2)))
    }

    function OmegaP_cal(t){
        return (E * Math.cos(alpha) * (Math.exp(-((t0 - ts)**2)/(gs)**2)))
               + (E * Math.cos(beta) * (Math.exp(-((t0 - tp)**2)/(gp)**2)))
    }

    function f(t, x){
        var x0 = [cos(alpha), 0, 0, 0, sin(alpha), 0]
        return function F(t0, tf, x, OmegaS_cal, OmegaP_cal, delta){
            delta * x[1] - (OmegaP_cal * x[3])/2,
            -(delta * x[0]) + (OmegaP_cal * x[2])/2,
            -(OmegaP_cal * x[1])/2 - (OmegaS_cal * x[5])/2,
            (OmegaP_cal * x[0])/2 + (OmegaS_cal * x[4])/2,
            delta * x[5] - (OmegaS_cal * x[3])/2,
            -(delta * x[4]) - (OmegaS_cal * x[2])/2
    }   
    sol = numeric.dopri(t0, tf, x0, f, 1e-8, 2000)
    y_line = sol.y
    time = sol.x
    }
}




// Pulse Function
//      function pulse(t, envelope, w0, ...)
//          return envelope(t)*Math.sin(w0*t)

// Envolope Function for Pulse
//      function envolope(t, t0, alpha, as, ap mu)
//          var exs = Math.exp(alpha*(t-t0))
//          return  (1/mu)*alpha*(as-ap)*ex/
//                  ((1+ex*Math.sqrt((1-as+(1-ap)*ex)(as+ap*ex)));

// Gaussian Creation Function
//      function gaussian(t0, Os, Op, gs, gp, ts, tp, alpha, beta, ss, sp)
//          return{Pump_Light: Math.sin(alpha) * Math.exp(-(t0-ts)**2/(ss)**2 + Math.sin(beta) * Math.exp(-(t0-tp)**2/(sp)**2)), 
//
//                  Stokes_Light: Math.cos(alpha) * Math.exp(-(t0-ts)**2/(ss)**2 + Math.cos(beta) * Math.exp(-(t0-tp)**2/(sp)**2)))
             







// Creates the constants for the diffrent Gaussians depending on the light being used
// 0s/p - The strength of the light/Amplitude  ts/p - Distance from center(orgin)   gs/p - Duration of the entire Gaussian(FW@HM)   ws/p - Frequency of the laser 
function Gaussian_Values(tp, alpha, beta, ll){
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

        alpha: 1/3 * Math.PI,
        beta: 1/4 * Math.PI,
    }
    return args
}


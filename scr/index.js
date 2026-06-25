//Convertion of units from SI to atomic units
const fsperau = (2.488843*(10**-17))/(1.0*(10**-15));
const auI = 3.50944*(10**16);
const evperAU = 27.2114079527*(10**0);

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

// Creates the System of Equations being used instead of doing Matrix Multiplication
//function System_Equations(d12, d23, Os, Op, t, ts, tp){
    //db_1 = making_omegas.delta_12 * a1 - (making_omegas.Omega_P * a2)/2
    //da_1 = making_omegas.delta_12 * b1 - (making_omegas.Omega_P * b2)/2

    //db_2 = -(making_omegas.Omega_P * a1)/2 - (making_omegas.Omega_S * a3)/2
    //da_2 = -(making_omegas.Omega_P * b1)/2 - (making_omegas.Omega_S * b3)/2

    //db_3 = making_omegas.delta_23 * a3 - (making_omegas.Omega_S * a2)/2
    //da_3 = making_omegas.delta_23 * b3 - (making_omegas.Omega_S * b2)/2

    // return   d12 * ts - (Os * tp)/2
    //          d12 * 
//}

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

// System of Equation  Functions F(x)
//      function System(Os, Op, t, delta)6yyvggbxdxdcvgvgvcvggggggggggggggggggggggggggggggggggggbgvvvvvvvvvvvvvvvvvvv   
\'
\\\

\
\\
\

\                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       N       '

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


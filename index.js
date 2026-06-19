//Calculates the Coupler Terms Omega_P and Omega_S
function field(t,args){
    return args['Omega_p']*Math.exp(-((t-args['tp'])/args['gp']))**2*(Math.cos(args['wp']*(t-args['tp']))) + 
           args['Omega_s']*Math.exp(-((t-args['ts'])/args['gs']))**2*(Math.cos(args['ws']*t-args['wp']*args['tp']))
}
function hamiltonian(t,args){
    var omega_p = -args['mu12']*field(t,args);
    var omega_s = -args['mu23']*field(t.args);
}




//Pratice Graphing 
const xArray = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0, 3.5, 4.0];
const yArray = [22, 21, 10, 19, 18, 17, 16, 15, 14, 13, 12];

const data = [{
    x:xArray,
    y:yArray,
    mode:"line+markers",
    line:{color: "blue"},
    marker:{size:8}
}];

const layout = {
    title: "Car Engine Size vs. Fuel Efficiency",
    xaxis: {range:[1,4.5], title: "Enginge Size(L)"},
    yaxis: {range:[10,25], title: "Fuel Efficiency(MPG)"}
};
Plotly.newPlot("myPraticePlot", data, layout)

const a = np.array([[1, 2], [3, 4]]);
const d = a.multiply(2);
console.log(d);
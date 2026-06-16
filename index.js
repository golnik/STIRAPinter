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
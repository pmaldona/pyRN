function generateGraph(containerName, ntwrk) {
    console.log("generate graph");
    console.log(containerName);
    let data = {
        nodes: ntwrk.nodes ? ntwrk.nodes : [],
        edges: ntwrk.edges ? ntwrk.edges : [],
    };
    return drawVis(containerName, data, ntwrk.options);
}

function drawVis(containerName, data, options) {
    console.log(containerName)
    let container = document.getElementById(containerName); 
    console.log(container)
    let _network = new vis.Network(container, data, options);
    _network.setOptions(options);
    return _network;
}

const GraphType = {
    Network: 'network',
    Hasse: 'hasse',
    OrgNetwork: 'org_network'
}

function setVisObject(type, nodes, edges, options, physics, layout) {
    console.log(nodes);
    console.log(vis);
    nodes = new vis.DataSet(nodes);
    edges = new vis.DataSet(edges);
    options = JSON.parse(options);
    options.edges.font = {color:'#000000', strokeColor: '#000000', strokeWidth: 0.5};
    if(physics == true) {
        options.physics = {enabled: false};
    }
    if(layout) {
        options["layout"]= layout;
    }
    let network = generateGraph(type, {nodes: nodes, edges: edges, options: options});
    return network;
}

export { GraphType, setVisObject, drawVis, generateGraph };
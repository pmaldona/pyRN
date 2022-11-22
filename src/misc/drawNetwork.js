function generateGraph(containerName, ntwrk) {
    console.log("generate graph");
    console.log(containerName);
    let container = document.getElementById(containerName);
    let data = {
        nodes: ntwrk.nodes ? ntwrk.nodes : [],
        edges: ntwrk.edges ? ntwrk.edges : [],
    };

    let _network = new vis.Network(container, data, ntwrk.options);
    _network.setOptions(ntwrk.options);
    return _network;
}

const GraphType = {
    Network: 'network',
    Hasse: 'hasse',
    OrgNetwork: 'org_network'
}

function setVisObject(type, nodes, edges, options, physics, layout) {
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

export { GraphType, setVisObject };
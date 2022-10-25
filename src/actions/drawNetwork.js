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

export { generateGraph };
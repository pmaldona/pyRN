<script>
    //import { Network } from 'vis-network';
    //import { DataSet } from 'vis-data/peer';
    export let genNetwork;
    export let initialValues;

    let nodes = initialValues.nodes;
    let edges = initialValues.edges;
    let options = initialValues.hasse ? initialValues.hasse.options : {};

    function generateGraph() {
        let container = document.getElementById("network");
        let data = {
            nodes: nodes ? nodes : [],
            edges: edges ? edges : [],
        };
    
        let network = new vis.Network(container, data, options);
        console.log(options);
        network.setOptions(options);
        console.log(network)
    }

    if(initialValues.name && nodes == undefined) {
        genNetwork().then(result => {
            nodes = new vis.DataSet(result.nodes);
            let e = [];
            for(let i = 0; i < result.edges.length; i++) {
                let edge = result.edges[i];
                console.log(edge);
                if(edge.title){
                    edge.label = edge.title;
                }
                
            }
            console.log(result);
            edges = new vis.DataSet(result.edges);
            options = JSON.parse(result.options);
            options.edges.font = {color:'#000000', strokeColor: '#000000', strokeWidth: 0.5};
            generateGraph();
        });
        
    }
    
</script>

<main>
    <div id="network"></div>
    {#if nodes == undefined}
        <div id="loader">
            <div id="circle">
                <div class="preloader-wrapper big active">
                    <div class="spinner-layer spinner-blue-only">
                        <div class="circle-clipper left">
                            <div class="circle"></div>
                        </div>
                        <div class="gap-patch">
                            <div class="circle"></div>
                        </div>
                        <div class="circle-clipper right">
                            <div class="circle"></div>
                        </div>
                    </div>
                </div>
            </div>
            <h5>
                Creating Network
            </h5>
        </div>
    {:else if initialValues.name == undefined}
        <h5>
            Please open a network file.
        </h5>
    {/if}
	
</main>

<style>
    #network {
        width: 450px;
        height: 450px;
        border: 1px solid lightgray;
        margin: 2px;
    }

    #loader {
        position: absolute;
        display: block;
        top: 50%;
        left: 50%;
        transform: translateX(-50%) translateY(-50%);
    }
    #circle {
        margin-left: 33%;
    }
</style>
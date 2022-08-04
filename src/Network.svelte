<script>
    //import { Network } from 'vis-network';
    //import { DataSet } from 'vis-data/peer';
    import { network } from './store/NetworkStore';
    import { filename } from './store/FileStore';
    export let genNetwork;

    let nodes;
    let edges;
    let options = {};

    if($filename != "" && nodes == undefined) {
        genNetwork().then(result => {
            nodes = new vis.DataSet(result.nodes);
            for(let i = 0; i < result.edges.length; i++) {
                let edge = result.edges[i];
                //console.log(edge);
                if(edge.title){
                    edge.label = edge.title;
                }
                
            }
            
            edges = new vis.DataSet(result.edges);
            options = JSON.parse(result.options);
            options.edges.font = {color:'#000000', strokeColor: '#000000', strokeWidth: 0.5};
            network.set_network({nodes: nodes, edges: edges, options: options});
        });   
    }
    
</script>

<main>
    <div id="network"></div>
    {#if $filename == ""}
        <h5 style="position: absolute; top: 55px; left: 10px;">
            Please open a network file.
        </h5>
    {:else if nodes == undefined}
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
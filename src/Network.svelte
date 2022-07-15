<script>
    import { Network, DataSet } from 'vis-network';
    export let genNetwork;
    export let initialValues;

    let nodes = initialValues.nodes;
    let edges = initialValues.edges;

    function generateGraph() {
        let container = document.getElementById("network");
        let data = {
            nodes: nodes ? nodes : [],
            edges: edges ? edges : [],
        };
        console.log(nodes);
        console.log(data);

        let options = {};
    
        let network = new Network(container, data, options);
    }

    if(initialValues.name && nodes == undefined) {
        genNetwork();
    }
    
</script>

<main>
    {#if nodes != undefined}
        <div id="network"></div>
    {:else if initialValues.name}
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
    {:else}
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
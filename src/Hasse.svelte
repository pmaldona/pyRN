<script>
    import { Network } from 'vis-network';
    import { DataSet } from 'vis-data/peer';
    export let genHasse;
    export let initialValues;
    
    let nodes = initialValues.hasse ? initialValues.hasse.nodes : undefined;
    let edges = initialValues.hasse ? initialValues.hasse.edges : undefined;
    let options = initialValues.hasse ? initialValues.hasse.options : {};
    
    function generateHasse() {
        let container = document.getElementById("network");
        let data = {
            nodes: nodes ? nodes : new DataSet([]),
            edges: edges ? edges : new DataSet([]),
        };
    
        let hasse = new Network(container, data, options);
    }

    if(initialValues.name && nodes == undefined) {
        genHasse().then(result => {
            console.log(result);
            nodes = new DataSet(result.nodes);
            edges = new DataSet(result.edges);
            options = JSON.parse(result.options);
            console.log(options);
            generateHasse();
        });
    }

    if(nodes != undefined) {
        console.log(nodes);
        generateHasse();
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
                Creating Hasse
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
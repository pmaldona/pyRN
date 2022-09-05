<script>
    //import { Network } from 'vis-network';
    //import { DataSet } from 'vis-data/peer';
    import { network } from './store/NetworkStore';
    import { filename } from './store/FileStore';
    import { protoSyn } from './store/ProtoSynStrStore';
    export let genNetwork;
    export let genProtoSyn;

    let nodes;
    let edges;
    let options = {};

    let graphs = [
		{ id: 1, text: `Reaction Network` },
		{ id: 2, text: `Protosynergies` },
	];

	let selected = graphs[0];

    let statistics = {
        species_count: 0,
        reaction_count: 0,
        generators: 0,
        transitions: 0
    };

    if($filename != "" && nodes == undefined) {
          drawNetwork();
    }

    function genRNStr() {
        genNetwork().then(result => {
            nodes = new vis.DataSet(result.nodes);
            for(let i = 0; i < result.edges.length; i++) {
                let edge = result.edges[i];
                // console.log(edge);
                if(edge.title){
                    edge.label = edge.title;
                }
                edge.smooth = {enabled: false};
                for(let j = 0; j < result.edges.length; j++) {
                    if(j != i) {
                        if(edge.from == result.edges[j].to && edge.to == result.edges[j].from) {
                            edge.smooth = {enabled: true, type: "discrete", roundness: 0.7};
                        }
                    }
                }
                edge.physics = false;
            }
            
            edges = new vis.DataSet(result.edges);
            options = JSON.parse(result.options);
            options.edges.font = {color:'#000000', strokeColor: '#000000', strokeWidth: 0.5}
            options.physics = {enabled: false};
            statistics.species_count = result.species_count;
            statistics.reaction_count = result.reaction_count;
            network.set_network({nodes: nodes, edges: edges, options: options, species_count: result.species_count, reaction_count: result.reaction_count});
        }); 
    }

    function genProSyn() {
        genProtoSyn().then(result => {
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
            statistics.generators = result.generators;
            statistics.transitions = result.transitions;
            protoSyn.set_proto_syn({nodes: nodes, edges: edges, options: options});
        }); 
    }

    function drawNetwork() {
        if(selected.id == 2) {
            genProSyn();
        } else {
            genRNStr();
        }
    }
    
</script>

<main>
    <div style="display: flex;">
        <div id="network"></div>
        <div style="margin: 15px;">
            <select bind:value={selected} on:change="{() => console.log("change")}" style="display:block;">
                {#each graphs as graph}
                    <option value={graph}>
                        {graph.text}
                    </option>
                {/each}
            </select>
            <!-- svelte-ignore a11y-missing-attribute -->
            <a class="waves-effect waves-light btn" style="margin-top: 5px;" on:click={drawNetwork}>redraw lattice</a>
        </div>
    </div>
    <div style="position: absolute; bottom: 51px; right: 5px;">
        <h6>Network statistics:</h6>
        {#if selected.id == 1}
            <textarea readonly style="resize: none; width: 300px; height: 100px;">
#Species: {statistics.species_count}
#Reactions: {statistics.reaction_count}
            </textarea>
        {:else}
            <textarea readonly style="resize: none; width: 300px; height: 100px;">
#Partition-Generators: {statistics.generators}
#Generative-Transitions: {statistics.transitions}
            </textarea>
        {/if}
       
    </div>
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
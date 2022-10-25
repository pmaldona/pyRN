<script>
    import { hasse } from './store/HasseStore';
    import { layout, synStr } from './store/SynStrStore'; 
    import { filename } from './store/FileStore';
    import Slider from '@bulatdashiev/svelte-slider';
    export let genHasse;
    export let genSynergeticStructure;

    let nodes;
    let edges;
    let options = {};
    
    let graphs = [
		{ id: 1, text: `Reactive Organisation Structures` },
		{ id: 2, text: `Synergistic Structures` }
	];

	let selected = graphs[0];

    let statistics = {
        species_count: 0,
        reaction_count: 0,
        basic: 0,
        closed: 0,
        basic_closed: 0,
        ssm: 0,
        orgs: 0,
        union: 0,
        syn_union: 0,
        conn_orgs: 0
    };

    let x_space = [75, 200];
    let y_space = [75, 200];

    let useLayout = true;

    if($filename && nodes == undefined) {
        drawLattice();
    }

    function genOrgStr() {
        genHasse().then(result => {
            console.log(result);
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
            options["layout"]= {"improvedLayout": false};
            statistics.species_count = result.species_count;
            statistics.reaction_count = result.reaction_count;
            hasse.set_hasse({nodes: nodes, edges: edges, options: options});
        });
    }

    function genSynStr() {
        genSynergeticStructure().then(result => {
            console.log(result);
            let n = useLayout ? layout(result.nodes, x_space[0], y_space[0]) : result.nodes;
            nodes = new vis.DataSet(n);
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
            statistics.basic = result.basic;
            statistics.closed = result.closed;
            statistics.basic_closed = result.basic_closed;
            statistics.ssm = result.ssm;
            statistics.orgs = result.orgs;
            statistics.union = result.union;
            statistics.syn_union = result.syn_union;
            statistics.conn_orgs = result.conn_orgs;
            synStr.set_syn_str({nodes: nodes, edges: edges, options: options});
        });
    }

    function drawLattice() {
        console.log(selected.id);
        if(selected.id == 2) {
            genSynStr();
        } else {
            genOrgStr();
        }
    }
    // if(nodes != undefined) {
    //     console.log(nodes);
    //     generateHasse();
    // }
    
</script>

<main>
    <div style="display: flex;">
        <div id="hasse"></div>
        <div style="margin: 15px;">
            <select bind:value={selected} on:change="{() => console.log(selected)}" style="display:block">
                {#each graphs as graph}
                    <option value={graph}>
                        {graph.text}
                    </option>
                {/each}
            </select>
            X-Space: {x_space[0]}
            <Slider bind:value = {x_space} max="200"/>

            Y-Space: {y_space[0]}
            <Slider bind:value = {y_space} max="200"/>

            Use Layout:
            <label>
                <input type=checkbox bind:value={useLayout}>
            </label>
            <br>
            <!-- svelte-ignore a11y-missing-attribute -->
            <a class="waves-effect waves-light btn" style="margin-top: 5px;" on:click={drawLattice}>redraw lattice</a>
        </div>
    </div>
    <div style="position: absolute; bottom: 51px; right: 5px;">
        <h6>Network statistics:</h6>
        {#if selected.id == 1}
            <textarea readonly style="resize: none; width: 300px; height: 100px;">
Avg. #Species: {statistics.species_count}
Avg. #Reactions: {statistics.reaction_count}
            </textarea>
        {:else}
            <textarea readonly style="resize: none; width: 300px; height: 100px;">
#basic sets: {statistics.basic}
#closed reactive non basic sets: {statistics.closed}
#closed reactive sets only: {statistics.basic_closed}
#semi-self-maintaining sets: {statistics.ssm}
#organizations: {statistics.orgs}
#spurious union closure: {statistics.union}
#synergistic union closure: {statistics.syn_union}
#connected organizations: {statistics.conn_orgs}
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
                Creating Hasse
            </h5>
        </div>
    {/if}
	
</main>

<style>
    #hasse {
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
    div {
		--progress-bg: transparent;
		--track-bg: #ff9355;
	}
</style>
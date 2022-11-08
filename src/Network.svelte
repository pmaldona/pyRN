<script>
    //import { Network } from 'vis-network';
    //import { DataSet } from 'vis-data/peer';
    import { network } from './store/NetworkStore';
    import { filename } from './store/FileStore';
    import { protoSyn } from './store/ProtoSynStrStore';
    import Modal,{getModal} from './Modal.svelte'
    export let genNetwork;
    export let genProtoSyn;
    export let addSpecies;
    export let expNetwork;
    export let randInflow;
    export let randOutflow;

    let with_inflow = false;
	let new_extra = 1;
	let percentage_of_species = 0.1;
	let extra = 1;
	let percentage_of_reactions = 0.1;
	let label_of_species = "x";
	let extra_inflow = 0.1;
	let extra_outflow = 0.1;

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

    C2S.prototype.circle = CanvasRenderingContext2D.prototype.circle;
    C2S.prototype.square = CanvasRenderingContext2D.prototype.square;
    C2S.prototype.triangle = CanvasRenderingContext2D.prototype.triangle;
    C2S.prototype.triangleDown = CanvasRenderingContext2D.prototype.triangleDown;
    C2S.prototype.star = CanvasRenderingContext2D.prototype.star;
    C2S.prototype.diamond = CanvasRenderingContext2D.prototype.diamond;
    C2S.prototype.roundRect = CanvasRenderingContext2D.prototype.roundRect;
    C2S.prototype.ellipse_vis = CanvasRenderingContext2D.prototype.ellipse_vis;
    C2S.prototype.database = CanvasRenderingContext2D.prototype.database;
    C2S.prototype.arrowEndpoint = CanvasRenderingContext2D.prototype.arrowEndpoint;
    C2S.prototype.circleEndpoint = CanvasRenderingContext2D.prototype.circleEndpoint;
    C2S.prototype.dashedLine = CanvasRenderingContext2D.prototype.dashedLine;

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
                edge.color = {color: '#848484', highlight: '#848484', hover: '#848484', inherit: false, opacity: 1.0};
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
            console.log(network.get_network());
            console.log(nodes);
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
                edge.color = {color: '#848484', highlight: '#848484', hover: '#848484', inherit: false, opacity: 1.0};
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

    function action(e) {
		e.preventDefault();
        let labels = label_of_species.split(',');
        for(let i = 0; i < labels.length; i++) {
            labels[i] = labels[i].trim();
        }
        if(labels.length != new_extra) {
            labels = [labels[0]];
        }

		let form_obj = {
            Nse: new_extra,
	        p: new_extra==0 ? percentage_of_species : 0,
	        extra: extra,
	        m: extra==0 ? percentage_of_reactions : 0,
	        l: labels
	    };
        console.log(form_obj);
        getModal().close();
		let promise = addSpecies(form_obj);
        promise.then(result => {
            if(result == true) {
                drawNetwork();
            }
        });
	}

    function exportSvg()
    {
        var networkContainer = network.get_network().body.container;
        var ctx = new C2S({width: networkContainer.clientWidth, height: networkContainer.clientWidth, embedImages: true});

        var canvasProto = network.get_network().canvas.__proto__;
        var currentGetContext = canvasProto.getContext;
        canvasProto.getContext = function()
        {
            return ctx;
        }
        var svgOptions = {
            nodes: {
                shapeProperties: {
                    interpolation: false //so images are not scaled svg will get full image
                },
                scaling: { label: { drawThreshold : 0} },
                font:{color:'#000000'}
            },
            edges: {
                scaling: { label: { drawThreshold : 0} }
            }
        };
        network.get_network().setOptions(svgOptions);
        network.get_network().redraw();
        network.get_network().setOptions(options);
        canvasProto.getContext = currentGetContext;
        ctx.waitForComplete(function()
            {
                var svg = ctx.getSerializedSvg();
                showSvg(svg);
            });
    }
    function showSvg(svg)
    {
        var svgBlob = new Blob([svg], {type: 'image/svg+xml'});
        openBlob(svgBlob, "network.svg");
    }

    function openBlob(blob, fileName)
	  {
		if(window.navigator && window.navigator.msSaveOrOpenBlob)
        {

            //blobToDataURL(blob, function(dataurl){window.open(dataurl);});
            window.navigator.msSaveOrOpenBlob(blob,fileName);
        }
        else
        {
			var a = document.getElementById("blobLink");
			if(!a)
			{
				a = document.createElement("a");
				document.body.appendChild(a);
				a.setAttribute("id", "blobLink");
				a.style = "display: none";
			}
			var data = window.URL.createObjectURL(blob);
			a.href = data;
			a.download = fileName;
			a.click();
			setTimeout(function()
				{
				// For Firefox it is necessary to delay revoking the ObjectURL
				window.URL.revokeObjectURL(data);
				}
				, 100);
        }
    }

    function exportNetwork() {
        expNetwork();
    }

    function addInflow() {
        randInflow(extra_inflow);
        drawNetwork();
    }

    function addOutflow() {
        randOutflow(extra_outflow);
        drawNetwork();
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
            <div style="margin-top: 5px; height: 250px; overflow-y: auto;">
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn" style="margin-top: 5px;" on:click={drawNetwork}>redraw</a>
                <div>
                    {#if $filename != ""}
                        {#if selected.id == 1}
                            <!-- svelte-ignore a11y-missing-attribute -->
                            <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => getModal().open()}>Add Extra Species</a>
                            <div>
                                <input type="number" id="perc_species" min="2" max="1" bind:value={extra_inflow} style="width: 20%;">
                                <!-- svelte-ignore a11y-missing-attribute -->
                                <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => addInflow}>Add Inflow</a>
                            </div>
                            <div>
                                <input type="number" id="perc_species" min="0" max="1" bind:value={extra_outflow} style="width: 20%;">
                                <!-- svelte-ignore a11y-missing-attribute -->
                                <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => addOutflow}>Add Outflow</a>
                            </div>
                        {/if}
                    {/if}
                </div>
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => exportSvg()}>Export to SVG</a>
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => exportNetwork()}>Export Network</a>
            </div>
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

    <Modal>
        <h4>Add Species</h4>
        <div style="height: 250px; overflow:auto;">
            <form onsubmit="action(e)">
                <table>
                    <tr>
                        <td>Number of added species</td>
                        <td><input type="number" id="species" step="1" bind:value={new_extra} style="width: 20%;"></td> 
                    </tr>
                    {#if new_extra==0}
                        <tr>
                            <td>Percentage of species added to reactions</td>
                            <td>
                                <input type="number" id="perc_species" min="2" step="1" bind:value={percentage_of_species} style="width: 20%;">
                            </td>
                        </tr>
                    {/if}
                    <tr>
                        <td>Number of reactions to add the species to</td>
                        <td>
                            <input type="number" id="reactions" step="1" bind:value={extra} style="width: 20%;">
                        </td>
                        <!-- <td><input type="number" id="reactions" min="1" step="1" bind:value={random_reactions} style="width: 20%;"></td>  -->
                    </tr>
                    {#if extra==0}
                        <tr>
                            <td>Percentage of reactions to add the species to</td>
                            <td>
                                <input type="number" id="perc_reactions" bind:value={percentage_of_reactions} style="width: 20%;">
                            </td>
                        </tr>
                    {/if}
                    <tr>
                        <td>Label of the new species</td>
                        <td>
                            <input type="text" id="reactions" step="1" bind:value={label_of_species} style="width: 100%;">
                        </td>
                    </tr>
                    <tr>
                        <td>Species will be split by ','</td>
                    </tr>
                </table>
                <input type="submit" value="Generate" on:click={action}>
            </form>
        </div>
    </Modal>
</main>

<style>
    #network {
        width: 650px;
        height: 650px;
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

    table, tr, td {
		border: 0px solid black;
		margin: 5px;
		padding: 0px;
	}
</style>
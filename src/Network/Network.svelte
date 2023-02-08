<script>
    import { GraphType, setVisObject } from '../misc/drawNetwork';
    import exportSvg from '../misc/svg';
    import { genRNStr } from './Network';

    export let is_loaded;

    let network = {};

    let graphs = [
		{ id: 1, text: `Reaction Network` },
		{ id: 2, text: `Protosynergies` },
	];

	let selected = graphs[0];

    let has_file_open = false;
    
    let statistics = {
        species_count: 0,
        reaction_count: 0,
        generators: 0,
        transitions: 0
    };

    drawNetwork();

    function _genRNStr() {
        is_loaded().then(result => {
            has_file_open = result;
        });
        return genRNStr(GraphType.Network);
    }

    function genProSyn() {
        is_loaded().then(result => {
            has_file_open = result;
        });
        let promise = eel.gen_protosynergetic()();
		let result = promise.then(result => {
            if (result == null) {
                return;
            }
            for(let i = 0; i < result.edges.length; i++) {
                let edge = result.edges[i];
                //console.log(edge);
                if(edge.title){
                    edge.label = edge.title;
                }
                edge.color = {color: '#848484', highlight: '#848484', hover: '#848484', inherit: false, opacity: 1.0};
            }
            
            return setVisObject(GraphType.Network, result.nodes, result.edges, result.options, false);
        }); 
        return result;
    }

    async function drawNetwork() {
        if(selected.id == 2) {
            let res = await genProSyn();
            if(res != undefined){
                network = res;
            }
        } else {
            let res = await _genRNStr();
            if(res != undefined){
                network = res;
            }
        }
    }

    function exportNetwork() {
        let exported = window.electron.save().then(result => {
            let path = result;
            let success = eel.export_network(path)();
            let ret = success.then(result => {
				return result
            });
            return ret;
        });
		if(exported) {
			return true;
		}
		return false;
    }
    
</script>

<main>
    <div style="display: flex;">
        <div id="network"></div>
        <div style="margin: 15px;">
            <select bind:value={selected} on:change="{(e) => {drawNetwork();}}" style="display:block;">
                {#each graphs as graph}
                    <option value={graph}>
                        {graph.text}
                    </option>
                {/each}
            </select>
            <div style="margin-top: 5px; height: 650px; overflow-y: auto;">
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn" style="margin-top: 5px;" on:click={drawNetwork}>Redraw</a>
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => exportSvg(network)}>Export to SVG</a>
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => exportNetwork()}>Export Network</a>
            </div>
        </div>
    </div>
    
    <!-- <div style="position: absolute; bottom: 51px; right: 5px;">
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
       
    </div> -->
    {#if has_file_open == false}
        <h5 style="position: absolute; top: 55px; left: 10px;">
            Please open a network file.
        </h5>
    {:else if network.body == undefined}
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
</style>
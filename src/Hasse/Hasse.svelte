<script>
    import OrgModal from './OrgModal.svelte';
    import Slider from '@bulatdashiev/svelte-slider';
    import Modal,{getModal} from '../components/Modal.svelte'
    import { GraphType, setVisObject } from '../misc/drawNetwork';
    import exportSvg from '../misc/svg';

    export let is_loaded;
    
    let modal;
    let cur_org = "";

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

    let has_file_open = false;
    
    let network = {};

    drawLattice();

    function genOrgStr() {
        is_loaded().then(result => {
            has_file_open = result;
        });
        let promise = eel.calculate_orgs()();
        promise.then(result => {
            if(result == null) {
                return;
            }
            for(let i = 0; i < result.edges.length; i++) {
                let edge = result.edges[i];
                //console.log(edge);
                if(edge.title){
                    edge.label = edge.title;
                }
                
            }
            console.log(result.nodes);
            network = setVisObject(GraphType.Hasse, result.nodes, result.edges, result.options, false, {"improvedLayout": false});

            network.on('click', async function(properties) {
                var ids = properties.nodes;
                console.log(network);
                var clickedNode = network.body.nodes[ids]; //nodes.get(ids)[0];
                if (clickedNode) {
                    console.log('clicked nodes:', clickedNode);
                    cur_org = clickedNode.options;
                    await modal.init(cur_org);
                    //modal.setReactions(reactions);
                    getModal().open();
                }
            });
        });
    }

    function genSynStr() {
        is_loaded().then(result => {
            has_file_open = result;
        });
        let promise = eel.gen_synergetic()();
        promise.then(result => {
            let nodes = useLayout ? layout(result.nodes, x_space[0], y_space[0]) : result.nodes;
            // nodes = new vis.DataSet(n);
            if(result == null) {
                return;
            }
            for(let i = 0; i < result.edges.length; i++) {
                let edge = result.edges[i];
                //console.log(edge);
                if(edge.title){
                    edge.label = edge.title;
                }
                
            }
            network = setVisObject(GraphType.Hasse, nodes, result.edges, result.options, false, undefined);
        });
    }

    function drawLattice() {
        if(selected.id == 2) {
            genSynStr();
        } else {
            genOrgStr();
        }
    }

    function layout(nodes, x_space = 75, y_space = 75) {
        let length_dict = {}
        nodes.forEach(node => {
            let len = node.title.split(" ").length;
            if(length_dict[len]){
                length_dict[len]["length"] = length_dict[len]["length"] + 1;
            }
            else{
                length_dict[len] = {};
                length_dict[len]["c"] = 0;
                length_dict[len]["length"] = 1;
            }
        });

        console.log(length_dict);

        nodes.forEach(node => {
            console.log(node.title + " => " + node.title.split(" ").length);
            let len = node.title.split(" ").length;
            let count = length_dict[len]["length"];
            let _y = -y_space * (len-1);
            let _x = 0;
            if (count > 1){
                _x = -x_space*((count-1)/2) + length_dict[len]["c"]*x_space
                length_dict[len]["c"] = length_dict[len]["c"]+1;
            }
            
            node["y"] = _y;
            node["x"] = _x;
            node["fixed"] = {x: false, y: true}
            console.log(length_dict);
        });

        return nodes;
    }
    
</script>

<main>
    <div style="display: flex;">
        <div id="hasse"></div>
        <div style="margin: 15px; height: 650px">
            <select bind:value={selected} on:change="{() => drawLattice()}" style="display:block">
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
            <br>
            <!-- svelte-ignore a11y-missing-attribute -->
            <a class="waves-effect waves-light btn" style="margin-top: 20px;" on:click={() => exportSvg()}>Export to SVG</a>
        </div>
    </div>
    <!-- <div style="position: absolute; bottom: 51px; right: 5px;">
        <h6>Network statistics:</h6>
        {#if selected.id == 1}
            <textarea readonly style="resize: none; width: 300px; height: 100px;">
Avg. #Species: {statistics.species_count}
Avg. #Reactions: {statistics.reaction_count}
Selected Org: {cur_org}
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
       
    </div> -->
    {#if has_file_open == ""}
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
                Creating Hasse
            </h5>
        </div>
    {/if}
	<Modal>
        <OrgModal bind:this={modal}></OrgModal>
    </Modal>
</main>

<style>
    #hasse {
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
    div {
		--progress-bg: transparent;
		--track-bg: #ff9355;
	}
</style>
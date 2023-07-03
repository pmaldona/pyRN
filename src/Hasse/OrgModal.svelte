<script>
    import { GraphType, generateGraph, setVisObject } from '../misc/drawNetwork';
    import { genRNStr } from '../Network/Network';
    import Modal,{getModal} from '../components/Modal.svelte'
    import ImageModal from '../components/ImageModal.svelte';
    import { saveNodePositions, setNodePositions } from '../actions/Base';

    let modal;

    let network = {};
    let network_array_length = 0;
    let network_index = 0;
    let organization = {};
    let species_ids = [];
    let deplete_species = {};
    let can_perturb = {};
    let reaction_ids = ["id1"];
    let inflowReactions = [];
    let outflowReactions = [];
    let catalystSpecies = [];
    let image_source;
    let org_indx = 0;

    let overprod = false;

    export async function init(org) {
        network_array_length = 0;
        organization = org;
        
        console.log(organization);
        console.log(organization.title);
        //"['cows' 'dung' 'farmer' 'fertilizer' 'grass' 'infrastructure' 'milk'
 //'money' 'water' 'worms']"
        species_ids = organization.title.substring(1,organization.title.length-1).replaceAll("'", "").replaceAll("\n","").replaceAll("\r","").split(" ");
        species_ids.forEach(species => {
            deplete_species[species] = false;
            can_perturb[species] = false;
        });
        org_indx = await eel.get_selected_org(organization.title)();
        let pert_spec = await eel.get_perturbing_species(org_indx)();
        console.log(pert_spec);
        pert_spec.substring(1,pert_spec.length-1).replaceAll("'", "").replaceAll("\n","").replaceAll("\r","").split(" ").forEach(species => {
            can_perturb[species] = true;
        });
        console.log("Before Draw!");
        await drawNetworkFromIndex(0);
        console.log("After Draw!");
        let nodes = Object.values(network.body.nodes);
        let reaction_str = await eel.get_reactions()();
        let reactions = [];
        let included_reactions = [];
        console.log(species_ids);
        nodes.forEach(node => {
            if(!species_ids.includes(node.id)){
                let included = true
                if(node.options.shape=="square") {
                    node.edges.forEach(edge => {
                        let id = edge.toId;
                        if(node.id == edge.toId) {
                            id = edge.fromId;
                        }

                        if(species_ids.includes(id) == false) {
                            included = false;
                        }
                    });
                }else {
                    included = false;
                }
                if(included == false) {
                    node.options.color = {
                        background: "#E8E8E8",  
                        border: "#E8E8E8",
                        highlight: {
                            background: "#E8E8E8",
                            border: "#E8E8E8"
                        },
                        hover: {
                            background: "#E8E8E8",
                            border: "#E8E8E8"
                        }
                    }
                    node.title = "";
                    node.labelModule.elementOptions.label = "";
                    node.labelModule.lineCount = 0;
                    node.labelModule.lines = [];
                    if(node.options.shape=="square") {
                        
                        node.edges.forEach(edge => {
                            edge.options.color = {color: '#E8E8E8', highlight: '#E8E8E8', hover: '#E8E8E8', inherit: false, opacity: 1.0};;
                            edge.title = "";
                            edge.labelModule.elementOptions.label = "";
                            edge.labelModule.lineCount = 0;
                            edge.labelModule.lines = [];
                        });
                    }
                }
                else {
                    if(node.options.shape=="square") {
                        if(!reaction_ids.includes(node.id)) {
                            let reac = "";
                            reaction_str.forEach(str => {
                                if (str.substring(0,node.id.length)== node.id) {
                                    reac = str.substring(node.id.length+2);
                                }
                            });
                            included_reactions.push(node);
                            reactions.push({id: node.id, str: reac});
                        }
                    }
                }
            }
        });
        reaction_ids = reactions;
        catalystSpecies = getCatalists();
        getInOutflow(included_reactions);
        getBasicSets();
        console.log(inflowReactions);
        network.body.nodes = setNodePositions(network.body.nodes);
        network.redraw();
        saveNodePositions(network.body.nodes);
    }

    async function drawNetworkFromIndex(index) {
        let _network = await eel.get_network_graph(index)();
        console.log(_network);
        // _network = JSON.parse(_network);
        // _network = toLower(_network);
        for(let i = 0; i < _network.edges.length; i++) {
            let edge = _network.edges[i];
            if(edge.title){
                edge.label = edge.title;
            }
            edge.color = {color: '#848484', highlight: '#848484', hover: '#848484', inherit: false, opacity: 1.0};
            edge.smooth = {enabled: false};
            for(let j = 0; j < _network.edges.length; j++) {
                if(j != i) {
                    if(edge.from == _network.edges[j].to && edge.to == _network.edges[j].from) {
                        edge.smooth = {enabled: true, type: "discrete", roundness: 0.7};
                    }
                }
            }
            edge.physics = false;
        }
        console.log(_network.nodes);
        let res = setVisObject(GraphType.OrgNetwork, _network.nodes, _network.edges, _network.options, true);
        if(res != undefined){
            network = res;
            console.log("NETWORK");
            network.on("dragEnd", function(properties) {
                var ids = properties.nodes;
                if(ids.length == 0) {
                    return;
                }
                saveNodePositions(network.body.nodes);
            });
            network.body.nodes = setNodePositions(network.body.nodes);
            saveNodePositions(network.body.nodes);
            network.fit();
        }
    }

    function getCatalists() {
        let nodes = Object.values(network.body.nodes);
        let ids = reaction_ids.map(x => x.id);
        let catalysts = [];
        nodes.forEach(node => {
            if(ids.includes(node.id)){
                let from = [];
                let to = [];
                node.edges.forEach(edge => {
                    if(!from.includes(edge.fromId) && node.id != edge.fromId) {
                        from.push(edge.fromId)
                    }
                    if(!to.includes(edge.toId) && node.id != edge.toId) {
                        to.push(edge.toId);
                    }
                });
                from.forEach(f => {
                    if(to.includes(f)) {
                        let _node = {};
                        nodes.forEach(_n => {
                            if(_n.id == f) {
                                _node = _n;
                            }
                        });
                        
                        let infl = 0;
                        let outfl = 0;
                        _node.edges.forEach(edge => {
                            if(edge.toId == f) {
                                infl = Number.parseFloat(edge.title) ? Number.parseFloat(edge.title) : 0;
                            }
                            if(edge.fromId == f) {
                                outfl = Number.parseFloat(edge.title) ? Number.parseFloat(edge.title) : 0;
                            }
                        });
                        if(infl - outfl >= 0) {
                            if(!catalysts.includes(f)){
                                catalysts.push(f);
                            }
                        }
                    }
                });
            }
        });
        return catalysts;
    }

    function getInOutflow(included) {
        console.log(included)
        let outflow = [];
        let inflow = [];
        included.forEach(reaction => {
            if(reaction.edges.length == 1) {
                if(species_ids.includes(reaction.edges[0].fromId)) {
                    outflow.push(reaction.id);
                }
                if(species_ids.includes(reaction.edges[0].toId)) {
                    inflow.push(reaction.id);
                }
            }
        });
        inflowReactions = inflow;
        outflowReactions = outflow;
    }

    function getBasicSets() {

        let promise = eel.get_basics_from_set(species_ids)();
        promise.then(result => {
            if(result == null) {
                return;
            }
            image_source = result;
        });
    }

    async function onImgClick() {
        //modal.setReactions(reactions);
        getModal('inner').open();
        await modal.init(image_source);
    }

    function toggleOverrep() {
        console.log("Overrep!");
    }

    async function depletion() {
        let spec = [];
        Object.keys(deplete_species).forEach(key => {
            if(deplete_species[key]) {
                spec.push(key);
            }
        });
        let result = await eel.perturb(spec, org_indx)();
        if(result != 0){
            network_array_length = result;
        }
        // console.log(network_object);
        // for(let node in network_object["1"].nodes){
        //     console.log(network_object["1"].nodes[node].id);
        // }
        // console.log(network_object["1"].nodes)
        // for(let network_key in Object.keys(network_object)) {
        //     let nodes = new vis.DataSet(network_object[network_key].nodes);
        //     let edges = new vis.DataSet(network_object[network_key].edges);
        //     let options = JSON.parse(network_object[network_key].options);

        //     let net = generateGraph(GraphType.OrgNetwork, {nodes: nodes, edges: edges, options: options});
        //     console.log(net);
        //     net.body.nodes = setNodePositions(net.body.nodes);
        //     network_array.push(net);
        // }
    }

    async function changeNetwork() {
        saveNodePositions(network.body.nodes);
        await drawNetworkFromIndex(network_index);
        network.redraw();
        network.body.nodes = setNodePositions(network.body.nodes);
    }
</script>

<main>
    <h4>Organization: {organization.label}</h4>
	<div style="display: flex; width:850px;">
        <div>
            <div id="org_network"></div>
            <div id= "slider">
                {#if network_array_length > 0}
                    <label>
                        Network Image: <input type="range" bind:value={network_index} min="0" max="{network_array_length}" on:change={changeNetwork} />
                    </label>
                {/if}
            </div>
        </div>
        <div style="margin: 15px;">
            <div style="margin-top: 5px; height: 500px; width:250px; overflow-y: auto;">
                <label style={"color: black; font-size: 14px;"}>
                    {"Show overproduced: "}
                    <input type="checkbox" bind:checked={overprod} style="opacity: 1; position: relative;" on:change={toggleOverrep}>
                </label>
                <hr class="solid">
                <p>Includes Species:</p>
                <style type="text/css">
                    .tg  {border-collapse:collapse;border-spacing:0;}
                    .tr  {border-bottom:none;}
                    .tg td{border-color:black;border-style:none;overflow:hidden;padding:0px 0px;word-break:normal;}
                    .tg th{border-color:black;border-style:none;overflow:hidden;padding:0px 0px;word-break:normal;}
                    .tg .tg-0lax{text-align:left;vertical-align:center;}
                </style>
                <table class="tg">
                    <thead>
                        {#each species_ids as species}
                            <tr class = "tr">
                                <td class="tg-0lax">{species}</td>
                                <td class="tg-0lax">
                                    {#if can_perturb[species]}
                                        <!-- svelte-ignore a11y-missing-attribute -->
                                        <a class="waves-effect waves-light btn-flat" style="margin: 4px; heigth: 5px;" on:click={() => {deplete_species[species] = !deplete_species[species]}}>{deplete_species[species] ? "[x]" : "[ ]"}</a>
                                    {/if}
                                </td>
                            </tr>
                        {/each}
                        <tr class = "tr">
                            <td><p>Total: {species_ids.length}</p></td>
                            <td>
                                <!-- svelte-ignore a11y-missing-attribute -->
                                <a class="waves-effect waves-light btn" style="margin-top: 5px;" on:click={depletion}>Deplete marked</a>
                            </td>
                        </tr>
                    </thead>
                </table>
                
                <hr class="solid">
                <p>Included Reactions:</p>
                <ul>
                    {#each reaction_ids as {id, str}}
                        <li>
                            {id}: {str}
                        </li>
                    {/each}
                </ul>
                <p>Total: {reaction_ids.length}</p>
                <hr class="solid">
                <p>Number of Catalists: {catalystSpecies.length}</p>
                <hr class="solid">
                <p>Size of inflow: {inflowReactions.length}</p>
                <hr class="solid">
                <p>Size of outflow: {outflowReactions.length}</p>
                <hr class="solid">
                {#if image_source != undefined}
                    <!-- svelte-ignore a11y-missing-attribute -->
                    <img src='data:image/png;base64,{image_source}' style="max-width: 100%; height: auto; object-fit: contain;" on:click={onImgClick}>
                {/if}
                <hr class="solid">
                <p>Number of Synergies: {reaction_ids.length}</p>
            </div>
        </div>
    </div>
    <Modal id="inner">
        <ImageModal bind:this={modal}></ImageModal>
    </Modal>
</main>

<style>
    #org_network {
        width: 550px;
        height: 550px;
        border: 1px solid lightgray;
        margin: 2px;
    }
    .solid {
        border-top: 1px solid #bbb;
    }
</style>
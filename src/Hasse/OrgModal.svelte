<script>
    import { GraphType } from '../misc/drawNetwork';
    import { genRNStr } from '../Network/Network';
    import Modal,{getModal} from '../components/Modal.svelte'
    import ImageModal from '../components/ImageModal.svelte';

    let modal;

    let network = {};
    let organization = {};
    let species_ids = [];
    let reaction_ids = ["id1"];
    let inflowReactions = [];
    let outflowReactions = [];
    let catalystSpecies = [];
    let image_source;

    export async function init(org) {
        organization = org;
        console.log(organization);
        console.log(organization.title);
        //"['cows' 'dung' 'farmer' 'fertilizer' 'grass' 'infrastructure' 'milk'
 //'money' 'water' 'worms']"
        species_ids = organization.title.substring(1,organization.title.length-1).replaceAll("'", "").replaceAll("\n","").replaceAll("\r","").split(" ");
        network = await genRNStr(GraphType.OrgNetwork);
        console.log(network);
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
        network.redraw();
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
</script>

<main>
    <h4>Organization: {organization.label}</h4>
	<div style="display: flex; width:850px;">
        <div id="org_network"></div>
        <div style="margin: 15px;">
            <div style="margin-top: 5px; height: 500px; width:250px; overflow-y: auto;">
                <p>Includes Species:</p>
                <ul>
                    {#each species_ids as species}
                        <li>
                            {species}
                        </li>
                    {/each}
                </ul>
                <p>Total: {species_ids.length}</p>
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
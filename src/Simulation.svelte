<script>
    import { filename } from './store/FileStore';
    export let genBasicsSp;
    export let genBasicsR;
    export let genStoich;

    let plots = [
        { id: 1, text: `Stoichiometry` },
		{ id: 2, text: `Basics with all Species` },
		//{ id: 2, text: `Basics with all Reactions` }
	];

    let selected = plots[0];

    let image_source;

    if($filename != "" && image_source == undefined) {
          plot();
    }

    function pltBasicsSp() {
        genBasicsSp().then(result => {
            image_source = result;
        }); 
    }

    function pltBasicsR() {

    }

    function pltStoich() {
        genStoich().then(result => {
            image_source = result;
        }); 
    }

    function plot() {
        if(selected.id == 2) {
            pltBasicsSp();
        } else if (selected.id == 3){
            pltBasicsR();
        } else {
            pltStoich();
        }
    }
</script>

<main>
    <div style="display: flex;">
        <div id="plot">
            {#if image_source != undefined}
                <img src='data:image/png;base64,{image_source}' style="max-width: 100%; height: auto; object-fit: contain;">
            {/if}
        </div>
        <div style="margin: 15px;">
            <select bind:value={selected} on:change="{() => console.log(selected)}" style="display:block">
                {#each plots as plot}
                    <option value={plot}>
                        {plot.text}
                    </option>
                {/each}
            </select>
            <!-- svelte-ignore a11y-missing-attribute -->
            <a class="waves-effect waves-light btn" on:click={plot}>replot</a>
        </div>
    </div>
    {#if $filename == ""}
        <h5 style="position: absolute; top: 55px; left: 10px;">
            Please open a network file.
        </h5>
    {:else if image_source == undefined}
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
                Creating Plot
            </h5>
        </div>
    {/if}
	
</main>

<style>
    #plot {
        width: 450px;
        height: 450px;
        border: 1px solid lightgray;
        margin: 2px;
        display: flex;
        justify-content: center;
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
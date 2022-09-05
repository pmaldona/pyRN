<script>
    import { filename } from './store/FileStore';
    export let genBasicsSp;
    export let genBasicsR;
    export let genStoich;
    export let genConcentrations;
    export let genRates;

    let plots = [
        { id: 1, text: `Stoichiometry` },
        { id: 2, text: `Concentrations`},
        { id: 3, text: `Rates`},
		//{ id: 4, text: `Basics with all Species` },
		//{ id: 5, text: `Basics with all Reactions` }
	];

    let selected = plots[0];

    let image_source;

    let timeStart = 0;
    let timeFinal = 50;
    let steps = 100;
    let cutoff = 0.1;
    let i_sp = undefined;
    let rt = undefined;

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

    function pltConcentrations() {
        genConcentrations(timeStart, timeFinal, steps, cutoff).then(result => {
            image_source = result;
        });
    }

    function pltRates() {
        genRates(timeStart, timeFinal, steps, cutoff).then(result => {
            image_source = result;
        });
    }

    function plot() {
        if(selected.id == 2) {
            pltConcentrations();
        } else if (selected.id == 3){
            pltRates();
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
        <div style="margin: 15px; width: 250px;">
            <select bind:value={selected} on:change="{() => console.log(selected)}" style="display:block;">
                {#each plots as plot}
                    <option value={plot}>
                        {plot.text}
                    </option>
                {/each}
            </select>
            <label>
                Start:
                <input type=number bind:value={timeStart}>
            </label>
            <label>
                End: 
                <input type=number bind:value={timeFinal}>
            </label>
            <label>
                Steps: 
                <input type=number bind:value={steps}>
            </label>
            <label>
                Cutoff: 
                <input type=number bind:value={cutoff}>
            </label>
            <div>
                Initial Concentrations:
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn-flat" on:click={plot}>Init</a>
                <br>
                Rate Constants:
                <!-- svelte-ignore a11y-missing-attribute -->
                <a class="waves-effect waves-light btn-flat" on:click={plot}>Rate</a>
            </div>
            <!-- svelte-ignore a11y-missing-attribute -->
            <a class="waves-effect waves-light btn" style="margin-top: 5px;" on:click={plot}>replot</a>
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
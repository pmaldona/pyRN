<script>
    import { text } from 'svelte/internal';
    import { filename } from './store/FileStore';
    export let genBasicsSp;
    export let genBasicsR;
    export let genStoich;
    export let genConcentrations;
    export let genRates;
    export let randSimpleWalk;
    export let plotSimpleWalk;
    export let plotAbstraction;
    export let plotTrajectory;

    let plots = [
        { id: 1, text: `Stoichiometry`, disabled: false},
        { id: 2, text: `Concentrations`, disabled: true},
        { id: 3, text: `Rates`, disabled: true},
        { id: 4, text: `Simple Random Walk`, disabled: false},
        { id: 5, text: `Abstraction Size`, disabled: true},
        { id: 6, text: `Trajectory Hasse`, disabled: true},
		//{ id: 4, text: `Basics with all Species` },
		//{ id: 5, text: `Basics with all Reactions` }
	];
    let keys = [];

    let selected = plots[0];
    let selected_key = undefined;

    let image_source;

    let timeStart = 0;
    let timeFinal = 50;
    let steps = 100;
    let cutoff = 0.1;
    let i_sp = undefined;
    let rt = undefined;
    let w = 10;
    let l = 10;
    let d = 1;
    let nmin = 3;
    let fname = "rand_walk.json";
    let has_random_walk = false;
    let convPert = false;

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

    function start_random_walk() {
        randSimpleWalk(w, l, d, nmin, fname).then(result => {
            if (result != false) {
                console.log(result);
                keys = result;
                if(keys.length > 0) {
                    selected_key = keys[0];
                    plotSimpleWalk(selected_key).then(result => {
                        image_source = result;
                    });
                    has_random_walk = true;
                    plots[4].disabled = false;
                    plots[5].disabled = false;
                }
            }
        });
    }

    function plot_simple_random_walk() {
        plotSimpleWalk(selected_key).then(result => {
            image_source = result;
        });
    }

    function plot_abstraction() {
        plotAbstraction(selected_key).then(result => {
            image_source = result;
        });
    }

    function plot_trajectory() {
        plotTrajectory(selected_key, convPert).then(result => {
            image_source = result;
        });
    }

</script>

<main>
    <div style="display: flex;">
        <div id="plot">
            {#if image_source != undefined}
                <!-- svelte-ignore a11y-missing-attribute -->
                <img src='data:image/png;base64,{image_source}' style="max-width: 100%; height: auto; object-fit: contain;">
            {/if}
        </div>
        <div style="margin: 15px; width: 250px;">
            <select bind:value={selected} on:change="{() => console.log(selected)}" style="display:block;">
                {#each plots as plot}
                    {#if plot.disabled == true}
                        <option value={plot} disabled>
                            {plot.text}
                        </option>
                    {:else}
                        <option value={plot}>
                            {plot.text}
                        </option>
                    {/if}
                {/each}
            </select>
            {#if selected.id == 2 || selected.id == 3}
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
            {/if}
            {#if selected.id == 4}
                <label>
                    Number of walks:
                    <input type=number min="1" step="1" bind:value={w}>
                </label>
                <label>
                    Steps per walk: 
                    <input type=number min="1" step="1" bind:value={l}>
                </label>
                <label>
                    Change in active species: 
                    <input type=number min="0" step="1" bind:value={d}>
                </label>
                <label>
                    Minimal active species: 
                    <input type=number min="0" step="1" bind:value={nmin}>
                </label>
                <label>
                    Save walks as: 
                    <input type=text bind:value={fname}>
                </label>
                <div>
                    <!-- svelte-ignore a11y-missing-attribute -->
                    <a class="waves-effect waves-light btn-flat" on:click={start_random_walk}>Simple Random Walk</a>
                </div>
                {#if keys.length == 0}
                    <select value={0} disabled style="display:block;"></select>
                    <div>
                        <!-- svelte-ignore a11y-missing-attribute -->
                        <a class="waves-effect waves-light btn-flat" disabled>Plot Random Walk</a>
                    </div>
                {:else}
                    <select bind:value={selected_key} on:change="{() => console.log(selected_key)}" style="display:block;">
                        {#each keys as key}
                            <option value={key}>
                                {key}
                            </option>
                        {/each}
                    </select>
                    <!-- svelte-ignore a11y-missing-attribute -->
                    <div>
                        <a class="waves-effect waves-light btn-flat" on:click={plot_simple_random_walk}>Plot Random Walk</a>
                    </div>
                {/if}
            {/if}
            {#if selected.id == 5}
                <select bind:value={selected_key} on:change="{() => console.log(selected_key)}" style="display:block;">
                    {#each keys as key}
                        <option value={key}>
                            {key}
                        </option>
                    {/each}
                </select>
                <!-- svelte-ignore a11y-missing-attribute -->
                <div>
                    <a class="waves-effect waves-light btn-flat" on:click={plot_abstraction}>Plot Abstraction Size</a>
                </div>
            {/if}
            {#if selected.id == 6}
                <select bind:value={selected_key} on:change="{() => console.log(selected_key)}" style="display:block;">
                    {#each keys as key}
                        <option value={key}>
                            {key}
                        </option>
                    {/each}
                </select>
                <label>
                    Show Convergence and Peretubations: 
                    <input type="checkbox" bind:checked={convPert} style="opacity: 1; position: relative;">
                </label>
                <!-- svelte-ignore a11y-missing-attribute -->
                <div>
                    <a class="waves-effect waves-light btn-flat" on:click={plot_trajectory}>Plot Trajectory</a>
                </div>
            {/if}
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
        width: 650px;
        height: 650px;
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
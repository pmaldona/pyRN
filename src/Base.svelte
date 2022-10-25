<script>
    //import { ipcRenderer } from "electron";
	import Modal,{getModal} from './Modal.svelte'

	export let initialValues;
    export let setName;
	export let genRandNet;

	let with_inflow = false;
	let random_species = 2;
	let random_reactions = 2;
	let distributions = ["log scaled"];
	let distribution = distributions[0];
	let extra = 0.4;
	let pr = 0;
	let pp = 0;
	let inflow = 0.1;
	let outflow = 0.1;

    window.eel.expose(say_hello_js, 'say_hello_js');
	function say_hello_js(x) {
		console.log("Hello from " + x);
	}

    function openFile() {
        window.electron.open().then(result => {
            let _name = result;
            let success = eel.openFile(result)();
            success.then(result => {
                console.log(result);
                if(result == true) {
                    setName(_name);
                }
            });
            
        });
    }

	function genNetwork() {
		getModal().open();
		//genRandNet();
	}

	function action(e) {
		e.preventDefault();
		if(random_reactions == random_species) {
			let form_obj = {
				has_inflow: with_inflow,
				random_species: random_species,
				random_reactions: random_reactions,
				distribution: distribution,
				pr: pr,
				pp: pp,
				inflow: inflow,
				outflow: outflow
			};
			getModal().close();
			genRandNet(form_obj);
		}
	}
</script>

<main>
	<h1>CRNS UI</h1>
    <!-- svelte-ignore a11y-missing-attribute -->
    <a class="waves-effect waves-light btn" style="margin: 20px;" on:click={openFile}>Open File</a>
	<br>
	<!-- svelte-ignore a11y-missing-attribute -->
	<a class="waves-effect waves-light btn" style="margin: 20px;" on:click={genNetwork}>Generate Random Network</a>
</main>


<Modal>
	<h4>Initialise Random Network</h4>
	<div style="height: 300px; overflow:auto;">
		<form onsubmit="action(e)">
			<label style="color: black; font-size: 14px;">
				<input type="checkbox" bind:checked={with_inflow} style="opacity: 1; position: relative;">
				with inflow?
			</label>
			<table>
				<tr>
					<td>Number of Species</td>
					<td><input type="number" id="species" min="2" step="1" bind:value={random_species} style="width: 20%;"></td> 
				</tr>
				<tr>
					<td>Number of Reactions</td>
					<td>
						<input type="number" id="reactions" min="2" step="1" bind:value={random_reactions} style="width: 20%;">
					</td>
					<td>
						{#if random_reactions != random_species}
							<label for="" style="color: red; font-size: 8px;">
								Reactions and species must be equal.
							</label>
						{/if}
					</td>
				</tr>
				<tr>
					<td>Distribution</td>
					<td>
						<select bind:value={distribution} on:change="{() => console.log(distribution)}" style="display:block;">
							{#each distributions as dist}
								<option value={dist}>
									{dist}
								</option>
							{/each}
						</select>
					</td>
					<!-- <td><input type="number" id="reactions" min="1" step="1" bind:value={random_reactions} style="width: 20%;"></td>  -->
				</tr>
				{#if distribution = distributions[0]}
					<tr>
						<td>Scale penalties</td>
						<td><input type="number" id="reactions" bind:value={pr} style="width: 20%;"></td>
						<td><input type="number" id="reactions" bind:value={pp} style="width: 20%;"></td> 
					</tr>
				{/if}
				{#if with_inflow}
					<tr>
						<td>Inflow</td>
						<td><input type="number" id="inflow" min="0" max="1" bind:value={inflow} style="width: 20%;"></td> 
					</tr>
					<tr>
						<td>Outflow</td>
						<td><input type="number" id="outflow" min="0" max="1" bind:value={outflow} style="width: 20%;"></td> 
					</tr>
				{/if}
			</table>
			<input type="submit" value="Generate" on:click={action}>
		</form>
	</div>
</Modal>

<style>
	main {
		text-align: center;
		padding: 1em;
		max-width: 240px;
		margin: 0 auto;
	}

	h1 {
		color: #ff3e00;
		text-transform: uppercase;
		font-size: 4em;
		font-weight: 100;
	}

	@media (min-width: 640px) {
		main {
			max-width: none;
		}
	}

	table, tr, td {
		border: 0px solid black;
		margin: 5px;
		padding: 0px;
	}
</style>
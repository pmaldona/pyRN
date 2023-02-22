<script>
	import Button from '../components/Button.svelte';
	import { openFile, genNetwork, action } from '../actions/Base';
	import Modal from '../components/Modal.svelte'
    import { Form, Check, Table, Row, Cell, Number, Select } from '../components/Form/index';

	export let openPage;

	let random_species = 12;
	let random_vector = [];
	let inflow = 0;
	let outflow = 0;
	let transformation = 0;
	let synthesis = 0;
	let decompostion = 0;
	let single_replacement = 0;
	let double_replacement = 0;
	let advanced = false;

	const submit=(e) => {
		e.preventDefault();
		if(advanced) {
			random_vector= [inflow, outflow, transformation, synthesis, decompostion, single_replacement, double_replacement];
		} else {
			random_vector = undefined;
		}
		action( 
			random_species,
			random_vector
		).then(result => {
			if(result == true) {
				openPage(1);
			}
		});
	}
</script>

<main>
	<h1>CRNS UI</h1>
	<Button text={"Open File"} action={openFile.bind(null, openPage.bind(null, 1))}/>
	<br>
	<Button text={"Generate Random Network"} action={genNetwork} />
</main>


<Modal heading="Initialise Random Network" contentStyle="height: 300px; overflow:auto;">
	<form onsubmit="submit(e)">
		<label style={"color: black; font-size: 14px;"}>
			{"Advanced: "}
			<input type="checkbox" bind:checked={advanced} style="opacity: 1; position: relative;">
		</label>
		<Table>
			<Row>
				<Cell content={"Number of species"} />
				<Cell>
					<input type="number" min={2} max={undefined} step={1} bind:value={random_species} style={undefined}>
				</Cell>
			</Row>
			{#if advanced}
				<Row>
					<Cell content={"Number of inflow reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={inflow} style={undefined}>
					</Cell>
				</Row>
				<Row>
					<Cell content={"Number of outflow reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={outflow} style={undefined}>
					</Cell>
				</Row>
				<Row>
					<Cell content={"Number of transformation reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={transformation} style={undefined}>
					</Cell>
				</Row>
				<Row>
					<Cell content={"Number of synthesis reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={synthesis} style={undefined}>
					</Cell>
				</Row>
				<Row>
					<Cell content={"Number of decompostion reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={decompostion} style={undefined}>
					</Cell>
				</Row>
				<Row>
					<Cell content={"Number of single replacment reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={single_replacement} style={undefined}>
					</Cell>
				</Row>
				<Row>
					<Cell content={"Number of double replacment reactions"} />
					<Cell>
						<input type="number" min={0} max={undefined} step={1} bind:value={double_replacement} style={undefined}>
					</Cell>
				</Row>
			{/if}
		</Table>
		<input type="submit" value={"Generate"} on:click={submit}>
	</form>
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
</style>
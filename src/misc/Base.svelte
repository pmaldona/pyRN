<script>
	import Button from '../components/Button.svelte';
	import { openFile, genNetwork, action } from '../actions/Base';
	import Modal from '../components/Modal.svelte'
    import { Form, Check, Table, Row, Cell, Number, Select } from '../components/Form/index';

	export let openPage;

	let with_inflow = false;
	let random_species = 2;
	let random_reactions = 2;
	let distributions = ["x*0+1"];
	let distribution = distributions[0];
	let extra = 0.4;
	let pr = 0;
	let pp = 0;
	let inflow = 0.1;
	let outflow = 0.1;
</script>

<main>
	<h1>CRNS UI</h1>
	<Button text={"Open File"} action={openFile.bind(null, openPage.bind(null, 1))}/>
	<br>
	<Button text={"Generate Random Network"} action={genNetwork} />
</main>


<Modal heading="Initialise Random Network" contentStyle="height: 300px; overflow:auto;">
	<Form value={"Generate"} action={(e) => {
		e.preventDefault();
		action(
			with_inflow, 
			random_species, 
			random_reactions, 
			extra, 
			distribution, 
			pr, 
			pp, 
			inflow, 
			outflow, 
			openPage.bind(null, 1));
		}
	}>
		<Check labelText={"with inflow: "} labelStyle={"color: black; font-size: 14px;"} bind={with_inflow} />
		<Table>
			<Row>
				<Cell content={"Number of species"} />
				<Cell><Number bind={random_species} min="2" step="1"/></Cell>
			</Row>
			<Row>
				<Cell content={"Number of reactions"} />
				<Cell><Number bind={random_reactions} min="2" step="1"/></Cell>
			</Row>
			<Row>
				<Cell content={"Distribution"} />
				<Cell>
					<Select 
						bind={distribution} 
						callback={console.log.bind(null, distribution)} 
						elements={distributions}
					/></Cell>
			</Row>
			<Row>
				<Cell content={"Scale penalties"} />
				<Cell><Number bind={pr} /></Cell>
				<Cell><Number bind={pp} /></Cell>
			</Row>
			{#if with_inflow}
				<Row>
					<Cell content={"Inflow"} />
					<Cell><Number bind={inflow} min="0" max="1"/></Cell>
				</Row>
				<Row>
					<Cell content={"Outflow"} />
					<Cell><Number bind={outflow} min="0" max="1"/></Cell>
				</Row>
			{/if}
		</Table>
	</Form>
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
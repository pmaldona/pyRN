<script>
	import { DataSet } from 'vis-data/peer';
	import Base from './Base.svelte';
  	import Network from './Network.svelte';
	import Hasse from './Hasse.svelte';

  	const pages = [Base, Network, Hasse];

	let page = 0;

	let state = {};

	const setName = (name) => {
		state.name = name;
		state.hasse = {};
	}

	const openPage = (_p) => {
		page = _p;
		console.log(page);
	}

	const generateHasse = async () => {
		let promise = eel.calculate_orgs()();
		let hasse = await promise.then(result => {
			//console.log(result);
			return result;
		});
		state.hasse = {};
		state.hasse.nodes = new DataSet(hasse.nodes);
		state.hasse.edges = new DataSet(hasse.edges);
		state.hasse.options = new DataSet(hasse.options);
		return hasse;
	}

	const generateNetwork = async () => {
		console.log("genNetwork");
		let promise = eel.gen_network()();
		let network = await promise.then(result => {
			return result;
		});
		state.network_raw = network;
		state.network = {};
		state.network.nodes = new DataSet(network.nodes);
		state.network.edges = new DataSet(network.edges);
		//state.network.options = new DataSet(network.Options);

		return network;
	}
</script>

<nav>
	<div class="nav-wrapper">
		<ul id="nav-mobile" class="left">
			<li on:click={() => {openPage(0)}}><a>Main</a></li>
			<li on:click={() => {openPage(1)}}><a>Reaction Network</a></li>
			<li on:click={() => {openPage(2)}}><a>Hasse Diagramm</a></li>
			<li on:click={() => {openPage(3)}}><a>Simulation</a></li>
		</ul>
	</div>
</nav>

<svelte:component
  	this={pages[page]}
	{setName}
	genHasse={generateHasse}
	genNetwork={generateNetwork}
  	initialValues={state}
/>

<footer>
	<p>{state.name || "No file opened."}</p>
</footer>

<style>
	footer {
		position: absolute;
		bottom: 0px;
	}
</style>